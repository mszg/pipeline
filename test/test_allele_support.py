"""Regression cases exercise actual BAM CIGARs and exact allele sequences."""
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from allele_support import aligned_window, classify_indel, count_indel_support, count_snv_support, indel_window, observations


def read(sequence, cigar, start=0, flag=0, name="read"):
    record = pysam.AlignedSegment()
    record.query_name = name
    record.query_sequence = sequence
    record.reference_id = 0
    record.reference_start = start
    record.mapping_quality = 60
    record.flag = flag
    record.cigarstring = cigar
    record.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    return record


class IndelSupportTests(unittest.TestCase):
    def setUp(self):
        self.reference = "GGACGTCC"
        self.deletion = indel_window(self.reference, 2, "AC", "A")
        self.insertion = indel_window(self.reference, 2, "A", "AT")

    def test_reference_requires_exact_deleted_bases(self):
        self.assertEqual(classify_indel(read(self.reference, "8M"), self.deletion), "REF")
        self.assertEqual(classify_indel(read("GGATGTCC", "8M"), self.deletion), "OTHER")

    def test_exact_deletion(self):
        self.assertEqual(classify_indel(read("GGAGTCC", "3M1D4M"), self.deletion), "ALT")

    def test_exact_insertion(self):
        self.assertEqual(classify_indel(read("GGATCGTCC", "3M1I5M"), self.insertion), "ALT")

    def test_wrong_anchor_never_exact_reference_or_alternate(self):
        for record, window in [
            (read("GGTCGTCC", "8M"), self.deletion),
            (read("GGTGTCC", "3M1D4M"), self.deletion),
            (read("GGTTCGTCC", "3M1I5M"), self.insertion),
        ]:
            with self.subTest(cigar=record.cigarstring):
                self.assertEqual(classify_indel(record, window), "OTHER")

    def test_missing_anchor_or_reference_skip(self):
        for record in [read("GGCGTCC", "2M1D5M"), read("GGAGTCC", "3M1N4M")]:
            self.assertEqual(classify_indel(record, self.deletion), "OTHER")

    def test_clipped_or_incomplete_windows(self):
        for record in [read("GGA", "3M"), read("GGAC", "4M"), read("GGATTT", "3M3S"), read("CGTCC", "5M", start=3)]:
            self.assertEqual(classify_indel(record, self.deletion), "OTHER")

    def test_insertion_needs_observed_right_flank(self):
        self.assertEqual(classify_indel(read("GGAT", "3M1I"), self.insertion), "OTHER")

    def test_ambiguous_base_and_extra_event(self):
        for record in [read("GGNCGTCC", "8M"), read("GGANGTCC", "8M"), read("GGATCGTCC", "3M1I5M"), read("GGATGTCC", "3M1D1I4M")]:
            self.assertEqual(classify_indel(record, self.deletion), "OTHER")

    def test_reverse_strand_uses_bam_orientation(self):
        self.assertEqual(classify_indel(read("GGAGTCC", "3M1D4M", flag=16), self.deletion), "ALT")

    def test_soft_clip_outside_window_is_allowed(self):
        self.assertEqual(classify_indel(read("TTGGAGTCCAA", "2S3M1D4M2S"), self.deletion), "ALT")

    def test_equal_and_mismatch_cigar_ops(self):
        self.assertEqual(classify_indel(read(self.reference, "3=1=4="), self.deletion), "REF")
        self.assertEqual(classify_indel(read("GGATGTCC", "3=1X4="), self.deletion), "OTHER")

    def test_repeat_shifted_deletion_is_alt_not_reference(self):
        window = indel_window("GATATATC", 0, "GAT", "G")
        self.assertEqual((window.ref, window.alt), ("GATATATC", "GATATC"))
        for cigar in ["1M2D5M", "3M2D3M", "5M2D1M"]:
            self.assertEqual(classify_indel(read("GATATC", cigar), window), "ALT")
        self.assertEqual(classify_indel(read("GATATATC", "8M"), window), "REF")
        self.assertEqual(classify_indel(read("GATAT", "5M"), window), "OTHER")

    def test_repeat_shifted_insertion_is_alt(self):
        window = indel_window("GATATC", 0, "G", "GAT")
        for cigar in ["1M2I5M", "3M2I3M", "5M2I1M"]:
            self.assertEqual(classify_indel(read("GATATATC", cigar), window), "ALT")

    def test_homopolymer_repeat_and_wrong_repeat_length(self):
        window = indel_window("CAAAAG", 0, "CA", "C")
        for cigar in ["1M1D4M", "4M1D1M"]:
            self.assertEqual(classify_indel(read("CAAAG", cigar), window), "ALT")
        self.assertEqual(classify_indel(read("CAAG", "1M2D3M"), window), "OTHER")

    def test_long_deletion_and_insertion(self):
        reference = "GA" + "C" * 12 + "T"
        window = indel_window(reference, 1, "A" + "C" * 10, "A")
        self.assertEqual(classify_indel(read("GACCT", "2M10D3M"), window), "ALT")
        window = indel_window("GACT", 1, "A", "A" + "T" * 12)
        self.assertEqual(classify_indel(read("GA" + "T" * 12 + "CT", "2M12I2M"), window), "ALT")

    def test_unsupported_or_no_reference_flank(self):
        self.assertIsNone(indel_window("AC", 0, "AC", "A"))
        self.assertIsNone(indel_window("ACN", 0, "AC", "A"))
        self.assertIsNone(indel_window("ACGT", 0, "AC", "AT"))
        with self.assertRaises(ValueError):
            indel_window("ACGT", 0, "TC", "T")

    def test_real_indexed_bam_filters_and_counts(self):
        records = [read(self.reference, "8M", name="ref"), read("GGAGTCC", "3M1D4M", name="alt"), read("GGTGTCC", "3M1D4M", name="bad")]
        records += [read("GGAGTCC", "3M1D4M", flag=flag, name=f"excluded{flag}") for flag in (256, 512, 1024, 2048)]
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "reads.bam")
            with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "ref", "LN": 8}]}) as bam:
                for record in records:
                    bam.write(record)
            pysam.index(path)
            with pysam.AlignmentFile(path) as bam:
                self.assertEqual(count_indel_support(bam, "ref", self.deletion), (1, 1, 1, True))


class PileupTests(unittest.TestCase):
    def test_snv_forward_reverse_and_placeholders(self):
        self.assertEqual(count_snv_support("A", "C", ".,AaCc*#<>NnGg"), (4, 2, 8, True))

    def test_adjacent_indel_is_not_simple_snv_support(self):
        self.assertEqual(count_snv_support("A", "C", ".+1tC-1t"), (0, 0, 2, True))

    def test_start_end_and_multi_digit_indels(self):
        parsed = observations("^].$c+12acgtacgtacgt,$")
        self.assertEqual(len(parsed), 3)
        self.assertEqual(parsed[1]["indels"], [("+", 12, "ACGTACGTACGT")])

    def test_malformed_pileup_fails(self):
        for bases in ["^", ".+t", ".+2a", ".+0", ".!", ".+1?"]:
            with self.subTest(bases=bases), self.assertRaises(ValueError):
                observations(bases)

    def test_anchor_only_indel_support_is_unsupported(self):
        self.assertEqual(count_snv_support("AC", "A", "*#<>"), (0, 0, 4, False))
        self.assertEqual(count_snv_support("A", "AT", "G+1t"), (0, 0, 1, False))


if __name__ == "__main__":
    unittest.main()
