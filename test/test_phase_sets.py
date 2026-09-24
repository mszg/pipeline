import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from check_phase_sets import check_phase_sets


class PhaseSetTests(unittest.TestCase):
    def check(self, genotypes, targets="chr1\t0\t100\n", sample="TEST"):
        with tempfile.TemporaryDirectory() as directory:
            vcf, bed = Path(directory) / "phase.vcf", Path(directory) / "target.bed"
            header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
            body = "".join(f"{chrom}\t{position}\t.\tA\tC\t30\tPASS\t.\tGT:PS\t{gt}:{ps}\n" for chrom, position, gt, ps in genotypes)
            vcf.write_text(header + body)
            bed.write_text(targets)
            return check_phase_sets(vcf, bed, "TEST")

    def test_single_connected_set_and_unphased_het(self):
        rows = self.check([("chr1", 10, "0|1", "10"), ("chr1", 50, "1|0", "10"), ("chr1", 80, "0/1", ".")])
        self.assertEqual((rows[0]["phase_set"], rows[0]["phased_heterozygotes"], rows[0]["unphased_heterozygotes"]), ("10", 2, 1))

    def test_multiple_phase_sets_fail(self):
        with self.assertRaisesRegex(ValueError, "found 2"):
            self.check([("chr1", 10, "0|1", "10"), ("chr1", 50, "1|0", "50")])

    def test_missing_ps_fails(self):
        with self.assertRaisesRegex(ValueError, "has no PS"):
            self.check([("chr1", 10, "0|1", ".")])

    def test_no_informative_phasing_fails(self):
        for records in [[], [("chr1", 10, "1/1", ".")], [("chr1", 10, "0/1", ".")]]:
            with self.subTest(records=records), self.assertRaisesRegex(ValueError, "found 0"):
                self.check(records)

    def test_different_contigs_may_have_different_phase_sets(self):
        rows = self.check([("chr1", 10, "0|1", "10"), ("chr2", 50, "1|0", "50")], "chr1\t0\t100\nchr2\t0\t100\n")
        self.assertEqual(len(rows), 2)

    def test_sample_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError, "Expected one VCF sample"):
            self.check([("chr1", 10, "0|1", "10")], sample="OTHER")


if __name__ == "__main__":
    unittest.main()
