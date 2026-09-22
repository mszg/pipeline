import csv
import runpy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace as NS

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/make_consensus_mask.py"


class ConsensusMaskTests(unittest.TestCase):
    def mask(self, hp=2, support=50, status="ACCEPT", target="chr1\t0\t20\n", uncertain="",
             low="chr1\t11\t12\n", ref="AC"):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contents = {"ref.fai": "chr1\t30\t0\t30\t31\n", "target.bed": target,
                        "uncertain.bed": uncertain, "low.bed": low}
            for filename, content in contents.items():
                (root / filename).write_text(content)
            row = dict(CHROM="chr1", POS="11", REF=ref, ALT="A", TYPE="DEL", STATUS=status,
                       PHASED_GT="0|1", HP1_ALT="0", HP2_ALT=str(support))
            with (root / "support.tsv").open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
                writer.writeheader(); writer.writerow(row)
            context = NS(input=NS(fai=root / "ref.fai", target_bed=root / "target.bed",
                                  uncertain_bed=root / "uncertain.bed", low_depth_bed=root / "low.bed", support_tsv=root / "support.tsv"),
                         output=NS(bed=root / "mask.bed"), params=NS(haplotype=hp, callable_min_depth=50))
            runpy.run_path(str(SCRIPT), init_globals={"snakemake": context})
            return {p for line in (root / "mask.bed").read_text().splitlines() for _, start, end in [line.split("\t")] for p in range(int(start), int(end))}

    def test_accepted_alt_deletion_rescued_at_depth_threshold(self):
        mask = self.mask()
        self.assertNotIn(10, mask)  # callable anchor allows the deletion to be applied
        self.assertNotIn(11, mask)  # deleted base is supported by 50 exact ALT reads

    def test_low_support_stays_masked(self):
        self.assertIn(11, self.mask(support=49))

    def test_reference_haplotype_stays_masked(self):
        self.assertIn(11, self.mask(hp=1))

    def test_unresolved_stays_masked(self):
        self.assertIn(11, self.mask(status="UNRESOLVED"))

    def test_uncertainty_mask_cannot_be_rescued(self):
        self.assertIn(11, self.mask(uncertain="chr1\t11\t12\n"))

    def test_outside_target_cannot_be_rescued(self):
        self.assertIn(11, self.mask(target="chr1\t0\t11\n"))

    def test_uncertain_anchor_prevents_deleted_base_rescue(self):
        mask = self.mask(uncertain="chr1\t10\t11\n")
        self.assertIn(10, mask)
        self.assertIn(11, mask)

    def test_outside_target_anchor_prevents_deleted_base_rescue(self):
        mask = self.mask(target="chr1\t11\t20\n")
        self.assertIn(10, mask)
        self.assertIn(11, mask)

    def test_low_depth_anchor_prevents_deleted_base_rescue(self):
        mask = self.mask(low="chr1\t10\t12\n")
        self.assertIn(10, mask)
        self.assertIn(11, mask)

    def test_partial_uncertainty_blocks_entire_deletion_rescue(self):
        mask = self.mask(ref="ACG", low="chr1\t11\t13\n", uncertain="chr1\t12\t13\n")
        self.assertNotIn(10, mask)
        self.assertIn(11, mask)  # other deleted base stays masked when variant is skipped
        self.assertIn(12, mask)


if __name__ == "__main__":
    unittest.main()
