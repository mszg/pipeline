import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from evaluation_policy import aggregate_assessments, evaluate_gene, load_manifest


def policy(criteria, sample_type="sample"):
    return {"sample_types": {sample_type: {"criteria": criteria}}}


class EvaluationPolicyTests(unittest.TestCase):
    def test_missing_zero_and_unconfigured_are_distinct(self):
        config = policy({"variants": {"min": 0}})
        self.assertEqual(evaluate_gene({"variants": 0}, config)["status"], "PASS")
        self.assertEqual(evaluate_gene({}, config)["status"], "NOT ASSESSED")
        result = evaluate_gene({"variants": 0}, policy({"variants": {}}))
        self.assertEqual(result["status"], "NOT ASSESSED")
        self.assertIn("no evaluation threshold", " ".join(result["reasons"]))

    def test_absent_policy_does_not_invent_thresholds(self):
        self.assertEqual(evaluate_gene({"variants": 0}, {})["status"], "NOT ASSESSED")

    def test_null_placeholder_is_unconfigured_not_zero(self):
        self.assertEqual(evaluate_gene({"depth": 0}, policy({"depth": {"min": None}}))["status"],
                         "NOT ASSESSED")

    def test_collector_can_supply_explicit_applicability(self):
        result = evaluate_gene({"depth": 50, "_not_applicable": ["phase_blocks"]},
                               policy({"depth": {"min": 20}, "phase_blocks": {"max": 1}}))
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["criteria"][1]["status"], "NOT APPLICABLE")

    def test_no_hets_is_not_applicable_not_failed(self):
        result = evaluate_gene({"heterozygous_variants": 0, "depth": 50}, policy({
            "phased_fraction": {"min": 1, "not_applicable_if": {"heterozygous_variants": 0}},
            "depth": {"min": 50},
        }))
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["criteria"][0]["status"], "NOT APPLICABLE")

    def test_missing_het_count_cannot_establish_not_applicable(self):
        result = evaluate_gene({}, policy({
            "phased_fraction": {"min": 1, "not_applicable_if": {"heterozygous_variants": 0}},
        }))
        self.assertEqual(result["status"], "NOT ASSESSED")

    def test_fail_missing_review_precedence_preserves_reasons(self):
        values = {"depth": 0, "imbalance": 0.8}
        rules = policy({"depth": {"min": 10}, "mapped": {"min": 1},
                        "imbalance": {"max": 0.7, "severity": "REVIEW"}})
        result = evaluate_gene(values, rules)
        self.assertEqual(result["status"], "FAIL")
        self.assertEqual(len(result["reasons"]), 3)
        del rules["sample_types"]["sample"]["criteria"]["depth"]
        self.assertEqual(evaluate_gene(values, rules)["status"], "NOT ASSESSED")

    def test_abo_pass_cannot_hide_kel_missing(self):
        abo = evaluate_gene({"depth": 50}, policy({"depth": {"min": 20}}), gene="ABO")
        kel = evaluate_gene({}, policy({"depth": {"min": 20}}), gene="KEL")
        result = aggregate_assessments([abo, kel])
        self.assertEqual(result["status"], "NOT ASSESSED")
        self.assertEqual(len(result["reasons"]), 2)

    def test_controls_do_not_inherit_sample_coverage(self):
        rules = policy({"depth": {"min": 20}})
        self.assertEqual(evaluate_gene({"depth": 0}, rules, "negative_control")["status"], "NOT ASSESSED")
        rules["sample_types"]["negative_control"] = {"criteria": {"retained_reads": {"max": 0}}}
        result = evaluate_gene({"depth": 0, "retained_reads": 0}, rules, "negative_control")
        self.assertEqual(result["status"], "PASS")
        self.assertEqual([x["metric"] for x in result["criteria"]], ["retained_reads"])

    def test_gene_threshold_overrides_are_local(self):
        rules = policy({"depth": {"min": 20}})
        rules["genes"] = {"ABO": {"sample_types": {"sample": {"criteria": {"depth": {"min": 50}}}}}}
        self.assertEqual(evaluate_gene({"depth": 30}, rules, gene="ABO")["status"], "FAIL")
        self.assertEqual(evaluate_gene({"depth": 30}, rules, gene="KEL")["status"], "PASS")

    def test_parse_errors_block_pass_and_are_visible(self):
        result = evaluate_gene({"depth": 50}, policy({"depth": {"min": 20}}),
                               evidence_issues=["VCF parsing error: truncated file"])
        self.assertEqual(result["status"], "NOT ASSESSED")
        self.assertIn("truncated file", " ".join(result["reasons"]))

    def test_bad_threshold_not_silently_ignored(self):
        with self.assertRaisesRegex(ValueError, "Unknown criterion"):
            evaluate_gene({"depth": 50}, policy({"depth": {"minimum": 20}}))

    def test_optional_missing_does_not_hide_required_pass(self):
        result = evaluate_gene({"depth": 20}, policy({"depth": {"min": 20},
            "optional_metric": {"min": 1, "required": False}}))
        self.assertEqual(result["status"], "PASS")
        self.assertIn("optional evidence", " ".join(result["reasons"]))


class ManifestIdentityTests(unittest.TestCase):
    def read(self, rows):
        columns = ["sample", "gene", "amplicon", "fastq_input", "reference", "sample_id",
                   "barcode", "run_id", "sample_type"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "samples.tsv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
                writer.writeheader()
                for row in rows:
                    writer.writerow({"gene": "ABO", "amplicon": "frag1", "fastq_input": "reads/",
                                     "reference": "ref.fa", **row})
            return load_manifest(path)

    def test_leading_zeros_external_strings_and_gene_grouping(self):
        groups = self.read([
            {"sample": "0001", "sample_id": "=SUM(1,1)", "barcode": "007", "gene": "ABO"},
            {"sample": "0001", "sample_id": "=SUM(1,1)", "barcode": "007", "gene": "KEL"},
        ])
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["sample_id"], "=SUM(1,1)")
        self.assertEqual(groups[0]["barcode"], "007")
        self.assertEqual(groups[0]["genes"], ["ABO", "KEL"])
        self.assertEqual(groups[0]["rows"][0]["analysis"], "0001__ABO")

    def test_repeated_control_labels_never_merge_different_barcodes_or_runs(self):
        groups = self.read([
            {"sample": "p1", "sample_id": "Positive", "barcode": "01", "run_id": "run1"},
            {"sample": "p2", "sample_id": "Positive", "barcode": "02", "run_id": "run1"},
            {"sample": "p3", "sample_id": "Positive", "barcode": "01", "run_id": "run2"},
        ])
        self.assertEqual(len(groups), 3)

    def test_pooled_barcodes_are_flagged_not_misattributed(self):
        groups = self.read([
            {"sample": "same", "barcode": "01", "amplicon": "frag1"},
            {"sample": "same", "barcode": "02", "amplicon": "frag2"},
        ])
        self.assertEqual(len(groups), 2)
        self.assertTrue(all(group["ambiguous_analyses"] == ["same__ABO"] for group in groups))
        self.assertTrue(all(group["rows"][0]["analysis_identity_ambiguous"] for group in groups))

    def test_missing_barcode_stays_explicitly_unprovided(self):
        group = self.read([{"sample": "0001"}])[0]
        self.assertEqual(group["barcode"], "")
        self.assertIn("Barcode not provided", " ".join(group["identity_notes"]))

    def test_duplicate_units_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate workflow unit"):
            self.read([{"sample": "Positive", "barcode": "01"},
                       {"sample": "Positive", "barcode": "01"}])

    def test_reused_legacy_units_keep_all_identities_for_partial_report(self):
        groups = self.read([{"sample": "Positive", "barcode": "01"},
                            {"sample": "Positive", "barcode": "02"}])
        self.assertEqual(len(groups), 2)
        self.assertTrue(all(group["ambiguous_analyses"] == ["Positive__ABO"] for group in groups))


if __name__ == "__main__":
    unittest.main()
