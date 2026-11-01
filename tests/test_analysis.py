"""Scientific and execution invariants for the CPPTRAJ job generator."""
import copy
import unittest

from mdwb.analysis import (build_batches, build_pooled_batches,
                           enabled_analyses, quote_mask, safe_name)


SELECTIONS = {
    "protein": ":1-20", "nucleic": ":21-40", "solute": ":1-40",
    "backbone": ":1-20@N,CA,C,O",
    "representative": "(:1-20@CA)|(:21-40@C1')",
}


def config(level="standard"):
    return {"sections": [{"name": "protein", "mask": "protein"},
                         {"name": "DNA backbone", "mask": ":21-40@P,O5',C5',C4',C3',O3'"}],
            "analysis": {"level": level},
            "interactions": [{"name": "protein DNA", "mask1": "protein", "mask2": "nucleic"}]}


class AnalysisTests(unittest.TestCase):
    def test_section_local_fits_do_not_change_later_observables(self):
        batches = build_batches(config("basic"), SELECTIONS)
        commands = batches[0]["commands"]
        rms = [command for command in commands if command.startswith("rms ")]
        self.assertEqual(len(rms), 4)
        self.assertTrue(all(" ref [COMMON] " in command for command in rms))
        self.assertTrue(all(" nomod " in command or " nofit " in command for command in rms))
        self.assertFalse(any(command.startswith(("strip ", "autoimage")) for command in commands))

    def test_pca_and_clusters_use_a_shared_job_and_preserve_frame_order(self):
        cfg = config("advanced")
        per_rep = build_batches(cfg, SELECTIONS)
        self.assertFalse(any("pca" in batch["name"] or "cluster" in batch["name"] for batch in per_rep))
        pooled = build_pooled_batches(cfg, SELECTIONS)
        self.assertEqual({batch["name"] for batch in pooled}, {"pca_pooled", "cluster_pooled"})
        for batch in pooled:
            self.assertEqual(batch["commands"][1], "createcrd A_POOLED")
            self.assertTrue(batch["commands"][0].startswith('strip "!('))
            self.assertEqual(batch["coordinate_mask"], SELECTIONS["representative"])
            self.assertFalse(any("filter" in cmd for cmd in batch["commands"]))
            self.assertTrue(all(a["metadata"]["pooled"] for a in batch["artifacts"]))
        pca = next(b for b in pooled if b["name"] == "pca_pooled")
        post = pca["post_commands"]
        covariance = next(i for i, c in enumerate(post) if "matrix covar" in c)
        diagonalize = next(i for i, c in enumerate(post) if "diagmatrix" in c)
        project = next(i for i, c in enumerate(post) if "projection" in c)
        self.assertLess(covariance, diagonalize)
        self.assertLess(diagonalize, project)

    def test_contact_sums_are_not_claimed_as_probabilities(self):
        artifacts = [art for b in build_batches(config(), SELECTIONS) for art in b["artifacts"]]
        maps = [a for a in artifacts if a["analysis"].startswith("contacts_") and a["kind"] == "matrix"]
        self.assertEqual(len(maps), 2)
        for art in maps:
            self.assertEqual(art["unit"], "mean atom-pair contacts")
            self.assertIn("may exceed one", art["metadata"]["definition"])
            self.assertNotIn("limits", art["metadata"])

    def test_nucleic_apostrophes_survive_quoting(self):
        mask = "(:21-40)&@C1',O4'"
        self.assertEqual(quote_mask(mask), '"' + mask + '"')
        batches = build_batches(config(), SELECTIONS)
        self.assertTrue(any("C5'" in command for b in batches for command in b["commands"]))

    def test_masks_cannot_inject_cpptraj_commands(self):
        for mask in (':1\nquit', ':1" out overwrite.dat', ':1\\', ':1# hidden', '$UNTRUSTED', '\x00'):
            with self.subTest(mask=mask), self.assertRaises(ValueError):
                quote_mask(mask)

    def test_absent_molecule_classes_skip_specialist_actions(self):
        cfg = {"sections": [{"name": "ligand", "mask": ":1"}], "analysis": {"level": "standard"}}
        batches = build_batches(cfg, {"solute": ":1", "representative": ":1&!@H="})
        self.assertFalse({"dssp", "nastruct"} & {b["name"] for b in batches})
        self.assertIn("structure", {b["name"] for b in batches})

    def test_bidirectional_interface_hydrogen_bonds(self):
        batch = next(b for b in build_batches(config(), SELECTIONS) if b["name"].startswith("interface_"))
        hb = [c for c in batch["commands"] if c.startswith("hbond ")]
        self.assertEqual(len(hb), 2)
        self.assertIn('donormask "(:1-20)&@/N,O,F"', hb[0])
        self.assertIn('donormask "(:21-40)&@/N,O,F"', hb[1])

    def test_explicit_analysis_list_overrides_level(self):
        cfg = config("advanced")
        cfg["analysis"]["enabled"] = ["rg"]
        cfg["advanced"] = {"pairwise_rmsd": True}
        self.assertEqual(enabled_analyses(cfg), {"rg"})
        self.assertEqual(build_pooled_batches(cfg, SELECTIONS), [])
        self.assertTrue(all(c.startswith("radgyr ") for b in build_batches(cfg, SELECTIONS) for c in b["commands"]))

    def test_unknown_analyses_and_selections_fail_before_execution(self):
        cfg = config()
        cfg["analysis"]["enabled"] = ["ramachandrn"]
        with self.assertRaisesRegex(ValueError, "Unknown analyses"):
            build_batches(cfg, SELECTIONS)
        cfg = config()
        cfg["sections"][0]["mask"] = "protien"
        with self.assertRaisesRegex(ValueError, "Unknown selection"):
            build_batches(cfg, SELECTIONS)

    def test_output_names_cannot_collide_after_sanitization(self):
        self.assertNotEqual(safe_name("a/b"), safe_name("a_b"))
        self.assertNotEqual(safe_name("a b"), safe_name("a/b"))
        cfg = config("advanced")
        artifacts = [art for b in build_batches(cfg, SELECTIONS) + build_pooled_batches(cfg, SELECTIONS) for art in b["artifacts"]]
        self.assertEqual(len(artifacts), len({art["path"] for art in artifacts}))

    def test_invalid_resource_parameters_fail(self):
        for field, value in (("cluster_sieve", 0), ("cluster_count", 3.1), ("pca_modes", float("nan"))):
            cfg = config("advanced")
            cfg["advanced"] = {field: value}
            with self.subTest(field=field), self.assertRaises(ValueError):
                build_pooled_batches(cfg, SELECTIONS)

    def test_input_configuration_is_not_mutated(self):
        cfg = config("advanced")
        expected = copy.deepcopy(cfg)
        build_batches(cfg, SELECTIONS)
        build_pooled_batches(cfg, SELECTIONS)
        self.assertEqual(cfg, expected)


if __name__ == "__main__":
    unittest.main()
