"""Contracts for scientific input grouping, masks and explicit time metadata.

The synthetic Amber snippets test metadata parsing only: they intentionally omit
force-field sections and must never be presented as runnable MD topologies.
"""
from __future__ import annotations

import copy
import gzip
import json
import os
import tempfile
import unittest
from pathlib import Path

from mdwb.config import discover, load_config, natural_key, normalize_config, resolve_selections
from mdwb.topology import _amber_flags, inspect_topology, residue_range


def amber_metadata(residues, periodic=0):
    """Produce fixed-width metadata with genuinely contiguous four-char fields."""
    labels, atoms, residue_starts = [], [], []
    for label, atom_names in residues:
        labels.append(label)
        residue_starts.append(len(atoms) + 1)
        atoms.extend(atom_names)
    pointers = [0] * 31
    pointers[0], pointers[11], pointers[27] = len(atoms), len(labels), periodic
    chunks = ["%VERSION  VERSION_STAMP = V0001.000  DATE = 09/05/26\n"]
    for flag, values, kind, width, per_line in (
        ("POINTERS", pointers, "I", 8, 10),
        ("ATOM_NAME", atoms, "a", 4, 20),
        ("RESIDUE_LABEL", labels, "a", 4, 20),
        ("RESIDUE_POINTER", residue_starts, "I", 8, 10),
    ):
        chunks.append(f"%FLAG {flag}\n%FORMAT({per_line}{kind}{width})\n")
        for offset in range(0, len(values), per_line):
            group = values[offset:offset + per_line]
            chunks.append("".join(f"{value:>{width}}" if kind == "I" else f"{value:<{width}}" for value in group) + "\n")
    # Unrelated floating-point flags must not leak into stored atom names.
    chunks.append("%FLAG CHARGE\n%FORMAT(5E16.8)\n  0.00000000E+00\n")
    return "".join(chunks)


class WorkspaceCase(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="mdwb_metadata_")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def topology(self, residues=None, periodic=0, filename="system.prmtop"):
        path = self.root / filename
        path.write_text(amber_metadata(residues or [("ALA", ["N", "CA", "C", "O", "H"])] , periodic), encoding="ascii")
        return path

    def config(self, **updates):
        self.topology()
        (self.root / "rep1.nc").write_bytes(b"coordinate placeholder")
        config = {"topology": "system.prmtop", "replicas": [{"name": "rep1", "trajectories": ["rep1.nc"]}]}
        config.update(updates)
        return config


class TopologyTests(WorkspaceCase):
    def test_fixed_width_amber_fields_are_not_whitespace_split(self):
        residues = [("NALA", ["N", "CA", "C", "O", "1HB"]), ("CHID", ["N", "CA", "C", "O", "HD1"])]
        path = self.topology(residues)
        flags = _amber_flags(path)
        self.assertEqual(flags["RESIDUE_LABEL"], ["NALA", "CHID"])
        self.assertEqual(flags["RESIDUE_POINTER"], ["1", "6"])
        self.assertEqual(flags["ATOM_NAME"], [atom for _, names in residues for atom in names])
        self.assertNotIn("CHARGE", flags)
        info = inspect_topology(path)
        self.assertEqual(info["natom"], 10)
        self.assertEqual(info["residues"][1]["atoms"], residues[1][1])
        self.assertEqual([r["index"] for r in info["residues"]], [1, 2])

    def test_protein_dna_rna_termini_water_and_ions_keep_topology_indices(self):
        residues = [
            ("NALA", ["N", "CA", "C", "O"]),
            ("CHID", ["N", "CA", "C", "O"]),
            ("DA5", ["C4'", "C1'", "N9"]),
            ("DT3", ["P", "C4'", "C1'", "N1"]),
            ("RA5", ["C4'", "C1'", "O2'", "N9"]),
            ("RU3", ["C4'", "C1'", "O2'", "N1"]),
            ("WAT", ["O", "H1", "H2"]),
            ("HOH", ["O", "H1", "H2"]),
            ("Na+", ["Na+"]),
            ("CL", ["Cl-"]),
            ("LIG", ["C1", "O1"]),
        ]
        info = inspect_topology(self.topology(residues, periodic=2))
        self.assertTrue(info["periodic"])
        self.assertEqual([r["kind"] for r in info["residues"]], ["protein", "protein", "dna", "dna", "rna", "rna", "water", "water", "ions", "ions", "other"])
        self.assertEqual(info["selections"]["protein"], ":1-2")
        self.assertEqual(info["selections"]["dna"], ":3-4")
        self.assertEqual(info["selections"]["rna"], ":5-6")
        self.assertEqual(info["selections"]["complex"], ":1-6")
        self.assertEqual(info["selections"]["solute"], ":1-6,11")
        self.assertEqual(info["selections"]["water"], ":7-8")
        self.assertIn("@N,CA,C,O", info["selections"]["backbone"])
        self.assertIn("@C4'", info["selections"]["representative"])
        self.assertIn("11", " ".join(info["warnings"]))

    def test_unknown_backbone_residues_are_inferred_with_a_warning(self):
        residues = [("MOD", ["N", "CA", "C", "O"]), ("8OG", ["C1'", "C2'", "C3'", "C4'", "O4'", "N9"])]
        info = inspect_topology(self.topology(residues))
        self.assertEqual([r["kind"] for r in info["residues"]], ["protein", "dna"])
        self.assertTrue(any("inferred" in warning.lower() and "1-2" in warning for warning in info["warnings"]))

    def test_nonperiodic_amber_is_distinguished_from_unknown_format(self):
        self.assertIs(inspect_topology(self.topology())["periodic"], False)
        path = self.root / "system.psf"
        path.write_text("PSF\n", encoding="ascii")
        info = inspect_topology(path)
        self.assertIsNone(info["periodic"])
        self.assertEqual(info["selections"], {})
        self.assertTrue(info["warnings"])

    def test_incomplete_amber_metadata_is_rejected(self):
        path = self.root / "incomplete.prmtop"
        path.write_text("%VERSION VERSION_STAMP = V0001.000\n%FLAG POINTERS\n%FORMAT(10I8)\n       1\n", encoding="ascii")
        with self.assertRaisesRegex(ValueError, "Incomplete Amber"):
            inspect_topology(path)

    def test_gzip_amber_detection_preserves_atom_identity(self):
        plain = self.topology()
        compressed = self.root / "system.prmtop.gz"
        with gzip.open(compressed, "wt", encoding="ascii") as handle:
            handle.write(plain.read_text(encoding="ascii"))
        self.assertEqual(inspect_topology(compressed), inspect_topology(plain))

    def test_legacy_sugar_atom_names_generate_exact_numbered_masks(self):
        info = inspect_topology(self.topology([("DA", ["P", "C4*", "C1*", "N9", "C4", "H4*"])]))
        self.assertEqual(info["selections"]["nucleic_backbone"], "@1-2")
        self.assertEqual(info["selections"]["representative"], "(@2)")
        self.assertIn("Legacy", " ".join(info["warnings"]))

    def test_inconsistent_natom_is_rejected_before_masks_are_suggested(self):
        path = self.topology()
        payload = path.read_text(encoding="ascii")
        payload = payload.replace("       5", "       9", 1)
        path.write_text(payload, encoding="ascii")
        with self.assertRaisesRegex(ValueError, "Inconsistent Amber"):
            inspect_topology(path)

    def test_residue_range_compacts_sorted_unique_indices(self):
        self.assertEqual(residue_range([8, 3, 2, 8, 4, 1, 10]), "1-4,8,10")
        self.assertEqual(residue_range([]), "")


class DiscoveryTests(WorkspaceCase):
    def test_discovery_naturally_sorts_segments_without_claiming_replicas(self):
        for name in ("segment10.nc", "segment2.nc", "segment1.nc", "SYSTEM.PRMTOP"):
            (self.root / name).write_bytes(b"placeholder")
        ignored = self.root / ".git"
        ignored.mkdir()
        (ignored / "hidden.nc").write_bytes(b"placeholder")
        found = discover(self.root)
        self.assertEqual([p.name for p in found["trajectories"]], ["segment1.nc", "segment2.nc", "segment10.nc"])
        self.assertEqual([p.name for p in found["topologies"]], ["SYSTEM.PRMTOP"])
        self.assertEqual(set(found), {"topologies", "trajectories"})
        self.assertEqual(sorted(["rep10", "rep2", "rep1"], key=natural_key), ["rep1", "rep2", "rep10"])


class ConfigurationTests(WorkspaceCase):
    def test_default_time_is_unknown_and_input_data_is_not_mutated(self):
        original = self.config()
        before = copy.deepcopy(original)
        normalized = normalize_config(original, self.root)
        self.assertIsNone(normalized["frames"]["dt_ps"])
        self.assertEqual(original, before)
        self.assertTrue(Path(normalized["topology"]).is_absolute())

    def test_declared_interval_and_stride_remain_separate(self):
        config = self.config(frames={"dt_ps": 2.5, "start": 11, "stop": 61, "stride": 5})
        frames = normalize_config(config, self.root)["frames"]
        self.assertEqual(frames["dt_ps"], 2.5)
        self.assertEqual(frames["stride"], 5)
        self.assertEqual(frames["start"], 11)

    def test_nonphysical_frame_intervals_are_rejected(self):
        for interval in (0, -1, float("nan"), float("inf"), "2 ps", True):
            with self.subTest(interval=interval):
                with self.assertRaises(ValueError):
                    normalize_config(self.config(frames={"dt_ps": interval}), self.root)

    def test_invalid_frame_ranges_and_boolean_stride_are_rejected(self):
        for frames in ({"start": 0}, {"start": 8, "stop": 7}, {"stride": True}, {"stride": 1.5}):
            with self.subTest(frames=frames):
                with self.assertRaises(ValueError):
                    normalize_config(self.config(frames=frames), self.root)

    def test_explicit_segment_order_is_preserved_not_resorted(self):
        (self.root / "segment10.nc").write_bytes(b"placeholder")
        config = self.config(replicas=[{"name": "rep1", "trajectories": ["segment10.nc", "rep1.nc"]}])
        normalized = normalize_config(config, self.root)
        self.assertEqual([Path(p).name for p in normalized["replicas"][0]["trajectories"]], ["segment10.nc", "rep1.nc"])

    def test_duplicate_paths_are_rejected_within_and_across_replicas(self):
        combinations = [
            [{"name": "rep1", "trajectories": ["rep1.nc", "./rep1.nc"]}],
            [{"name": "rep1", "trajectories": ["rep1.nc"]}, {"name": "rep2", "trajectories": [str(self.root / "rep1.nc")]}],
        ]
        for replicas in combinations:
            with self.subTest(replicas=replicas):
                with self.assertRaisesRegex(ValueError, "more than once"):
                    normalize_config(self.config(replicas=replicas), self.root)

    def test_duplicate_replica_names_are_rejected(self):
        (self.root / "rep2.nc").write_bytes(b"placeholder")
        replicas = [{"name": "rep1", "trajectories": ["rep1.nc"]}, {"name": "rep1", "trajectories": ["rep2.nc"]}]
        with self.assertRaisesRegex(ValueError, "Duplicate replica"):
            normalize_config(self.config(replicas=replicas), self.root)

    @unittest.skipUnless(os.name == "nt", "Windows paths are case-insensitive")
    def test_windows_case_alias_does_not_create_a_second_replica(self):
        replicas = [{"name": "rep1", "trajectories": ["rep1.nc"]}, {"name": "rep2", "trajectories": ["REP1.NC"]}]
        with self.assertRaisesRegex(ValueError, "more than once"):
            normalize_config(self.config(replicas=replicas), self.root)

    def test_output_cannot_contain_topology_or_trajectory(self):
        for output in (".", "..", "system.prmtop", "rep1.nc"):
            with self.subTest(output=output):
                with self.assertRaisesRegex(ValueError, "must not contain input"):
                    normalize_config(self.config(output=output), self.root)

    def test_output_cannot_contain_reference(self):
        results = self.root / "analysis_results"
        results.mkdir()
        (results / "common.rst7").write_bytes(b"reference placeholder")
        with self.assertRaisesRegex(ValueError, "must not contain input"):
            normalize_config(self.config(reference="analysis_results/common.rst7"), self.root)

    def test_separate_output_subdirectory_is_allowed(self):
        normalized = normalize_config(self.config(output="results/rmsd"), self.root)
        self.assertEqual(Path(normalized["output"]), self.root / "results" / "rmsd")
        self.assertFalse(Path(normalized["output"]).exists())

    def test_config_paths_resolve_relative_to_json_location(self):
        payload = self.config()
        config_path = self.root / "settings.json"
        config_path.write_text(json.dumps(payload), encoding="utf-8")
        normalized = load_config(config_path)
        self.assertEqual(Path(normalized["topology"]), self.root / "system.prmtop")
        self.assertEqual(Path(normalized["replicas"][0]["trajectories"][0]), self.root / "rep1.nc")

    def test_missing_trajectory_does_not_silently_skip_samples(self):
        with self.assertRaisesRegex(ValueError, "Input file not found"):
            normalize_config(self.config(replicas=[{"name": "rep1", "trajectories": ["missing.nc"]}]), self.root)

    def test_custom_residue_section_does_not_reindex_topology(self):
        config = normalize_config(self.config(sections=[{"name": "domain", "mask": ":12-25@CA"}], fit_mask="ca"), self.root)
        selections, info = resolve_selections(config)
        self.assertEqual(config["sections"][0]["mask"], ":12-25@CA")
        self.assertEqual(config["fit_mask"], selections["ca"])
        self.assertEqual(config["imaging"]["mode"], "off")
        self.assertEqual(info["natom"], 5)

    def test_unknown_matrix_alias_fails_before_cpptraj(self):
        config = normalize_config(self.config(advanced={"matrix_mask": "representativ_typo"}), self.root)
        with self.assertRaisesRegex(ValueError, "Unknown selection"):
            resolve_selections(config)

    def test_unknown_nested_field_does_not_silently_use_a_default(self):
        with self.assertRaisesRegex(ValueError, "Unknown frames"):
            normalize_config(self.config(frames={"dt": 2}), self.root)

    def test_unrecognized_topology_with_explicit_masks_remains_usable(self):
        raw = self.config()
        raw.update(topology="system.psf", fit_mask=":1-4@CA", sections=[{"name": "domain", "mask": ":1-4"}], imaging={"mode": "off"})
        (self.root / "system.psf").write_text("PSF\n", encoding="ascii")
        config = normalize_config(raw, self.root)
        resolve_selections(config)
        self.assertEqual(config["advanced"]["matrix_mask"], ":1-4@CA")

    def test_explicit_analysis_and_nucleic_options_are_accepted(self):
        config = self.config(analysis={"enabled": ["rmsd", "nastruct"]}, nucleic={"resrange": "1-12,15-26", "resmap": {"8OG": "G"}}, interactions=[{"name": "interface", "mask1": ":1-5", "mask2": ":6-10", "contact_cutoff": 4.0}])
        normalized = normalize_config(config, self.root)
        self.assertEqual(normalized["nucleic"]["resmap"], {"8OG": "G"})
        self.assertEqual(normalized["analysis"]["enabled"], ["rmsd", "nastruct"])


if __name__ == "__main__":
    unittest.main()
