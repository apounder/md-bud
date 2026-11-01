# Examples

Open [rendered/index.html](rendered/index.html) in a local browser; download the
file first when browsing this repository on GitHub. No server is needed.

The report contains three unequal **synthetic** replicas, circular torsions,
state fractions and matrix data. Click **Load 3D structure** on the overview to
inspect experimental crambin (PDB 1CRN). Its structure is unrelated to the
synthetic curves. Attribution is in [THIRD_PARTY.md](../THIRD_PARTY.md).

Rebuild the example using an existing Python environment with NumPy:

```text
python examples/build_demo.py
```

The JSON files are configuration templates: replace their topology, trajectory,
residue and output paths before use. Start with `protein.json`, or use the wizard
to construct selections appropriate for your own system.
