# Third-party assets

The workbench source is MIT licensed. Bundled third-party files retain their own
licenses; no ownership of these files is claimed.

## 3Dmol.js 2.5.3

- Upstream: https://github.com/3dmol/3Dmol.js/tree/2.5.3
- Download: https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.3/3Dmol-min.js
- File: `mdwb/vendor/3Dmol-2.5.3-min.js` (524,009 bytes)
- SHA-256: `bcf422137a34b4e206413a8a83a67c24a1ace84e282bf3b164f6bda6aeb42a20`
- License: [upstream license and incorporated notices](mdwb/vendor/3Dmol-LICENSE.txt),
  retrieved from the pinned tag on 2026-09-05.

The pinned asset is embedded only when a reference structure is available and
executes only after the viewer button is clicked. Generated reports include its
license. No CDN connection is needed to open a report.

## Crambin example structure

`examples/1CRN.pdb` is the experimental structure [PDB 1CRN](https://www.rcsb.org/structure/1CRN),
deposited by W. A. Hendrickson and M. M. Teeter.
Structure identifier: https://doi.org/10.2210/pdb1CRN/pdb.
Downloaded from https://files.rcsb.org/download/1CRN.pdb on 2026-09-05.
The coordinate archive is provided under the wwPDB's open-access data policy:
https://www.wwpdb.org/about/usage-policies.

This structure demonstrates the viewer only. The example analysis curves are
synthetic and do not describe a simulation of crambin.
