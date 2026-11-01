# Pi-stacking from residue numbers

Choose **Pi-stacking: discover and name rings from residue numbers** in the
wizard's additional-analysis menu, after selecting a basic/standard/advanced
workflow:

```text
python mdworkbench.py wizard --root /path/to/simulations
python mdworkbench.py run analysis_config.json
```

The portable `mdworkbench.pyz` offers the same workflow. Use the configuration
filename actually saved by your wizard.

## Confirming the systems

1. Enter one or more target **topology residue numbers**, for example `24,30-32`.
   These are the one-based indices in the workbench residue map, not necessarily
   PDB author residue numbers or MOL2 substructure IDs.
2. Review the discovered rings: residue, atom labels, topology atom indices,
   ring/fused-system type, and detection evidence are displayed. Choose the ones
   of interest and give each target a name. You do not enter atom names or masks.
3. Choose whether to include contacts between distinct rings in the same residue
   (off by default). Select the partner residue scope: all residues or a range.
4. Confirm the partner systems. Their automatic names can be customized too.
   Fused systems are preferred by default; individual constituent rings remain
   selectable. Overlapping atom sets are never paired against one another.
5. Set centroid distance, plane-angle tolerance, and geometric/mass centers.

Confirmed atom indices, names, residue numbers and criteria are saved in JSON.
Subsequent runs use exactly those selections; discovery is not repeated silently.
The topology must retain identical atom ordering across replicas. The run checks
that the masks select the expected atoms. Editing the saved topology or selecting
a different system requires renewed review.

## Discovery method and limits

Discovery follows the **bond graph**, without residue or atom-name templates:

- Amber prmtop/parm7: explicit bond records plus atomic numbers/masses and
  force-field atom types.
- MOL2: explicit bond records and element/atom-type information.
- PDB: CONECT bonds where provided. For residues without explicit bonds,
  within-residue heavy-atom connectivity is inferred from covalent distances
  and is explicitly labeled unverified. Partially specified residue connectivity
  is not silently repaired. Resolve alternate locations first.

The search finds small chordless 5–7-member rings containing C/N/O/S and connected
fused systems. Definite saturated carbon/nitrogen types screen out rings such as
ordinary sugars. Guanidinium-like C(N)3 connectivity is offered separately as an
**acyclic candidate**, reflecting the ARG groups in the supplied donor script.
Hydrogen atom names, residue aliases and modified-residue names do not define
membership. Large macrocycles, unusual elements and arbitrary acyclic conjugated
systems are outside this automatic search. Ring search limits bound work on large
or unusually dense residue graphs.

Connectivity and force-field type hints identify **candidates**, not a formal
aromaticity assignment. Bond order/electronic information in MD topologies can be
incomplete. Confirm uncertain candidates, especially coordinate-inferred PDB
rings; a planar-looking cycle alone is not evidence of aromatic chemistry. The
workflow does not install or execute a cheminformatics toolkit. See the
[RDKit aromaticity discussion](https://www.rdkit.org/docs/RDKit_Book.html#aromaticity)
for why ring membership and aromaticity are separate questions.

## Geometry and relationship to the donor script

The default contact definition follows `stacking-donor-analysis.sh`:

- centroid distance **<= 5 Å**;
- folded plane-normal angle **<= 30°**, where folded angle is
  `min(angle, 180° - angle)`.

This includes parallel and antiparallel normals and is independent of normal
sign. Both comparisons include the boundary. **Geometric centers** are the new
default; choose **mass centers** to reproduce the donor script's distance
convention. Thresholds are configurable. This defines a geometric contact;
there is no interaction-energy, ring-overlap/lateral-offset, or T-shaped-contact
criterion, and occupancy is not a lifetime or proof of an attractive interaction.

CPPTRAJ fits each confirmed system's plane and calculates pair distances with its
normal periodic-imaging behavior. Stacking batches preserve the periodic cell
orientation by using `nomod` for the global RMS check. Plane normals still require
whole, consistently imaged rings; a ring split across the periodic boundary must
be repaired in trajectory preparation. With no valid box, distances are Cartesian.

All confirmed target–partner pairs are evaluated over **every selected frame**.
There is no first-frame proximity screen or aromatic-residue allowlist, so a contact
that develops later is included and a never-stacked pair reports zero occupancy.
Pairs are unique; target–target contacts are not double-counted. Thirty-two pairs
share each CPPTRAJ trajectory pass, with ring vectors reused inside the batch.
`stacking.max_pairs` defaults to 2000 (maximum configurable 10000), and retained
vector/scalar arrays are checked against the existing memory budget.

The saved start/stop/stride and replica/segment grouping apply unchanged. Raw
geometry is retained, so saved `stacking.distance_cutoff` and `angle_cutoff` can be
adjusted before regenerating reports without rerunning coordinate analysis.
Changing centers or ring membership requires rerunning CPPTRAJ.

Command references: [plane vectors](https://amberhub.chpc.utah.edu/vector/),
[normal angles](https://amberhub.chpc.utah.edu/vectormath/),
[distances](https://amberhub.chpc.utah.edu/distance/).

## Occupancies and replica statistics

Open **Replica agreement and variability** in the HTML report. Dedicated outputs
in `reports/stacking/` are:

| File | Contents |
| --- | --- |
| `per_replica.dat` | Pair, target and residue occupancy fractions/percentages, frame counts and validity |
| `replica_summary.dat` | Equal-replica average occupancy and sample SD of replica occupancies |
| `pairs.dat`, `rings.dat` | Pair names and the confirmed system/atom definitions |
| `criteria.dat` | Thresholds, center convention and geometry definition |
| `any_contact_timeseries.dat` | Per-target/per-residue binary any-partner traces |
| `omissions.dat` | Missing or invalid replica geometry |

For replica `r`, occupancy is stacked frames divided by analyzed frames. Across
complete replicas, the mean weights each replica equally, regardless of length.
The SD uses `ddof=1` and describes **between-replica occupancy variation**. For
example, replicas with 20% and 80% occupancy have a mean of 50% and SD of about
42.43 **percentage points**. One replica has no estimable between-replica SD,
so that field is blank, not zero. This SD is not a framewise fluctuation SD or a
confidence interval.

A target-ring occupancy counts a frame once if **any** selected partner is
stacked. A residue occupancy similarly takes the union across its selected target
rings and partners. Simultaneous contacts and overlapping fused/constituent
representations therefore cannot make these union occupancies exceed 100%.

Missing/invalid geometry is not treated as an unstacked frame. Per-replica rows
identify incomplete data and show any fraction calculated over valid frames;
**incomplete replicas are excluded from aggregate mean/SD**. Requested and complete
replica counts are explicit. A union requires every contributing partner dataset
on the same frame grid. Zero observed contacts with complete geometry is a valid
zero result and remains in the mean/SD.

Pair occupancy, target/residue union occupancy, centroid distance and folded-angle
traces enter the ordinary plotting/distribution/convergence workflow. Binary
occupancy is a numeric indicator, so its running mean is a contact probability.
The earlier convergence reference choices and resource caps apply. Stable zero
occupancy or a stable running fraction alone does not establish adequate sampling.

## Validation

Synthetic tests cover arbitrary atom names, fused rings, saturated-ring exclusion,
PDB/MOL2 discovery, guanidinium candidates, cutoff/normal-sign behavior, all-pair
planning, unequal replica lengths, union occupancies, missing data and HTML/
convergence integration. Native CPPTRAJ execution against real trajectories is
still unverified in the current environment; see [validation](VALIDATION.md).
