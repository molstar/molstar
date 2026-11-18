# Plan: bond-order perception benchmark over the CCD

Status: proposed, not started. A prerequisite of any change to the perception algorithm
([aromatic Kekulé](aromatic-kekulization.md) first).

## 1. Goal

Measure the algorithm, not test the code: every CCD component, perceived from coordinates, compared with its CCD
definition up to resonance and tautomerism. Report per tier and per component; detect regressions as a diff against a
recorded baseline; inspect a failure in a viewer. Unit tests (quick lane) keep only small hand-picked fixtures of
regressions found here.

Constraints:

- Nothing reaches users of the extension: no RDKit, no data, no benchmark code in the published package.
- Little in the Mol\* repository: the benchmark's code, its manifest of known cases and a small baseline. Every download
  is cached outside git.
- TypeScript only. A Python comparator (full RDKit: standardizer, resonance enumeration, xyz2mol) was weighed and
  dropped: InChI, the criterion, is the same library in rdkit-esm.
- mmCIF (text) only: PDB format files leave the archive in June 2027, and recent components have no PDB-format entry.
  Text over BinaryCIF: the cached files can be read when a component fails; the excerpts are small either way.

## 2. Packaging

`packages/extensions/bond-orders/benchmark/`: a private package with its own `package.json` and lockfile.

- **Not a workspace member**: `pnpm-workspace.yaml`'s `packages/extensions/*` matches direct children only, so the root
  `pnpm install` never installs RDKit. Installed by hand (`pnpm install` in `benchmark/`).
- **Not published**: the extension's `files` is `lib` and `src`.
- **Depends on** `@molstar/bond-orders-extension` through `link:..` and runs against its built `lib` (`npx tsc -b`
  first), like the heavy-lane CLIs; `rdkit-esm` as its own dependency.
- **Committed**: the code, a README (how to run), `manifest.json` (exclusions and known cases with reasons), `baseline/`
  (§6), a `.gitignore` for `cache/` and `out/`.
- To check at implementation: the root `tsconfig.json`, Biome and `pnpm check:workspace` (import graph, inventory)
  either cover the folder correctly or exclude it.
- Docs when it lands: `_spec §1.10` ("No benchmark here yet") rewritten to point at `benchmark/`, `_spec §2` (layout)
  and `§4` (tests: the quick lane keeps hand-picked fixtures only); the example's README (§7).

## 3. Inputs

### 3.1 The CCD

`components.cif.gz` (wwPDB), downloaded once into `cache/`. Per component: `chem_comp_atom` (ideal and model
coordinates, both with H), `chem_comp_bond` (the reference), `pdbx_chem_comp_descriptor` (the CCD's InChI),
`chem_comp.pdbx_model_coordinates_db_code` (the representative entry) and `pdbx_model_coordinates_missing_flag`.

### 3.2 The representative instance

One instance per component: the entry named by `pdbx_model_coordinates_db_code`, fetched as mmCIF from the ModelServer
(`encoding=cif`) with its surroundings (`residueSurroundings`, a few Å: links and neighbours, not the whole entry),
cached. The first complete copy of the component in that entry is used; the copy is recorded in the report (entry,
`label_asym_id`, `label_seq_id` or `auth_seq_id`).

**Complete** (every tier): every heavy atom of the CCD definition is present except leaving atoms
(`chem_comp_atom.pdbx_leaving_atom_flag`). Missing leaving atoms are expected (OXT, OP3, the anomeric O1 of a linked
residue) and handled by the comparison (§4.2); any other missing heavy atom excludes the component from the tier. In the
with-H tiers, H must be complete as well.

### 3.3 What perception must not see

Removed from every input before the model is built:

- `chem_comp_bond` (the CCD file and the entry): connectivity then comes from distances, orders are all unknown;
- `atom_site.pdbx_formal_charge` / `chem_comp_atom.charge`: charges would leak the reference into the SDF (§4.1).

`struct_conn` stays in the entry: inter-residue links are part of the realistic input. Its `pdbx_value_order`, where
stated, is what real files carry; it orders the link bonds only, which are not compared (§4.2).

The removal is done on the parsed block, never on the file: a new `CifBlock` without `chem_comp_bond`, the atom category
wrapped so that its charge field reads as absent; then `trajectoryFromCCD` / `trajectoryFromMmCIF`. The viewer (§7)
keeps its own copy of this function (it is in the workspace, the benchmark is not); both follow this section.

## 4. Pipeline

Two steps, split at files on disk.

**perceive** (no RDKit), per component and tier:

1. Build the model from the CCD (the chosen coordinates) or the entry, with the removals of §3.3, H removed in the
   H-stripped tiers.
2. Perceive (`forceCompute`).
3. Write a V2000 SDF of the residue (§4.1) to `out/<tier>/<comp>.sdf`: an artifact to open when a component fails.

**compare** (rdkit-esm), per SDF: parse and sanitize, no autofix; compute the InChI and the canonical SMILES; compare
with the reference (§5) and put the component in one bucket. Writes `out/<tier>.json` and the HTML report (§4.3).

### 4.1 SDF writer

Ported from `tandemviz-web` (`utils/writer/sdf.ts`, `sdfCharges.ts`) to v6 imports, inside the benchmark. Charges from
`calculateHydrogensCharge` (`@molstar/model/props/computed/chemistry/valence-model`) with `assignCharge: 'always'` (with
`'auto'` a file charge would win). Leaving the input charges out (§3.3) keeps the reference out of the perceived side. A
molecule RDKit cannot sanitize (a C with five bonds) is a result: its own bucket with RDKit's message, never autofixed.

### 4.2 Reference, and the atoms compared

Both sides are compared on the same atoms: those of the component present in the input, plus one cap per link. Decided
from the input, not from the component's type: no list of polymer types.

- **Reference**: the CCD graph (bonds, charges), restricted to the atoms present in the input. The CCD's InChI
  (`pdbx_chem_comp_descriptor`) is used only for a complete, unlinked residue with no leaving atom missing; otherwise
  RDKit computes the reference InChI from the restricted graph, through the same writer.
- **Caps**: each covalent link to another residue (`struct_conn`, the polymer bonds) is one cap atom on both sides: the
  partner atom's element (N across a peptide C, C across a glycosidic O, P across a phosphodiester O3'), neutral, its H
  filled by RDKit. The cap's bond is single on the reference side unless `struct_conn.pdbx_value_order` says otherwise
  (a PLP-Lys Schiff base: C=N); on the perceived side, the order the model has. The tandemviz writer's attachments
  (`includeConnected`) are written as `*`, which InChI rejects: they become cap atoms.
- **Perceived side**: the residue's atoms with the orders perceived in context (perception reads across links through
  the inter-unit neighbours), plus the same caps.

A residue without its leaving atoms and without links (model coordinates of a free amino acid without OXT) is compared,
on both sides, as the residue without them and no cap.

Limit: a cap is not the real neighbour. A failure caused beyond the first atom across a link does not show at L1/L2; the
L0 diff by atom name still lists it.

### 4.3 HTML report

`out/<tier>/index.html`, static, opened from disk. One row per component not passing (failures only: a few thousand rows
at most), with:

- the component, its bucket and level, regression or improvement against the baseline (§6);
- the expected and perceived 2D depictions side by side (RDKit SVG, inlined; the perceived one aligned on the expected
  one), the bonds that differ at L0 highlighted;
- a link to the viewer (§7) with the component's input: `?comp=XYZ&tier=C&entry=1ABC&asym=B&seq=301`. The viewer's base
  URL is an option of `compare`, the example's dev server by default.

The 2D drawings are the run's; the viewer is live (the current code).

## 5. Tiers and comparison

| Tier | Coordinates                               | H                          | Pass criterion                     |
| ---- | ----------------------------------------- | -------------------------- | ---------------------------------- |
| 0    | CCD ideal, CCD bonds kept (no perception) | with                       | InChI equals the CCD's             |
| A-H  | CCD ideal                                 | with                       | standard InChI                     |
| A    | CCD ideal                                 | stripped                   | InChI without `/q` and `/p` layers |
| B-H  | CCD model                                 | with                       | standard InChI                     |
| B    | CCD model                                 | stripped                   | InChI without `/q` and `/p` layers |
| C    | representative instance, in context       | as deposited, usually none | InChI without `/q` and `/p` layers |

- **Tier 0** checks the pipeline (writer, charge model, RDKit version) before perception is measured: a mismatch there
  is a benchmark bug or a CCD inconsistency, never a perception failure.
- **With H**, charges and protonation are determined: the full standard InChI must match. **Without H**, protonation and
  charge are unknowable: `/q` and `/p` are dropped. Standard InChI merges mobile-H tautomers (amidines, 2-pyridone /
  2-hydroxypyridine, carboxylic acid / carboxylate) but not keto/enol, a heavy-atom difference the geometry decides.
- B and C share coordinates (the CCD's model coordinates come from the representative entry): C adds the context (links,
  neighbours, altlocs, occupancy).
- **Linked residues** (polymer residues, glycans, covalent ligands) are compared on the atom set of §4.2. Tier A has the
  free molecule (leaving atoms present, no link); tiers B and C often lack the leaving atoms, and C has the links. Each
  component carries a `linked` flag, and the counts are reported separately for linked and unlinked residues: polymer
  residues are most of the linked ones, and the easier case.
- Reported at three levels for every pass and failure: **L0** bond by bond by atom name (diagnostics: which bond),
  **L1** canonical SMILES (symmetric resonance: carboxylate, nitro, Kekulé forms), **L2** InChI as above (the
  criterion).

### 5.1 Buckets

`pass` (with its level: L0, L1 or L2), `fail`, `unsanitizable`, `excluded` (with the reason). Exclusions, from the CCD
alone, before any download:

- organometallics: any metal atom (HEM, chlorophylls, clusters), out of scope;
- single heavy atom, ions;
- atoms of unknown element;
- per tier: the tier's coordinates not complete (§3.2): ideal coordinates missing (A), model coordinates missing a heavy
  atom other than a leaving atom (B; checked per atom rather than by `pdbx_model_coordinates_missing_flag`, which may be
  set by a missing leaving atom: to verify), no complete copy in the representative entry (C); H incomplete (with-H
  tiers).

`manifest.json` lists further cases by hand, each with a reason (e.g. a geometry the perception is allowed to fail on).
These are reported as `known`, not counted as failures.

### 5.2 Diagnostics (no RDKit)

Cheap checks on every perceived residue, reported alongside, independent of the reference:

- an atom's orders plus its H exceed its valence;
- an atom with an aromatic bond has more than one double;
- a set of bonds flagged `Aromatic` whose written Kekulé form is not 4n+2, or whose 1e atoms do not each carry exactly
  one in-ring double (what G1 and G2 of the aromatic plan break).

## 6. Runs, baseline, regressions

- `run --tier A,B,C [--filter COMP...] [--subset labute2005]`: `perceive` then `compare`; resumable, downloads cached,
  concurrency bounded for the RCSB services.
- `baseline/<tier>.json`, committed: the counts per bucket and the IDs of components not passing, with their bucket. Not
  the passing IDs (about 45k).
- `diff`: the current run against the baseline. Components whose bucket got worse are regressions, those that got better
  are improvements; both listed with their L0 diff (and marked in the HTML report). A change to the algorithm comes with
  its diff, and the baseline is updated in the same commit.
- **Labute 2005 subset**: the 177 (entry, component) pairs of the v5 manifest (`origin/bonds-perception`,
  `manifest.labute2005.json`), run through tier C's pipeline on their own entries (not the representative ones) for a
  comparison with the paper and the v5 counts. Its `skipBonds` reasons move to `manifest.json`.

## 7. Viewer: `examples/bond-orders`

The existing example becomes the comparison viewer. In the workspace, built from the sources
(`pnpm dev:examples -- bond-orders`): an edit to the algorithm shows on reload. No RDKit.

- **Dual view by default**: two Prism viewers side by side, cameras linked (`syncCameras`, `@molstar/prism-viewer`).

  |       | Left: default settings      | Right: perception                                      |
  | ----- | --------------------------- | ------------------------------------------------------ |
  | Data  | the file as served          | the same file, with the removals of §3.3               |
  | Bonds | built-in (`chem_comp_bond`) | `bond-order-perception`, `forceCompute`                |
  | Reads | the expected bonds          | what the benchmark measures (connectivity by distance) |

  The right side must load the stripped data, not only switch the provider: on the full file `forceCompute` resets the
  orders but keeps the table's connectivity, which the benchmark does not have.

- **Inputs**: the URL of a report link (§4.3; the CCD component at the tier's coordinates, or the entry focused on the
  copy) or the header (any entry, as today). mmCIF only: the PDB-format loading goes. The Default / auto / forceCompute
  buttons go: the two sides replace them.
- **Deferred**: highlighting the differing bonds in 3D (the HTML report shows them in 2D); a perception trace
  (hybridization, functional groups, Hückel typing, the step that placed each double), which needs an opt-in hook in the
  extension's `State`.

## 8. Open questions

- Representative entry incomplete: exclude the component (current proposal), or fall back to another entry?
- Geometry quality in tier C: label failures with the wwPDB validation (RSCC, Mogul outliers) as geometry-limited, or
  leave this to `manifest.json`?
- Run time and download volume over the whole CCD: measure on a sample first (1,000 components).
- The perception trace (§7): whether, and in what form.
