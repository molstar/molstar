# Plan: bond-order perception benchmark over the CCD

Status: proposed, not started. A prerequisite of any change to the perception algorithm
([aromatic Kekulé](aromatic-kekulization.md) first).

## 1. Goal

Measure the algorithm, not test the code: every CCD component, perceived from coordinates, compared with its CCD
definition up to resonance and tautomerism. Report per tier and per component; detect regressions as a diff against a
recorded baseline; inspect a failure in a viewer. Unit tests (quick lane) keep only small hand-picked fixtures of
regressions found here.

The whole CCD is the primary target. A run can be narrowed to a named dataset from the literature (§6): Labute 2005
first, required for a publication.

Constraints:

- Nothing reaches users of the extension or of Mol\*: no RDKit, no data, no benchmark code in any install (§2).
- Little in the Mol\* repository: the benchmark's code, its manifest of known cases, the datasets and a small baseline.
  Every download is cached outside git.
- TypeScript only. A Python comparator (full RDKit: standardizer, resonance enumeration, xyz2mol) was weighed and
  dropped: InChI, the criterion, is the same library in `rdkit-esm`.
- mmCIF / BinaryCIF only: the newest components have 5-character codes, which the PDB format cannot hold, and PDB format
  files leave the archive in June 2027. The example's PDB-format loading (§7) is for `auto`, which needs connectivity
  without orders; the benchmark does not.

## 2. Packaging

Requirement: a root `pnpm install`, and an install of `@molstar/bond-orders-extension` or of `molstar`, never brings
`rdkit-esm`. Only someone who installs the benchmark gets it.

`packages/extensions/bond-orders/benchmark/`: a subfolder of the extension, next to the code it measures. Not a Mol\*
package: not a root workspace member (`packages/extensions/*` matches direct children only), not in `inventory.json`,
not in the root `tsconfig.json`, never published (the extension's `files` is `lib` and `src`).

- **Installed on its own**: a workspace of its own, `benchmark/pnpm-workspace.yaml` (`packages: []`; pnpm takes the
  nearest one): `pnpm install` in `benchmark/` writes a committed `benchmark/pnpm-lock.yaml` and installs it alone
  (probed). Root commands (`pnpm install --frozen-lockfile`, `-r`, `--filter`) never see it. pnpm 12 would install it
  alone without the file too; the file is what marks the folder as another workspace for the repository checks.
- **`benchmark/package.json`**: private; `"@molstar/bond-orders-extension": "link:.."`, `rdkit-esm`, and a `link:` for
  each `@molstar/*` package it imports itself. `link:`, not `file:`: `file:` would resolve the extension's `workspace:*`
  and `catalog:` outside the root workspace and fail. Node resolves the extension's imports from its real path (the root
  install's `node_modules`) and loads `lib`: a full `npx tsc -b` first, like the heavy-lane CLIs.
- **Committed**: the code, a README (how to run), `manifest.json` (exclusions and known cases with reasons),
  `datasets/`, `baseline/` (§6), a `.gitignore` for `cache/` and `out/`.
- **Repository checks**: `scripts/workspace/check.mjs` (`packageManifestPaths`) and `scripts/workspace/import-graph.mjs`
  skip any folder with a `pnpm-workspace.yaml` of its own (`isOtherWorkspace`, a generic rule: no path of this package
  in the shared scripts); without it the benchmark would fail as missing from the inventory. No `*.test.ts` in
  `benchmark/`: the quick lane would import RDKit. `npx tsc -b` is unaffected (the extension's tsconfig includes
  `src/**`); Biome and Prettier format the folder like the rest of the repository.
- Rejected: a devDependency of the extension or a workspace member (every Mol\* developer gets RDKit), an optional
  dependency (installed by default), an optional peer (in the published manifest, though the extension does not use
  RDKit).

## 3. Inputs

### 3.1 The CCD

`components.cif.gz` (wwPDB), downloaded once into `cache/`. Per component: `chem_comp_atom` (ideal and model
coordinates), `chem_comp_bond` (the reference, and the connectivity of §3.3), `pdbx_chem_comp_descriptor` (the CCD's
InChI), `chem_comp.pdbx_model_coordinates_db_code` (the representative entry) and `pdbx_model_coordinates_missing_flag`.

### 3.2 The representative instance (tier C)

Tier B already has the representative entry's coordinates (the CCD's model coordinates come from it). What the entry
adds is context, and perception reads only one part of it: the covalent links (a bond leaving the residue counts in the
valence, with its partner's element). Non-bonded neighbours are never read. So tier C measures the **linked** components
only, in the entry named by `pdbx_model_coordinates_db_code`; an unlinked one is tier B's.

- **Candidates**: the polymer types (`chem_comp.type` with `LINKING`, the saccharides) and the non-polymer components
  listed in their entry's `rcsb_entry_info.nonpolymer_bound_components` (RCSB data API, GraphQL, `entries(entry_ids:)`
  up to 1,000 per request; cached). That list counts metal coordination as bound (an ion, a heme's Fe–His): it
  over-selects, and the file decides (next). Glycans are branched entities, mostly absent from it: taken by type.
- **Download**: the entry's static BinaryCIF (`models.rcsb.org/<id>.bcif`, about 0.2 s), cached whole and gzipped in
  `cache/entries/` (BinaryCIF compresses to about a quarter: about 2 GB for the candidates); a few at a time, backing
  off on 429 and 5xx. An entry not served (obsolete, withdrawn: about 3%) excludes the component ("entry not served").
  The ModelServer's queries (`query-many`) run one after another at about 1 s each and return every copy: not used.
- **Linked**: the copy has a `struct_conn` row of type `covale` or a polymer bond. A candidate without (coordinated
  only, or free in that entry) is excluded from tier C ("not linked"): tier B covers it.
- **Link against the definition**: the CCD describes the unbound molecule. A link must replace one of the linked atom's
  H, or a leaving neighbour missing from the copy (condensations: peptide, glycosidic, phosphodiester bonds). A link
  that replaces nothing excludes the component ("link replaces no H or leaving atom of the CCD"): mostly an addition at
  an atom without H (a boronate, a carbonyl C: the CCD's orders are the unbound form's), otherwise a spurious link or a
  leaving atom the CCD does not flag. An addition at an atom with H (an acrylamide's CH2) is not caught: H counts cannot
  tell it from a substitution.
- **The copy**: the first linked and complete copy (model 1; the first altloc), recorded in the report (entry,
  `label_asym_id`, `label_seq_id` or `auth_seq_id`).
- **Complete**: every heavy atom of the CCD definition is present except leaving atoms
  (`chem_comp_atom.pdbx_leaving_atom_flag`). Missing leaving atoms are expected (OXT, OP3, the anomeric O1 of a linked
  residue) and handled by the comparison (§4.2). An entry with no complete copy excludes the component ("representative
  entry incomplete"), with no fallback to another entry.
- **Served files carry `chem_comp_bond` with orders** (archive mmCIF, static BinaryCIF, ModelServer alike): it is the
  connectivity of §3.3 and the reference's orders (§4.2). A component revised after the entry's release may differ from
  the CCD (ALIS7 in 2BLF): the file's graph for the component is compared with the CCD's, and a difference excludes it
  ("entry dictionary differs from the CCD"). The files carry no charges (`pdbx_formal_charge` is `?`), and
  `struct_conn.pdbx_value_order` is mostly `?`.

### 3.3 What perception sees

**Connectivity**: the CCD is authoritative for the bond graph without orders, which atom pairs are bonded. Perception
measures orders given that graph, Labute 2005's setting. Bonding by distance is a separate problem, measured on its own
so that its failures never count as order failures.

- **`--connectivity ccd`** (the default): `chem_comp_bond` kept. `forceCompute` copies the built-in bonds' topology and
  resets every order and aromatic flag to unknown (`perceive.ts`), so the CCD's orders never reach perception. In tier
  C, the entry's own `chem_comp_bond` (§3.2), checked against the CCD's graph.
- **`--connectivity distance`**: `chem_comp_bond` removed (the CCD file and the entry): Mol\* bonds by distance. Before
  any order is compared, the graph is compared with the CCD's by atom-name pairs; a difference puts the component in the
  `connectivity` bucket with the missing and extra bonds, its orders not judged.

**Hydrogens** (`--hydrogens`, a run option like the connectivity, not a tier):

- **`stripped`** (the default): perception runs on H-suppressed structures. Every H is removed before the model is
  built: the CCD's, and in tier C the deposited ones, protonation artifacts included.
- **`kept`**: the H stay (the CCD's; in tier C the deposited ones). It measures the doubles placed when H are present:
  NAD+ against NADH, a phosphate's P=O off its protonated O. Protonation and charge are not judged: the compared layers
  are those of `stripped` (§5); the diagnostics (§5.2) show a double on an atom whose stated H leave it no room
  (`O: valence 3`), which a mobile-H InChI hides.

An H absent from the CCD, which Mol\* bonds by distance, is never a reason to reject a component. Mol\* bonds by
distance every H pair missing from the component dictionary, so an H of the CCD's model coordinates within bonding
distance of a second atom gets a second bond: the reference is then unsanitizable and the component excluded.
Protonation stays outside the benchmark (TODO).

**Charges**: `atom_site.pdbx_formal_charge` / `chem_comp_atom.charge` removed: they would leak the reference into the
SDF (§4.1).

`struct_conn` stays in the entry: inter-residue links are part of the realistic input. Its `pdbx_value_order`, where
stated, is what real files carry; it orders the link bonds only, which are not compared (§4.2).

The removals are done on the parsed block, never on the file: a new `CifBlock` without the removed categories, the atom
category wrapped so that its charge field reads as absent and its H rows are gone; then `trajectoryFromCCD` /
`trajectoryFromMmCIF`. The viewer (§7) keeps its own copy of this function (it is in the workspace, the benchmark is
not); both follow this section.

## 4. Pipeline

Two steps, split at files on disk.

**perceive** (no RDKit), per component and tier:

1. Build the model from the CCD (the chosen coordinates) or the entry, with the removals of §3.3.
2. Perceive (`forceCompute`).
3. Write a V2000 SDF of the residue (§4.1) to `out/<connectivity>/<hydrogens>/<tier>/<comp>.sdf`: an artifact to open
   when a component fails.

**compare** (rdkit-esm), per SDF: parse and sanitize with `{ throwOnError: true }`, no autofix; compute the InChI and
the canonical SMILES; compare with the reference (§5) and put the component in one bucket. Writes
`out/<connectivity>/<hydrogens>/<tier>.json` and the HTML report (§4.3).

### 4.1 SDF writer

A V2000 writer of its own (`benchmark/src/sdf.ts`; tandemviz's is tied to its UI). Charges from
`calculateHydrogensCharge` (`@molstar/model/props/computed/chemistry/valence-model`) with `assignCharge: 'always'` (with
`'auto'` a file charge would win), plus boron, which that model never charges: B with 4 bonds is B- (its octet, as in
the extension's valence table). Leaving the input charges out (§3.3) keeps the reference out of the perceived side. A
molecule RDKit cannot sanitize (a C with five bonds) is a result: the `unsanitizable` bucket, with the message of
RDKit's rejection, never autofixed.

### 4.2 Reference, and the atoms compared

Both sides are compared on the same atoms: the component's heavy atoms present in the input, plus one cap per link.
Decided from the input, not from the component's type: no list of polymer types.

- **Reference**: the same model and atoms with the CCD's orders (the built-in bonds, no perception), through the same
  writer and charge model; RDKit computes its InChI and SMILES. The CCD's own InChI (`pdbx_chem_comp_descriptor`) is
  tier 0's reference only: tier 0 shows the two agree (every component compared in the first 1,000 sample), and one path
  covers the residues missing atoms (tier B's leaving atoms) as well.
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

`out/<connectivity>/<hydrogens>/<tier>/index.html`, static, opened from disk (`report.ts`, run after `compare` or on its
own from the last run's files). One row per component not passing, and per passing one with a diagnostic (§5.2),
filtered by bucket, link kind and text; the excluded components in a table of their own. Covalent ligands (the
non-polymer components linked in their entry, from tier C's run: `covalentLigands` in `C.json`) are flagged in every
tier, with their counts and a filter: their CCD definition is the unbound molecule. Each row has:

- the component, its bucket and level, regression or improvement against the baseline (§6);
- the expected and perceived 2D depictions side by side (RDKit SVG files in `<tier>/svg/`, loaded lazily: inlined, a
  tier's thousands of rows make a page of hundreds of MB). Both mol blocks are the writer's on the same atoms in the
  same order (§4.1): the reference's 2D coordinates go to the perceived side by atom index, and both are drawn
  unsanitized, with the file's own orders (an unsanitizable one too). The bonds that differ at L0 are highlighted, by
  their atoms' positions in the mol block (recorded with the L0 diff);
- a link to the viewer (§7) with the component's input: `?ccd=XYZ&coords=model` (tiers A, B) or
  `?pdb=1ABC&asym=B&seq=301` (tier C; `asym` the `label_asym_id`, `seq` the `auth_seq_id`: a non-polymer has no
  `label_seq_id`), plus `removeH=1` as the run's `--hydrogens` and `connectivity=distance` for that variant. The
  viewer's base URL is `--viewer` of `run` and `report`, the example's dev server by default.

The 2D drawings are the run's; the viewer is live (the current code).

## 5. Tiers and comparison

| Tier | Coordinates                               | H                         | Pass criterion (InChI layers compared) |
| ---- | ----------------------------------------- | ------------------------- | -------------------------------------- |
| 0    | CCD ideal, CCD bonds kept (no perception) | with                      | formula, `/c`, `/h`, `/q`, `/p`        |
| A    | CCD ideal                                 | `--hydrogens`             | formula, `/c`, `/h`                    |
| B    | CCD model                                 | `--hydrogens`             | formula, `/c`, `/h`                    |
| C    | linked copy in its entry (§3.2)           | `--hydrogens` (deposited) | formula, `/c`, `/h`                    |

- **Tier 0** checks the pipeline (writer, charge model, RDKit version) before perception is measured: a mismatch there
  is a benchmark bug or a CCD inconsistency, never a perception failure.
- **Stereo is not compared** (`/b`, `/t`, `/m`, `/s`): the CCD's stereo does not always agree with its own coordinates
  (all 31 tier 0 mismatches of the first 1,000-component sample), and an order error shows in the formula or the `/h`
  layer anyway. Non-standard layers (`/f`, `/i`, `/r`, in some CCD InChIs) end the comparison.
- **Without H**, protonation and charge are unknowable: `/q` and `/p` are dropped on both sides. Standard InChI merges
  mobile-H tautomers (amidines, 2-pyridone / 2-hydroxypyridine, carboxylic acid / carboxylate) but not keto/enol or 1,5
  shifts (pyrazolone / hydroxypyrazole): those pass at L3.
- B and C share coordinates (the CCD's model coordinates come from the representative entry): C adds the links, the one
  part of the context perception reads, and measures the linked components only (§3.2).
- **Linked residues** (polymer residues, glycans, covalent ligands) are compared on the atom set of §4.2. Tier A has the
  free molecule (leaving atoms present, no link); tiers B and C often lack the leaving atoms, and C has the links. Tier
  C reports polymer residues (the easier case, most of the linked ones) and covalent ligands separately.
- Reported at four levels for every pass and failure: **L0** bond by bond by atom name (diagnostics: which bond), **L1**
  canonical SMILES (symmetric resonance: carboxylate, nitro, Kekulé forms), **L2** InChI as above (the criterion),
  **L3** InChI with the keto-enol and 1,5 tautomerism options (`KET`, `15T`), asked only when L2 differs: a tautomer of
  the reference is not an order failure. The geometry could tell some apart (a phenol perceived as a dienone loses its
  aromatic ring), so L3 is counted apart, for the Kekulé work to look at.

### 5.1 Buckets

`pass` (with its level: L0, L1, L2 or L3), `fail`, `unsanitizable`, `connectivity` (`--connectivity distance` only,
§3.3), `excluded` (with the reason). Exclusions, from the CCD alone, before any download:

- organometallics: any metal atom (HEM, chlorophylls, clusters), out of scope;
- single heavy atom, ions (counted on the atoms with the tier's coordinates);
- non-classical bonding in the CCD graph: B or C with more than 4 bonds, H with more than 1 (carboranes, bridging H);
- per tier, after `perceive`: a reference RDKit cannot sanitize (the CCD's own graph, e.g. a PF6 with a neutral P);
- atoms of unknown element;
- per tier: the tier's heavy-atom coordinates not complete (§3.2): ideal coordinates missing (A), model coordinates
  missing a heavy atom other than a leaving atom (B; checked per atom rather than by
  `pdbx_model_coordinates_missing_flag`, which may be set by a missing leaving atom: to verify), no complete copy in the
  representative entry (C).

`manifest.json` lists further cases by hand, each with a reason: a geometry the perception is allowed to fail on (no
wwPDB validation data is used; the manifest gives the explanation). These are reported as `known`, not counted as
failures.

### 5.2 Diagnostics (no RDKit)

Cheap checks on every perceived residue, reported alongside, independent of the reference:

- an atom's orders exceed its valence (H, B, C, N, O, F; the elements whose valence does not expand);
- an atom with an aromatic bond has more than one multiple bond;
- an aromatic C without a multiple bond: its ring left without its double (G1 of the aromatic plan), or an aromatic flag
  on a ring that is not aromatic (a saturated ring fused to an aromatic one), which neither InChI nor L0 sees.

Not yet: a set of `Aromatic` bonds whose written Kekulé form is not 4n+2 (G2).

## 6. Runs, datasets, baselines

- `run --tier A,B,C [--connectivity ccd|distance] [--hydrogens stripped|kept] [--filter ID...] [--dataset NAME]`:
  `perceive` then `compare`; resumable, downloads cached, concurrency bounded for the RCSB services. Its files go to
  `out/<dataset>/<connectivity>/<hydrogens>/` (`ccd` the dataset of the whole CCD), as the baselines; a `--filter` or
  `--sample` run to `<hydrogens>-filter/` or `<hydrogens>-sample/`, never over a full run's files.
- **Datasets**: `datasets/<name>.json`, the source cited and (entry, component) pairs. A dataset runs tier C's pipeline
  on its own entries, not the representative ones: the first complete copy, linked or not (a free ligand is what the
  literature measures), with one cap per link (§4.2). A row is `<entry>_<component>`: a component can come from several
  entries. Labute 2005 is the first: the 189 pairs (160 entries) of the v5 manifest (`origin/bonds-perception`,
  `manifest.labute2005.json`, read, never edited), for a comparison with the paper and the v5 counts; its `skipBonds`
  reasons move to `manifest.json`'s `known`, keyed by pair. Others from the literature follow the same format.
- **A pair judged at tier B** (`"tier": "B"` and a reason in the dataset file): the entry cannot be read as served, and
  the CCD's model coordinates come from that same entry (`pdbx_model_coordinates_db_code`), so tier B (§5) measures the
  same geometry without the context. For overlapping alternate copies, which the entry links to each other (spurious
  `struct_conn` rows, reported to RCSB): G37 in 1HTG, VAC in 4PHV. Nothing detects the case; the file states it, as
  `known` does. The run checks that the CCD's model entry is the pair's (an error otherwise: the CCD may change), and
  the report marks the row.
- `baseline/<dataset>/<connectivity>/<hydrogens>/<tier>.json`, committed: the counts per bucket. A dataset's also lists
  the IDs of its rows not passing, with their bucket (Labute 2005: under 1 KB, an example). The whole CCD's do not:
  thousands of IDs per tier (about 500 KB for the four) that few readers of the repository need. `baseline.ts` writes
  the committed file from a full run (`--filter` and `--sample` runs are refused), and a local copy with the rows in
  `out/baseline/` (ignored).
- **Diff**: each run against its baseline, when there is one, in `<tier>.json` (`diff`). A full run's counts against the
  baseline's. By row when the rows are at hand (committed, or the local copy if its counts are the committed ones): a
  row whose bucket got worse (pass < known < fail, unsanitizable < error) is a regression, better an improvement, in or
  out of `excluded` a change; listed with their L0 diff in the console, regressions marked in the HTML report. A
  `--filter` or `--sample` run is compared on its own rows, by row only. A change to the algorithm comes with its diff,
  and the baseline is updated in the same commit: whoever changes it diffs by row on their own machine, and the commit
  shows the counts.

## 7. Viewer: `examples/bond-providers`

The example is adapted to inspect failures (the maintainers' request). In the workspace, built from the sources
(`pnpm dev:examples -- bond-providers`): an edit to the algorithm shows on reload. No RDKit.

One page, two Prism viewers side by side, cameras linked (`syncCameras`, `@molstar/prism-viewer`). mmCIF / BinaryCIF
only (5-character codes); loading a local file comes later.

- **Header**: one field for a PDB ID or a CCD code. 4 characters starting with a digit read as a PDB ID, anything else
  as a CCD code; a PDB / CCD switch overrides the guess (the rule cannot tell every code). A CCD code has an ideal /
  model coordinates choice (tiers A, B). A PDB entry is the static BinaryCIF (`models.rcsb.org`), the CCD file
  `files.rcsb.org/ligands/download/<code>.cif`. Downloaded once, the bytes given to both viewers: the two sides read the
  same file.
- **Left, authoritative**: the file as served, built-in bonds. No controls.
- **Right, the provider under test**:
  - bond provider, `bond-order-perception` by default; the example's proof-of-concept providers stay in the list;
  - mode `auto` / `forceCompute` (perception only), `forceCompute` by default. `auto` is disabled, with a tooltip, when
    the model states its orders: `ComponentBond.Provider.get(model)` or `IndexPairBonds.Provider.get(model)`
    (`perceive.ts`). The bond table decides, not the format name; re-evaluated after each load. With BinaryCIF
    (`chem_comp_bond`) it is always disabled for now: `auto` is useful with files that state no orders;
  - Remove H (off by default: the provider is seen at work on what it is given) and Remove charges (on: the benchmark
    always removes them), applied on the parsed block (§3.3, the example's copy in `src/removals.ts`); model 1 and a
    residue's first altloc always;
  - the provider and the mode switch the bonds in place (`structure.bonds`); the removals re-parse the right side only,
    the camera kept.
- **URL**: follows the state, so a view is a link: `pdb=1ABC` or `ccd=XYZ`, `coords=model`, `asym`, `seq` (a copy to
  focus: `label_asym_id`, `auth_seq_id`), `bonds`, `mode`, `removeH=1`, `removeCharges=0`. Defaults are left out.
- **The report's links** (§4.3) fill it: `?ccd=197&coords=model&removeH=1`, `?pdb=1O4I&asym=B&seq=300&removeH=1`
  (`removeH` as the run's `--hydrogens`; perception in `forceCompute` and charges removed are the defaults). Not yet:
  `connectivity=distance`.
- **Deferred**: highlighting the differing bonds in 3D (the HTML report shows them in 2D); a perception trace
  (hybridization, functional groups, Hückel typing, the step that placed each double), which needs an opt-in hook in the
  extension's `State`.

## 8. Phases

One review per phase.

1. **Skeleton and pipeline check**: a pnpm probe in a scratch copy (root install, install in `benchmark/` with and
   without its workspace file, `link:..` resolution); `benchmark/` per §2, with the workspace-check skips;
   `components.cif.gz` downloaded and parsed; the exclusions of §5.1; the SDF writer; **tier 0** and **tier A** on a
   1,000-component sample, run time and memory recorded. Tier 0 mismatches are fixed or listed.
2. **Tier B**: complete coordinates checked per atom; buckets, levels L0 / L1 / L2, the diagnostics of §5.2; then the
   whole CCD.
3. **Tier C**: the candidates (data API), the downloads, the first linked and complete copy, the dictionary check, the
   caps (§4.2); a sample first, then everything.
4. **HTML report** (§4.3), with links to the viewer.
5. **Baselines, diff, datasets** (§6), Labute 2005 first; then `_spec §1.10`, `§2` and `§4`. From here the
   [aromatic Kekulé](aromatic-kekulization.md) plan can start.
6. **Viewer** (§7). It can go alongside phase 3 if inspecting failures in 3D is needed sooner.

## 9. Open questions

- The perception trace (§7): whether, and in what form.
