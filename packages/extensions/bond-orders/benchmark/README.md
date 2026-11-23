# Bond-order perception benchmark

Every CCD component perceived from coordinates and compared with its CCD definition by InChI (RDKit, `rdkit-esm`). The
plan: [../plans/benchmark.md](../plans/benchmark.md).

Not a workspace package and never published: a workspace of its own (`pnpm-workspace.yaml`), so a root `pnpm install`
does not install it and RDKit never reaches Mol\* developers or users (plan §2).

## Run

```sh
npx tsc -b                     # at the repository root: the benchmark runs the built lib of the extension and Mol*
cd packages/extensions/bond-orders/benchmark
pnpm install                   # this folder alone (its own lockfile); installs rdkit-esm
pnpm typecheck
node src/run.ts --sample 1000  # tiers 0, A and B on 1,000 components spread over the CCD
node src/run.ts --filter ATP,CHR --tier A
node src/run.ts --hydrogens kept --tier A,B   # the CCD's H kept in the input
node src/fetch.ts                              # tier C's entries ahead of a run (about 2 GB, gzipped)
node src/run.ts --tier C                       # the linked components in their representative entry
node src/report.ts --tier C                    # the HTML report again, from the last run's files
node src/run.ts --dataset labute2005           # a dataset's pairs, in their own entries
node src/baseline.ts --tier 0,A,B,C            # the last full runs become the baselines (commit them; rows kept in out/baseline/)
node src/baseline.ts --dataset labute2005
```

Options: `--tier` (`0,A,B` by default; `C` with a dataset), `--dataset` (`ccd`, the whole CCD, by default; or a file of
`datasets/`), `--hydrogens` (`stripped` by default, or `kept`), `--sample N`, `--filter ID,...`, `--cache` (`cache/`),
`--out` (`out/`), `--viewer` (the report's viewer links: the bond-providers example's dev server by default,
`pnpm dev:examples -- bond-providers`).

- `cache/components.cif`: the CCD, downloaded once from wwPDB and indexed by component (`components.index.json`).
- `cache/rcsb/bound.json`, `cache/entries/<ID>.bcif.gz`: tier C's inputs, the components each representative entry lists
  as bound (RCSB data API) and the candidates' entries (`models.rcsb.org`, whole files, gzipped).
- `out/<dataset>/<connectivity>/<hydrogens>/`, below `out/<…>/`: `ccd` the dataset of the whole CCD, `ccd` the
  connectivity (the CCD's, plan §3.3). A `--filter` or `--sample` run goes to `<hydrogens>-filter/` or
  `<hydrogens>-sample/` (`report --selection filter` draws it again), never over a full run.
- `out/<…>/<tier>/<ID>.sdf`, `<ID>.ref.sdf`: the perceived and reference mol blocks of each component (a dataset's row:
  `<entry>_<component>`).
- `out/<…>/<tier>.json`: counts per bucket and level, run time, peak memory, every component not passing (or passing
  with a diagnostic) with its reason, InChI, L0 diff and diagnostics, and the diff against the baseline.
- `out/<…>/<tier>/index.html`: the report, opened from disk. One row per component not passing (or passing with a
  diagnostic): the reference and perceived drawings on the same 2D coordinates (`svg/`), the bonds that differ at L0
  highlighted, a link to the viewer; filters by bucket, link kind, covalent ligand, diagnostics, regression and text.
  Excluded components in a table of their own.
- `baseline/<dataset>/<connectivity>/<hydrogens>/<tier>.json`, committed: the counts per bucket, and for a dataset the
  rows not passing with their bucket (the whole CCD's stay local: `out/baseline/`, written by `baseline.ts` too). A full
  run is compared with its baseline's counts, and every run by row when the rows are at hand: regressions (pass < known
  < fail, unsanitizable < error), improvements, and rows in or out of `excluded`, printed and in `<tier>.json`.
- `datasets/<name>.json`: a dataset from the literature, its source and its (entry, component) pairs, run in tier C's
  pipeline on the first complete copy, linked or not; a pair `"tier": "B"`, with its reason, judged on the CCD's model
  coordinates when they come from that entry (plan §6). `labute2005.json`: the 189 pairs of the v5 manifest.

## Tiers

| Tier | Input                                   | Pass: RDKit's InChI equals the CCD's on |
| ---- | --------------------------------------- | --------------------------------------- |
| 0    | CCD ideal coordinates and bonds, with H | formula, `/c`, `/h`, `/q`, `/p`         |
| A    | CCD ideal coordinates, no charges       | formula, `/c`, `/h`                     |
| B    | CCD model coordinates, no charges       | formula, `/c`, `/h`                     |
| C    | the linked copy in its entry, with caps | formula, `/c`, `/h`                     |

In tiers A and B, H are removed (`--hydrogens stripped`) or kept (`kept`: the doubles placed with H present, NAD+
against NADH; protonation and charge are not judged, the layers compared are the same). Stereo layers are not compared
(plan §5). `manifest.json` lists known cases with their reason (`known`, not failures), by component or by a dataset's
row (`<entry>_<component>`: the v5 `skipBonds` reasons of Labute 2005).

Tier 0 checks the pipeline (writer, charge model, RDKit) against the CCD's InChI: a mismatch there is a benchmark bug or
a CCD inconsistency. Tiers A and B compare with RDKit's identifiers of the same atoms with the CCD's orders
(`<ID>.ref.sdf`), at four levels: L0 every order equal, L1 the same canonical SMILES (another Kekulé form, a symmetric
resonance), L2 the same InChI (a mobile-H tautomer), L3 the same InChI with the `KET` and `15T` options (a keto-enol or
1,5 tautomer). Tier B leaves out the leaving atoms the model coordinates lack (OXT, OP3), counted apart
(`leavingMissing`). Each component also gets the diagnostics of plan §5.2 (no RDKit). Tier C measures the linked
components only (polymer residues, glycans, covalent ligands: an unlinked one is tier B's, on the same coordinates), on
their first linked and complete copy, with one cap atom per link; counts per link kind (`byKind`).

Requires Node 22.18 or later (TypeScript run without a build step).
