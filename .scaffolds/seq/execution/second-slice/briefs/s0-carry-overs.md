# Lane S0: the first slice's carry-overs

**Goal**: three small items the first slice left (`plan.md`, stage 0).
1. **The Canvas2D backend's header comment** (`packages/seq/viewer/src/canvas/backend.ts`, the comment at the top) still
   says the label pass is drawn on the base canvas, with two canvases. The backend now has three stacked canvases:
   base (background, data), overlay (the emphasis passes) and labels above both. Make the comment say what the code
   and the viewer `_spec` ("Canvas2D backend") say.
2. **The Storybook package's tests**: the Vite plugin's test (`packages/seq/storybook/_test/vite.test.ts`) and the Gene
   by region test (`_test/gene-by-region.test.ts`) sit outside `src`. Move them where the package's convention and the
   root `AGENTS.md` put tests (a `_test/` folder next to the code they test), keeping them in the quick lane.
   `vitest` is not in the package's dev dependencies: ask for it in the report (a hot file).
3. **The TP53 exon coordinates** in `packages/seq/core/src/resolve/_test/p53-geometry.ts` (and the copy in
   `encodings.test.ts`) were written from memory. Check them against RefSeq NM_000546 on GRCh38 (chr17, minus strand):
   the coding exons 2 to 11, 0-based, half-open, the stop codon included, lengths adding to 3 × 393 + 3.
   - You may read NCBI's public pages or E-utilities to check (no personal data in a request). Tests never reach the
     network: the values stay in the fixture.
   - If they are wrong, correct them in both places, and say in the fixture's comment where they were checked (the
     accession and version, e.g. NM_000546.6). Run the tests that use them.
   - If you cannot check them (no network), mark the fixture's comment "unverified" with what to check, and say so in
     the report.

**Read first**: the root `AGENTS.md` and `packages/seq/AGENTS.md`; the viewer `_spec` ("Canvas2D backend").

**Model**: fast (sonnet).

**Base**: `packages-seq` at `6e3b15296`. Branch `seq/s0-carry-overs`; commit on it only, at each small green step.

**Touches**:
- `packages/seq/viewer/src/canvas/backend.ts` (the comment only);
- `packages/seq/storybook/_test/**` and wherever the tests move to;
- `packages/seq/core/src/resolve/_test/{p53-geometry.ts,encodings.test.ts}`.

Not in scope: anything else, including behavior changes.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json` (ask in the report);
- anything in `.scaffolds/seq/` and in `packages/core/`.

**Headers**: a file you modify keeps the Mol\* header, its year range extended to 2026 and `@author David Sehnal
<david.sehnal@gmail.com>` present. No mention of `packages/sequences` or `.scaffolds/seq/` in code, comments or tests.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a permission is refused (report the
action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; for item 3 what you checked against (the
accession and version, the source) and whether the values changed; the hot-file changes you need; followups.
