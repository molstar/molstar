# The dataset lanes (DG, DA, DP, DT, DC): what each does

A dataset lane builds one dataset of [examples](../../../design/examples.md) to its section there (the table of
**What, Source, Licence, Size, Form**, the **Why** and **Anchors**, and its **Views** table) and the dataset's row in
"Datasets" and "Views". Where `packages/seq/examples/src/<dataset>/` already holds a stand-in (generated letters, a
synthetic table), the lane replaces it with the dataset as the section says, keeping the stories that already exist
working (their twins updated to the real anchors).

**Each lane**:
1. **Fixtures**: record the published entries the section names once, by hand, with
   `packages/seq/examples/scripts/record.ts` (add the dataset's recorder there; fetch, trim to what the views read,
   never reshape; `sources.json` rows with URL, date, licence, attribution, sizes and checksums). JSON fixtures are
   stored compact (the formatter skips them); a file in a format the core reads stays in it (FASTA, Stockholm,
   Newick, TSV, CSV, JSON, MessagePack, GFF3); otherwise the `seq` format. What the section says is derived or
   generated is code in the dataset's module, generated data on core's `PCG` with a fixed seed. Licences exactly as the
   section says; a licence not stated is not taken. Tests never reach the network.
2. **Builders** in `seqx-examples` (`<dataset>/*.ts`) that write each view's spec from the fixtures (the transforms as
   source params; tracks, layers, themes as the views table says).
3. **The views** of the section's Views table, each a story under `Seq* (next)/Examples/<Dataset>/<View>` with its
   headless twin (reading the same files through `node-assets.ts`), drawing against the data (the section's anchors
   checked in canvas draw lists and SVG, both themes), and **an editing view** as its table says (a patch, undone).
   The dataset's index story lists its views and the Features shown on its data (examples-and-stories).
4. **Performance**: a view at the scale the section names keeps the frame targets of
   [performance](../../../design/performance.md) (an `itWithin` bound in the heavy lane; a quick representative).

**Every lane's rules**: the language stays 0.4.0; no new core API without stopping to ask (a gap found is reported,
not built around); reuse what is built (FS's records runtime, GFF3, the matrix read; RC's row columns; LB's chips; HL's
highlights; IZ's zoom strip); minimal code.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; the heavy bound alone (`pnpm test:heavy <file>`).
Exit codes explicit (`; echo exit=$?`), non-zero test counts.

**Look**: run Storybook from your worktree on the port your lane brief names (background; stop it after) and look at
every view in light and dark; say what you saw, and what a reader of the figure would find wrong.

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE
file for more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, year range
to 2026, that author added if missing. No mention of `.scaffolds/seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix, or a public API change is needed, or a licence is unclear.

**Report**: commits; fixtures (URL, date, licence, bytes before and after); each view and what it shows; checks; what
you saw; gaps against the section; drift.
