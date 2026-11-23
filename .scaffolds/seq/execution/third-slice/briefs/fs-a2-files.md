# Lane FS-A2: GFF3, the matrix read, p53's entry and the Files stories

**Goal**: what [files.md](../../../design/files.md) (agreed 2026-11-19) leaves after FS-A1 (merged: the records runtime
in `packages/seq/core/src/data/records.ts`), plus the matrix read DC's PAE needs ([grids-and-rows](../../../design/grids-and-rows.md),
"A grid over a complex"):

1. **The GFF3 preset** of the records reader (files.md, The file's own syntax): nine columns, headerless, `#` and `##`
   lines past, `##FASTA` and after ignored, `.` missing, base 1, ends inclusive, attributes as fields (percent-decoded,
   `ID` the `id`). **The cyclic wrap** in the row reader: a row ending past a cyclic system's length becomes two rows of
   one item; on a linear one it is clipped (`.range`). A test on a feature across pUC19c's origin drawn as one item in
   canvas and SVG.
2. **A matrix read**: an array of arrays (JSON or MessagePack) found by `records`'s path served as a grid's values
   (the `matrix` shape the grid layer reads), its size checked against the axes' total (an issue on a mismatch), NaN
   for `null`; a typed flat array with its dimensions kept without a copy where the protocol's conventions allow. The
   source param that says the shape (`shape: matrix`, or the role a grid reads) is a param within 0.4.0; name it in
   the `_spec`'s Language and say why.
3. **`dataFilter`'s trailing `(not …)` rule** (files.md, One track a feature type): the viewer's filter is only a
   trailing `(not …)` segment; a leading item form stays the layer's base and its rejected rows are not counted as
   filtered (the "n filtered" chip). Test: a leading item form kept, not counted.
4. **p53's UniProt entry**: P04637's REST entry JSON (`rest.uniprot.org/uniprotkb/P04637.json`, CC BY 4.0) and the
   Proteins API's variants, recorded once by hand with `packages/seq/examples/scripts/record.ts` (now this lane's;
   create it if it doesn't exist: fetch, trim to what the views read, never reshape, write with its `sources.json` row:
   URL, date, licence, attribution), trimmed JSON fixtures in the examples package. **`uniprotEntry(b, uri)`** in
   `seqx-examples` writes the source node (its transform) and one ranges layer a feature type (a leading item form).
   Anchors on the read data: 393 residues; R175, R248, R273 inside the DNA-binding domain. Tests never reach the
   network.
5. **`Features/Files`** stories (examples-and-stories: features under `Seq* (next)/Features/<Feature>`): p53's entry
   from JSON; the same table as TSV, CSV, JSON and MessagePack in both layouts drawing the same; a GFF3 of the
   plasmid's features; each printing the file's head and the source node beside the view; each with a headless twin
   reading the same files through `node-assets.ts`.
6. The `_spec`s (seqx-core Data › Records; the viewer's filter) and files.md's interaction matrix rows say what is
   built.

**Not in scope**: held trees and pairs from a file (FS-B); the PAE's grid over a complex (DC). HL and IZ run beside
you (highlights, navigators): stay out of them.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/fs-a2-files
<that commit>`. **The language stays 0.4.0.**

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`. Exit codes explicit (`; echo exit=$?`), non-zero
test counts.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.scaffolds/seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json` (write needed entries in the report), anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix, or a public API past the params above is needed.

**Report**: commits; checks; the fixtures recorded (URL, date, bytes before and after trimming); the stories to look
at; drift.
