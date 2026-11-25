# Lane NB: the `numbering` kind

**Goal** (the author, 2026-11-19: "Add it now (0.4.0)"; lane DA stopped because the language cannot declare a
numbering): a `numbering` kind as [locations-and-queries](../../../design/locations-and-queries.md), "Numberings in
forms" says (a `numbering` child of a sequence: offset, segments, labels, a scheme), built with the `labels` method
only this slice:
- **The node**: `numbering` under a `sequence` (its id the numbering's name), params `method: labels`, `labels` (a
  list of strings, one per position, inline) or `source` + a field (a records source's column, FS-A1's runtime, a row a
  position or keyed by position), `scheme` (an open word: `kabat`, `imgt`, `chothia`, `uniprot`, …). Labels need no
  order (IMGT's out of index order, EC-7), carry insertion codes (`100A`, `111.1`), and may be absent at a position
  (`null`: no label there).
- **Resolve**: the model's `numberings` gains the declared ones beside `index` and `ordinal`; `numbering(system,
  name)` returns them (`resolve/model.ts` ~460, ~508, ~924; `resolve/query-model.ts` ~37), so the existing query,
  address and ruler code reads them: `(at :label "100A" :numbering kabat)`, `(range :from "95" :to "102"
  :numbering kabat)`, the ruler layer's `numbering`, hover showing the label, an address printed in a numbering. The
  label rules EC-1 to EC-8 as the design and `query/numbering.ts` already state.
- **Checks**: a label list's length equals the sequence's (else an issue), duplicate labels reported, an unknown scheme
  word accepted (open), the language description and builder updated; the language stays `0.4.0` (unreleased: an
  addition, recorded in its Versions).
- The core `_spec`'s Language and Queries, `features.md`'s row ("the language declares none yet"), `files.md`'s
  mentions ("no `numbering` kind yet") and `molecules.md` updated where they say the kind is missing.

**Tests** (fail on the parent): a sequence with Kabat and IMGT numberings: `(at :label …)` in each, a range across
insertion codes, IMGT out of index order, a missing label snapping as EC rules say, the ruler in each numbering in
canvas and SVG, hover's label, an address printed back in a numbering; labels from a records source.

**Not in scope**: other methods (offset, segments), the Antibody dataset (DA, after you), anything else. DG, DC, DT
run beside you.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/nb-numbering
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core size`. Exit codes
explicit (`; echo exit=$?`), non-zero test counts.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.scaffolds/seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix, or the kind needs more than the params above.

**Report**: commits; checks; sizes; drift.
