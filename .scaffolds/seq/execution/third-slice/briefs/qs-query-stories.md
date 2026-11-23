# Lane QS: Features/Query

**Goal**: [examples-and-stories](../../../design/examples-and-stories.md), "The query language in Features" (the author,
2026-11-19: "do we have storybook stories planned to showcase SeqQL functionality? Perhaps Features section?", then
"ok as recommended"). Today no story shows queries; the forms and addresses are built
([locations-and-queries](../../../design/locations-and-queries.md); `@molstar/seqx-core/query`: `resolveAddressIn`,
`addressFor`, the forms in `forms.ts`, the three result types; the viewer's `AddressModel`, `address.select`,
`address.highlight`).

Build `Seq* (next)/Features/Query` on p53's data (the merged UniProt entry and tracks):
1. **A query box story**: an input where an address or expression typed is parsed, checked and resolved as typed;
   its result selected (or highlighted, a toggle) in the view; the printed (canonical) form and the result's kind and
   size shown; issues with their hints when it doesn't resolve. Use the library's own resolve, print and issue paths,
   no private parser.
2. **Worked examples**, each a story with a short explanation and the query printed: items by category (the
   DNA-binding domain), a motif search in the letters, ranges and slices, a result carried from the protein onto its
   gene (p53's protein-on-gene view), `and`/`or`/`not` combined, a query over an alignment's rows and columns (the
   globins). Pick the forms that exist; where a design form is not built yet, leave it out and list it in the report.
3. Each story with its headless twin checking the result's location (positions, systems) and the printed form.
4. The examples Design page lists `Features/Query` beside the other Features; the core `_spec`'s Queries section names
   the stories as its examples.

**Not in scope**: new query forms or language changes (report what the stories needed and did not find). RC and LB
run beside you.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/qs-query-stories
<that commit>`. **The language stays 0.4.0.**

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`; `pnpm
format:check`. Exit codes explicit (`; echo exit=$?`), non-zero test counts.

**Look**: run Storybook from your worktree on port 6011 (`pnpm --filter @molstar/storybook exec storybook dev -p 6011
--ci --no-open`, background; stop it after) and try the query box in light and dark; say what you saw.

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE
file for more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, year range
to 2026, that author added if missing. No mention of `.seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks; the stories and what each shows; forms the design names that are not built; what you saw.
