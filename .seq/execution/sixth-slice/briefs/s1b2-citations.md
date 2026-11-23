# Lane S1E: comments cite specs, not designs (stage 1b)

**Model**: Sonnet (mechanical). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq`
first). Branch `seq/s1e`.

**Why** (`../plan.md`, stage 1b): code, comments, tests and specs never cite the scaffold. About 160 comments and test
headers in `packages/seq` cite a design by its file name ("molecules.md, Links"; "(sources.md, Consequences)"); five
test headers cite `.seq/sources/seqstar/edge-cases.md`; the core spec cites design files itself (around line 587);
`data/_test/bounds.test.ts` cites "Data, Bounds", a lead now named "An answer's bounds".

**Task**, a commit per package:

1. Find every citation: `git grep -n -E "\b[a-z0-9-]+\.md\b" -- packages/seq` and keep only those naming a file in
   `.seq/design/` or `.seq/sources/` (a package's own spec files, e.g. core's `util/_spec/assets.md` or a
   `_spec/README.md`, stay; `TODO.md`, `AGENTS.md`, `CHANGELOG.md` stay). Also `git grep -n "\.seq/" -- packages`
   except `packages/seq/AGENTS.md` (its one pointer stays until stage 3).
2. For each, read that design's `## In the specs` table: the row for the decision cited names the spec section. Rewrite
   the citation to it in the form the code already uses for spec citations: `(_spec/README.md, <Section>)` for the
   package's own spec, `(seqx-core _spec, <Section>)` for another package's. A citation whose row is design-only is
   dropped (keep the sentence, drop the parenthesis), unless the comment needs a reason, then say the reason in words.
3. An `EC-<n>` id stays (it names the edge-case catalog's item; `ec-ids.ts` reads it); only the path or file name goes.
4. Fix `data/_test/bounds.test.ts`'s "Data, Bounds" to the lead's new name, and the core spec's own design-file
   citations (rewrite to its own sections, or drop).
5. Add a check, `.seq/_checks/no-design-citations.ts` (in the shape of `defaults-designed.ts`): fails when a file under
   `packages/` names a design or source file of the scaffold (the list read from `.seq/design/` and `.seq/sources/` at
   run time), or `.seq/` outside `packages/seq/AGENTS.md`. It passes on your branch.

**Where**: a worktree on `seq/s1e`. Comments and spec text only: no code behaviour changes (a test's name may change
only if it holds a design file name; say which). Headers: every `.ts` file you modify gets
`@author David Sehnal <david.sehnal@gmail.com>` if missing and its copyright year range extended to 2026. Hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`), never `state.md`, `plan.md`, `loop.md`. Write only in your
worktree and your own temp space. No network.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/ packages/core/src/repr`, `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples packages/seq/ext-rondo`, `pnpm lint`, `pnpm format:check`, and after
`git add .seq packages`: `node .seq/_checks/no-design-citations.ts`, `node .seq/_checks/designs-in-specs.ts`,
`node .seq/oz/_checks/merge-gate.ts --package .seq/oz --instance .seq`.

**Report** (under 250 words): commits; citations rewritten, dropped, kept (counts per package); the check's result;
any citation you could not place.
