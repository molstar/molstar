# Lane FX7: six open items of fixes.md, clear of RG's files

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-7`. Model:
Sonnet. Follow `fx-fix-lane.md` (one commit per item, a test failing on its parent, a look in light and dark). Each
item's full text, the author's words and its tests are in `fixes.md` (Open) under the heading named here.

**Beside you**: RG's branch (not merged) changes `core/src/data/*` (store, tables, summaries, keys, protocol),
`render/{pairs*,grid*,matrix,layer-data,layer-type,tracks*,renderer,built-ins,config,index,issues,region-hover}.ts`,
`render/layout/{rings,warp}.ts`, `resolve/{model,model-region,unread}.ts`, `view/view.ts`, `viewer/src/describe*.ts`,
both `_spec/README.md`s. **Touch none of these**; an item that needs one stops and says so. A design draft lane runs
too (no code).

**Items**:

1. **"A long deletion painted as an SNV"**: the majority by width (`render/ranges.ts`, FX6's rule).
2. **"The VCF ClinVar story test fails under load"**: deterministic on a fake clock; 20 runs under load.
3. **"A click on one piece of a feature across a cyclic origin selects only that piece"** (both plasmid views).
4. **"A follower's rail turns a window over the cyclic origin into its complement"**.
5. **"A tree opened from its toggle has no way back on screen"**.
6. **"The blocks' Contacts track leaves ATP out"**.

Out of scope: the class items (hover contract, legends, empty well, drag), anything needing a design, `@molstar/core`,
`@molstar/io`, the language.

**Where**: a worktree on `seq/fx-7`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`.seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list edits. Stop any
Storybook server you start.

**Checks** (with numbers): `pnpm test:unit packages/seq/` (under 2 s a touched folder), `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers, `pnpm -s -C packages/seq/core size`, the `.seq/_checks` (features-planned with
`--plan .seq/execution/fifth-slice/plan.md`). No `pnpm test:heavy`.

**Report**: under 250 words: commits per item; checks; story ids; items stopped and why; the `.v6/changelogs/` text.
