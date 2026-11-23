# Lane RX: the genome ring, an extension (`@molstar/seqx-ext-rondo`)

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first): it holds the scaffold
(`b96e73eea`: `packages/seq/ext-rondo` with its `package.json`, tsconfigs, two empty entries, the inventory row, the root
reference, the storybook's dependency, and a `biome.json` rule allowing only `@molstar/seqx-core/extension` and
`@molstar/seqx-viewer`'s root). Branch `seq/rx`. Model: strong. Beside you: **GX** (the Genomes dataset's generator and
views in `seqx-examples/src/genomes/`, but not `genomes/ring/`), **FX19** (genotype tokens), **FX20** (`test-harness/`),
**FX21** (`formatGenomic` moved to `@molstar/core`). Stay out of their files.

**Why**: [ring-extension](../../../design/ring-extension.md), agreed 2026-11-21 with fifteen decisions. Read it whole
first (Decision, What it ships, stems, ribbons, stubs, colour, the ring zoom, the square, hover and the keyboard, other
layouts, the fallback, the builder helper, the hole, Defaults, the matrix's cells, Configuration, Consequences), then
[extension-openings](../../../design/extension-openings.md) and region-reads' pairs level as you need them.

**Amended 2026-11-23** (RX's first run stopped at step 1 on two findings, both verified): **this lane is the language
lane** while it runs: `section.center` (`title` | `none`; the design's Decided 9, GC2's brief named it, it never landed),
with the language's version bump as its `_spec` rule says, the history line, catalogs, config-keys tests, the builder;
no other lane writes in `lang/`. And **`PairsReader`** is this lane's, a Seq\* core refactor of
`render/pairs-region.ts` and `render/pairs-assembly.ts`: one reader of bins at a level and items in a region for held
and region pairs, the built-in `pairs` layer reading through it with its draw lists byte-equal before and after (a
test), its allocation and timing cases unchanged. The `gate` argument to `pairsLevels` with it.

**Task**, a green commit at each step:

0. **`section.center`**, then **`PairsReader`**, as amended above, each its own commit.
1. **The core's side** (Consequences, `seqx-core` and `seqx-viewer`): check each piece listed there exists
   (`section.include`, `section.center`, `navigator.switch`, the square bound to one section) and build what is missing
   and small: the `extension` entry's exports (`PairsReader`, the pairs bins' types, `pairsLevels`, `coarsenBins`,
   `placesOf`, beside `PolarWarp`); the three optional `LayerType` hooks `describe`, `summary`, `announce`, called by the
   hover card (as a built-in card provider's input: TT's providers, merged), the ARIA summary and the cursor, the
   viewer's own text after; extension bindings scoped to a layer type's hits, winning over the grammar's default there.
   A test for each (a hook's text reaches the card, the live region and the ARIA). The default and minimal viewers grow
   by these hooks only (the design's ~0.3 kB); report the rows. A missing piece that is not small stops with a report.
2. **The package**: `rondo` (`.`: `rondo:ribbons`, `rondo:cells`, `rondo:square`, `rondo:hues`, params, config
   `rondo.*`, `rondoRing`) and `rondoViewer` (`./viewer`: `rondo:partner`, `rondo:select-partners`, the stem menu, the
   `opens` double-click binding), stems as a nested coarser level, ribbons wide by count with fans, stubs inward, the
   hole's `center: none`; ribbons on a `tracks` row too. Its `_spec/README.md` (stems, ribbons, stubs, the fallback table,
   Configuration). Imports through the extension entry only: a gap is a step-1 item, never a reach into `src`.
3. **The fallback** (EC-235): the same spec drawn without the extension (canvas and SVG), the issue once, the version
   said; `viewer.extend` and its undo redraw.
4. **The showcase data**, `packages/seq/examples/src/genomes/ring/` only: GRCh38's 23 chromosomes (Y the 24th by
   `include`) as inline lengths typed from the GRC's assembly report (no file, no fetch), the seeded links generator (the
   examples design's step 7: 12,000 at load, 60 % within with distance decay, 40 % between, a dozen hub pairs, `type`,
   `q`; held and as BEDPE, one picture), a committed 1,000-link BEDPE for the reader's tests, `syntheticBands` for
   tests; the cytobands are the recorded `cytoBandIdeo` pages already in `fixtures/genomes/` (no fetch). Randomness from
   core's `PCG`/`hash3` only; byte-identical. Its views (the genome ring through `rondoRing`; one chromosome; a pair;
   the square beside) in a `RingViews` list of their own, so GX's `views.ts` is not yours.
5. **Stories**: a Design page "Genome ring"; a story per default of the design's Defaults (1 to 10 where data allows),
   light and dark, canvas and SVG, each with a headless twin on a fake clock (bins landing, then items in chr4's ring);
   the same story with the extension removed.

**Tests** (quick lane, `_test/`, the scaffold's exports `null` already): the Consequences list whole: drawing against
the data (every bin pair in view in exactly one ribbon, stub or fan; necks additive, a 300-partner hub; stubs only where
no place), the property tests (stems a partition, DPR 1, 2 and SVG alike; items only at or under `rondo.maxRibbons`; the
twist equal to GA's rule; the closed hit equal to the drawn cover, sampled; held and region pairs the same stems), the
fallback, EC-108, 169, 171, 206, 210, 218, 225, 235 in titles; the matrix's cells.
**Performance**: 23 blocks and 2,000 links laid out under 4 ms and drawn under 8 ms (quick representative, `itWithin`);
12,000 at the fit: layout under 16 ms once, a frame under 8 ms, a hover under 2 ms (heavy); the draw and hover loops at
0 B warm (an allocation case in the perf tier). Timing and allocation bounds are gates.

**Sizes**: the package under 30 kB minified, measured by a `_measure/rondo.ts` entry you write (and report the row for
`sizes.mjs`, which is the orchestrator's); minimal and default viewers not grown past the hooks.

**Where**: a worktree on `seq/rx`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each green
step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`): list what you need. **No language
change**: one needed stops the lane with a report. Anything in `@molstar/core` or `@molstar/io` stops the lane. A file
past 800 lines is a question.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples packages/seq/ext-rondo`, each `typecheck:test` (ext-rondo's too), the storybook typecheck,
`pnpm lint`, `pnpm format:check`, both perf tiers to completion, `pnpm -s -C packages/seq/core size` (read the rows for
FAIL), `pnpm test:medium packages/seq/`, `pnpm check:workspace`, the `.seq/_checks` (features-planned with `--plan
.seq/execution/fifth-slice/plan.md`, `defaults-designed.ts`).

**Report**: under 450 words: commits; checks with numbers; the core pieces found and built; story ids; sizes (the
package, the viewers); performance numbers; ECs covered; hot-file needs (the size row); the `.v6/changelogs/` text;
designs to amend (Consequences' list: extension-openings, layouts, themes, features, navigators).
