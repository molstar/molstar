# Lane CF: the chromosome's core findings

**Goal**: what CS found in `@molstar/seqx-core` (and one base-library addition) at 83 Mb, each with a test that fails
on its parent:

1. **A ruler's units** ([axis.md](../../../design/axis.md), "A ruler's units", decided): the ruler layer's `units`
   param, `none` (default) or `bases` (bp, kb, Mb by the tick step), in the language (`lang/language.ts`, the ruler's
   params beside `numbering`), resolve and the ruler's labels (`render/ruler.ts:226-234`, today `String(n)`); the
   formatting in `@molstar/core/util/number` as [base-library.md](../../../design/base-library.md) agrees ("genomic unit
   formatting and tick steps"): that is **the agreed surface**; nothing else in `@molstar/core`. Both consumers: canvas
   and SVG labels equal. The chromosome view (`packages/seq/examples/src/p53/chromosome.ts`) writes `units: 'bases'`:
   the one line in examples you may touch.
2. **Inline data tiles the whole system at level 1** (`resolve/model.ts:167` documents "one tile (the whole system)",
   but an inline node becomes a `seq` source of levels `[1]`, `seq-format.ts:281`, `resolve/model.ts:551-558`): an
   inline node over a long system costs 81,307 level-1 tile keys a frame at the fit (`store.ts:404-409`,
   `layer-data.ts:189`). Make the code do what the doc says (one tile, or the levels a small inline payload needs),
   and make `DataStore.keysOf`/`readTiles` stop past `maxTilesPerWant` instead of building every key (the
   `seq.store.too_many_tiles` issue stays). Test: an inline node over an 83 Mb system read at the fit costs a bounded
   number of keys.
3. **Shown means in view** ([dynamic-nodes.md](../../../design/dynamic-nodes.md), Lifetime, decided): a template's
   instance holds its node only while the node's region meets the view and its margin (`render/tracks.ts:1531`,
   `nodes.holdShown`; `view/templates.ts:49-64,99-125`); out of it the node is evicted as usual and its track goes.
   Test: browse five genes' exons and return to the fit: no exon tracks remain once their regions are out of view and
   evicted; a track whose node is on screen never vanishes (the reason for `holdShown`).

**Not in scope**: the hover card's strand (viewer, after RK); the ranges layout walking every row (core TODO; CS
answers density instead). RK (keyboard cursor) works in `packages/seq/viewer` and in `core/src/render/{rows,tracks,
grid,layer-type}.ts`: keep your `tracks.ts` change to the held-nodes hunk; CS works in `packages/seq/examples` and
`storybook`.

**Performance representative**: the chromosome's fit after browsing, frames under 16 ms (performance.md), as a heavy
`itWithin` with a quick representative (bounded keys a frame); `pnpm -C packages/seq/core test:perf` unchanged or
better.

**Model**: strong (opus). **Base**: `packages-seq` at the commit the orchestrator gives. Storybook port 6022.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, plus `pnpm test:unit packages/core/src/util` and `npx tsc -b
packages/core` for the base-library addition, and at the end `pnpm build` and `pnpm check:workspace` (the known
`glass-wordmark` failure is master's).

**Headers**, **hot files** and **stop when** (anything past the agreed surface, or a language change past `units`): as
every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; sizes; the perf numbers; drift; followups; a `.v6/changelogs/` entry's text (and the
`@molstar/core` addition's).
