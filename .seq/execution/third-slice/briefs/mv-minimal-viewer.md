# Lane MV: the minimal viewer, injected and opt in

**Goal** (the author, 2026-11-19: "Injectable, opt-in"): build
[viewer-composition.md](../../../design/viewer-composition.md), Creating a viewer, "Injected, opt in": the renderer's
layer types and layouts from the registry the spec lists, never a built-in table `Renderer` defaults to; the viewer's
history, story, navigator and rows pieces as behaviors or parts the default spec lists, never imported by
`createSeqViewer`; the `render` and `export` barrels re-exporting nothing a minimal composition does not ask for.
`DefaultSeqViewerSpec()` and `createSeqContext`'s default list everything: every existing story, example and test
composes as before, unchanged. What SW found keeps the minimal viewer large (plan, stage 7, M3): `BuiltInLayerTypes`
and the layouts table as `Renderer`'s default, `render/navigator.ts` via the `render` barrel and `export.ts`, the
viewer's static imports of history, story, navigator and rows.

**The agreed surface**: the registry kinds the design names (layer types, layouts) and where the default lists live;
`sideEffects: false` in `packages/seq/{core,viewer}/package.json` is a hot file: ask for it in the report. A breaking
change to an exported name: list it for `.v6/plans/breaking-v6-changes.md` in the report (Seq\* is unreleased 0.x: the
orchestrator records it). Anything past this: stop and ask.

**Tests** (fail on the parent): a viewer composed from an empty registry draws nothing it was not given and says so
(an issue for an unknown layer type, as an unregistered format does); the default spec draws every example as before
(the existing suite, unchanged); the size script's `viewer, minimal` entry's esbuild metafile holds none of the
unlisted modules (navigator, story, history, rows, the records runtime), as a test or a check the size script runs.

**Two consumers**: the canvas and `renderSvg` take the same registry's layer types; one test composes both from a
minimal registry.

**Sizes**: `pnpm -C packages/seq/core size` before and after, every entry; the minimal viewer's drop the measure of
success (SW: 959 against 1023 kB; after HC 1010.9 against 1085.5 kB).

**Model**: strong (opus). **Base**: `packages-seq` after RK merges (the orchestrator gives the commit). Storybook port
6021 (your own; stop it after).

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, with `pnpm build` and `pnpm check:workspace` at the end (the
known `glass-wordmark` failure is master's).

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; sizes before and after; the breaking names; drift; followups; a `.v6/changelogs/` entry's
text.
