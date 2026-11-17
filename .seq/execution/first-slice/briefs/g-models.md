# Lane G: the viewer's models

**Goal**: in `packages/seq/viewer/src/`: the DOM-free models of [viewer-models.md](../../../design/viewer-models.md) for
the first slice, on `@molstar/core/reactive`, tested in Node: `SpecModel` (the spec's immutable root, its revision as
the root's hash, a preview overlay; a `HistoryPort` whose checkpoint is the root), `ViewModel` (the view shown, a
viewport per section, scroll groups on equal axes, follow through a resolved view's follow path at scale 1, 3 or 1/3,
followers never lead, viewports surviving updates EC-111, intents before a width EC-112), `HoverModel` (`set(hit,
owner)` and `clear(owner)`, one owner at a time EC-141, no allocation per move; the hover's Location projected onto
linked sections lazily through `seqx-core/resolve` paths with closure and the copy rule), `StoreModel` (owning a
`seqx-core/data` store for the viewer), the `hover` and `navigate` behaviors as small policies over the models, the
gesture-binding table (`{ gesture, modifiers, target, action }`, the default preset of ux-foundations section 8,
resolving a gesture by hit target), and the commands' change-kind table. Time and frames through `@molstar/ui/dom`'s
`DomHost`, patched in tests. Where a behavior could reasonably go more than one way, make it a config item with the
recommended default.

**Read first**: `.seq/design/viewer-models.md`, `ux-foundations.md` (section 8), `changes.md`, `viewer-composition.md`
(only what models need), `locations-and-queries.md` (Location), `packages/core/src/reactive/_spec/README.md`,
`packages/ui/src/history/` and `packages/ui/src/commands/` (types), `packages/ui/src/dom/host.ts`, the core `_spec`
(Resolve, Resolved view, Data, Queries), `.seq/design/testing.md`.

**Model**: strong.

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/g-models`; commit on it only.

**Touches**: `packages/seq/viewer/src/**` (no Canvas2D backend: lane H), its `_test/`, and
`packages/seq/viewer/_spec/README.md` (sections for what you build). Read-only: `packages/seq/core`. Not in scope:
draw lists (lane F), `createSeqViewer` and ports (lane H), selection beyond what hover needs.

**Tests** (Node, no DOM, `DomHost` patched): viewport rules EC-111 to EC-114; hover ownership EC-141; projection of a
protein hover onto the gene through the codon map (a nucleotide lifting to its codon, both strands); bindings resolving
by target; undo through the history restoring the spec's root; `hover.set` reusing its records once warm.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`,
`vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, `packages/seq/*/package.json` (ask in the report for
the entry exports and `_test` nulls you need), anything in `.seq/` and in `packages/core/`. Ask in the report for
what you need there.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>`. No mention of `packages/sequences` or `.seq/` in code,
comments or tests.

**Checks**, each its own command, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/viewer` (your folder under 2 s, a test under 100 ms)
- `npx tsc -b packages/seq/viewer`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a permission is refused (report the
action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the quick lane's time for your folder; the
exports and `_test` nulls you need in `packages/seq/core/package.json`; any place the designs were unclear or wrong (as
drift: what, where, options, your recommendation); followups.
