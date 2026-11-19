# Lane PR: a drag ends wherever the button is released

**The bug** (the author, 2026-11-19: "when dragging the slider, if you release outside the canvas, it doesn't release
the event. this needs to be handled by creating a temporary event handler on the window/document. other controls
prolly have the same issue"): the mounts (`packages/seq/viewer/src/mount.ts`, `mount-navigator.ts`) end a drag only on
the container's `pointerup`, `pointercancel` and (in `mount.ts` only) `lostpointercapture`, relying on pointer capture;
where capture is refused or lost (an iframe boundary as in Storybook, a pointer already up), a release outside leaves
the drag running.

**Fix** ([navigators](../../../design/navigators.md), "A drag ends wherever the button is released"): every drag a
mount starts adds temporary `pointerup` and `pointercancel` listeners on the element's own window
(`ownerDocument.defaultView`), removed when the drag ends; also end it on `lostpointercapture` (the navigator mount
lacks it), a `pointermove` with `buttons === 0`, the window's `blur`, and `visibilitychange` to hidden. Each path ends
the drag once (`ViewerInput.up` and `NavigatorModel.up`/`cancel` idempotent; check). One small shared helper in the
viewer (a drag session: start on a taken press, end on any of these, dispose its listeners) used by both mounts and
by every storybook story that runs its own drag (`selection-root.ts`, `story.stories.ts`, any other found by searching
for `pointerdown`). Keep capture as now. No listener left behind after a drag or an unmount.

**Tests** (jsdom, fail on the parent): per mount, a drag released by a `pointerup` dispatched on the window, on another
element, by a buttonless move, by window blur, by a hidden document; listeners removed after the drag and on unmount;
a second press after a lost release starts a new drag cleanly.

**Not in scope**: anything else. HL, FS-A2 and RC run beside you.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/pr-pointer-release
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/viewer`;
`pnpm -C packages/seq/viewer typecheck:test`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`;
`pnpm format:check`. Exit codes explicit (`; echo exit=$?`), non-zero test counts.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.seq/` in code; the viewer `_spec`'s Input
says the rule.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks; how to try it (a navigator drag released outside the Storybook iframe).
