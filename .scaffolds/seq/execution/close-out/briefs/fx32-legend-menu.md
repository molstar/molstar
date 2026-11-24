# Lane FX32: the legend's layers menu styled as the header's, with Show all and Hide all

**Model**: Sonnet (one menu, decided). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard
packages-seq` first). Branch `seq/fx32`. No lane beside you.

**Why**: the use gate (the author, 2026-11-24, a screenshot of a legend row's eye menu over a merged row, "cdk2" and
"cyclin" each checked): "one weird thing i notices is the styling on the legend buttons. maybe tweak this to make it
consistent (and perhaps add "hide/show all"), then can consider the use gate done". In the screenshot each item is a
separate pill with a dark border, unlike the track header's Layers menu (`@molstar/ui`'s list rows).

**The decision** (write it into the designs and specs as written):
- The legend's layers menu looks and behaves as the track header's Layers menu: the same component, the same rows, in
  light and dark. A likely cause: `mount-legend.ts` opens `openContextMenu` on `doc.body`, the header
  (`header-menu-dom.ts`) on the viewer's container, inside the styles and theme it reads; find and fix the real cause.
  The legend's eye buttons themselves match the header's buttons too (size, hover, focus).
- A menu of several layers (the legend's and the header's Layers menu alike) ends with a separator and two entries:
  **Show all** and **Hide all**, each disabled with its reason when it would change nothing ("Every layer is shown",
  "Every layer is hidden"); each is one patch, one undo step, over the layers it lists. A track's last shown layer
  hidden by Hide all follows the rule a single hide follows (the row stays, the header says so).

**Read first**: the viewer spec's Legends and Track headers (the Layers menu); `.scaffolds/seq/design/themes.md`
(Legends) and `track-header-parts.md` (the menu); `viewer/src/mount-legend.ts`, `legend-eye-dom.ts`,
`header-menu-dom.ts`, `header-entries.ts` (the Layers entries); `@molstar/ui`'s `context-menu` (its items: `disabled`,
groups or a separator); how a patch of several edits is one undo step (`commands.ts`, the history).

**Items**, a commit at each:

1. **The look**: the cause found and fixed; the legend's menu and eye the header's. Proof: a DOM test that the
   legend's menu renders in the same container and with the same row classes as the header's, failing on its parent;
   the orchestrator's look in light and dark.
2. **Show all, Hide all**: in both menus, disabled with a reason when idle, one patch and one undo step. Proof: tests
   (each entry, disabled states, one undo step restoring every layer, the keyboard reaching them).
3. **The records**: `themes.md` (Legends) and `track-header-parts.md` (the Layers menu) and their Storybook pages, the
   viewer spec's Legends and Track headers. Proof: `designs-in-specs` green.

**How**: host neutral; nothing per frame; module size (the legend's menu stays in its lazy chunk).

**Where**: a worktree on `seq/fx32`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` (load under 3: `uptime` first);
`pnpm -s -C packages/seq/core size` (every row: spec 265.0 of 268, headless 995.0 of 1,000, default 1,587.5 of 1,592,
minimal 1,119.6 of 1,124); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold checks in
`.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 300 words): the cause; commits; each item with its proof (failing on its parent); the checks with
exit codes and counts; the size rows; which story to look at; the changelog text.
