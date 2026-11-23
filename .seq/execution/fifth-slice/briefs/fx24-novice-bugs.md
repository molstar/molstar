# Lane FX24: the novice pass's bugs

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx24`. Model:
Sonnet. Nothing beside you.

**Why**: the fifth slice's novice pass ([novice-findings.md](../novice-findings.md), items 1–4), and a flaky test found
at stage 4's gate ([fixes.md](../fixes.md), "The hover card story test failed once").

**Measure first, then fix, each with a test failing on the parent**:

1. **The ring's card shows a wrong position.** On Features › Genome Ring › Kinds of links
   (`seq-next-features-genome-ring--kinds-of-links`), the card for the t(1;17) ribbon ends with "chr1 58,720,257", then
   "q32.2 chr1 207,100,001–211,300,000" twice. On a chr14 contact stem it shows "chr14 58,720,257" with the same chr1
   band. Find where the card's position rows and band rows are built for a ring (the default card provider's position,
   the cytobands track's row) and say what they read. The likely causes are a stale or first block, an angle mapped to
   the wrong block, or the bands looked up at a cached position. The rows must name the block, position and band under
   the pointer, and the band row must appear once. Test headless: a hover at a known angle on two different blocks gives
   each its own block, position and band.
2. **The square's moves are missing from the view's history.** On The square beside
   (`seq-next-features-genome-ring--the-square-beside`), a click on chr4's cell moves the window, but Back (⌥←) stays
   `aria-disabled`. A drag or a resize does the same. A move by the square over blocks is one history step: a click, a
   drag (on release, not each frame) and a resize alike, as `navigator.ts` already makes a click on a block's context
   ("a jump to that block"). Back then returns to the previous window. Tests: after each gesture `canBack` is true, and
   Back restores the window. A drag's frames add one step, not many, and allocate nothing new (the drag's perf case
   stays 0 B warm).
3. **Check, then fix if real**:
   - a click on a menu item that has a submenu ("Show with a partner" in the stem menu) doesn't open the submenu, while
     the Right arrow does. Test it with real pointer events, headless, through the menu's model;
   - the wheel over the ring never scrolls the page. Say what the wheel does there (zoom or turn) and whether a wheel at
     the end of its range passes to the page as the linear views' does. Bring it in line, or say why not.
4. **The flaky test**: `packages/seq/storybook/src/_stories/features/hover-card/_test/hover-card-stories.test.ts`, "a
   host's adapter is told the viewer's card and place on rest", failed once in a quick run over `packages/seq/`,
   `cli/seqx-render/` and `packages/core/src/util` at load about 3. The adapter had been called, but
   `HoverCardProvider.shown.value` was undefined. Find what still reads the real clock under fake timers, or what state
   another file leaves behind (`isolate: false`), and fix the cause; a retry is not a fix. Show that it stays green
   under load: run it 20 times beside a busy process, and report how.

**Where**: a worktree on `seq/fx24`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); a commit per item;
hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). No language change. Anything in
`@molstar/core` or `@molstar/io` stops the lane. Never touch the network. Never overwrite files outside your worktree.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples packages/seq/ext-rondo`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers to completion, `pnpm -s -C packages/seq/core size` (read the rows), `pnpm
test:medium packages/seq/`.

**Report**: under 250 words: what each item measured (the cause), the fixes, commits, checks with numbers, the story ids
to look at, and the `.v6/changelogs/` text.
