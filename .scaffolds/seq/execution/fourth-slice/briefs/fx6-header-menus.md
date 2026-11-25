# Lane FX6: the header menus (the author's look, 2026-11-21)

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first; FX5 merged, `fd3b13b01`).
Branch `seq/fx-6`. Model: strong (the cause is likely in the viewer's input routing). The fix lane's rules
([fx-fix-lane.md](fx-fix-lane.md)): one commit an item, a test failing on the parent for each, headers updated, hot
files never staged, the look in light and dark (port 6013) before the report.

The author, on `Examples/Services` (and any story with track headers): open a track header's `⋯` menu (or the palette
slot's menu) over the view.

1. **The menu's actions do nothing**: Hide track, Move up, Reverse, Collapse… — clicking an item does not apply it.
2. **A click on a menu item goes to the canvas under it**: the click selects/acts on the data beneath the menu.
3. **The pointer over the menu still hovers the view under it**: the hover card shows and the hover emphasis is drawn
   for what is under the menu (a residue, a range item with its outline).
4. **`Color ▸` opens no submenu** (and so applies nothing).

All four likely share one cause: the viewer's pointer handling (`viewer/src/mount.ts`'s listeners, the input model,
`header-menu-dom.ts`, `header-parts-dom.ts`, `behaviors/headers.ts`, `header-entries.ts`, `@molstar/ui`'s menu
component) receives events meant for the menu — a document/window-level or capture-phase listener, a menu mounted under
the canvas or with `pointer-events: none`, a pointerdown that closes or re-renders the menu before its click, or a
hit-test that ignores the overlay. Find the cause first, with evidence; fix it at its root so a menu or any overlay
(the hover card, the keys panel, legends' popovers) owns its events: the view under it gets a `leave` (hover cleared,
card gone) while the pointer is over the overlay and never its clicks. Read track-header-parts.md (menus, submenus from
a command's select param: "Commands own choices"), ux-foundations (the grammar; the hover card), and the viewer
`_spec`.

Tests (jsdom, real event sequences: pointerdown, pointerup, click, pointermove, with the targets as a browser sends
them): each item applies its command (Hide hides, Move up moves, Reverse reverses); a click on an item selects nothing
in the view; pointermove over the menu clears the view's hover and shows no card; `Color ▸` opens its submenu by click,
hover and keyboard (ArrowRight), and a choice applies the scheme. Then **verify in the real browser** at port 6013 by
clicking (not only jsdom): on `Examples/Services/Composing Services` and `Features/Track headers`, light and dark;
report what you clicked and saw.

**Checks and report**: as [fx5-author-look.md](fx5-author-look.md) lists them (bounds 230 / 830 / 1330 / 940 kB); a
report under 250 words: the cause, commits, tests failing on the parent, what you clicked in the browser, the public
surface added, the `.v6/changelogs/` text.
