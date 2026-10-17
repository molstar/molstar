# Accessibility and touch

Required for every view: phones (compact mode), keyboard, screen readers, colorblind-safe and high-contrast colors. Best
practice (WCAG 2.2 AA, the WAI-ARIA patterns), checked per view in tests and Storybook. The code is the truth for how.

- **Contrast from the tokens**: the accent and status colors are text too, at least 4.5:1 on a panel; a high-contrast
  theme (black and white, a yellow accent, borders) and `prefers-contrast` (`styles/theme.css`). No meaning by color
  alone.
- **Targets**: at least 24px; 44px in compact mode (`--msp-touch-target`, `views/shell/shell.module.css`); the canvas's
  own buttons 38px there (`views/shell/canvas.module.css`).
- **Visible focus**, keyboard only (`:focus-visible`).
- **Keyboard**: every view usable without a pointer. The command menu is a modal combobox (Tab stays in, focus returns
  to the opener); trees are ARIA trees (`dom/tree-keys.ts`, `TreeNavigation`: Tab reaches a tree once, arrows move,
  open, close, Home and End jump, Enter does what a click does; focus roves, the active row the one tab stop). Trees are
  virtualized, so rows are flat siblings with `aria-level`, `aria-posinset`, `aria-setsize` and `aria-expanded` (a leaf
  none), and the keyboard moves by the rows' model: a move scrolls its row into view, then focuses it; a focused row
  scrolled away leaves the focus on the tree, which takes the keys until it is back; cards and dialogs trap and restore
  focus (`dom/focus.ts`); a closed panel gives focus to its tab; menu buttons follow the menu button pattern; the rail's
  tabs are a tab list. Keys can be remapped from the keyboard in Help (a recorder takes the next chord, Escape cancels;
  the result is said politely: [commands.md](commands.md) §6).
- **Every field has a name**: a params row's control is labelled by its label's id
  (`components/params/params-form.tsx`).
- **The canvas without sight or pointer**: picking, measuring and focusing by target (select, measure, focus take
  targets); the canvas has a text alternative (what is loaded and selected, in words; `viewer/src/ui/shell.tsx`).
- **Screen readers**: one polite live region for results (`ApiResult.did`), the latest message only
  (`managers/notifications.ts`, the canvas chrome's status line); logs are not live; the agent's reply is said once when
  finished, not as it streams; the sequence's keyboard cursor is said.
- **Touch**: a long press gives what hover gives a mouse (the hover label, the sequence's tooltip;
  `viewer/src/sources/canvas-host.ts`); ordered picks take one tap at a time, or targets.
- **Reduced motion**: `prefers-reduced-motion` turns the shell's transitions and the camera's animated moves into cuts
  (`viewer/src/api/camera.ts`, `viewer/src/sources/compose-plugin.ts`).
- **Data colors**: a colorblind-safe palette (Okabe–Ito) beside Prism's and Mol\*'s classic one, chosen in Settings (the
  viewer's `_spec/viewer-spec.md` §2.8).

## Automated checks

axe-core runs in the unit tests (jsdom), WCAG 2.0 to 2.2 A and AA plus its best practices, over every view that has a
test mounting it: `expectNoAxeViolations(root, { disable })` from `@molstar/prism-ui/testing/axe` (its own export, so
stories do not load it; a development dependency, MPL-2.0). It runs in the normal `pnpm test:unit`, so CI's native job
runs it; a new view is added to the list by one mount and one call. It sees into shadow roots, so the viewer's shell is
checked as a page has it.

- **Where**: the mock shell, the components and a card (`ui/src/views/shell/_test/accessibility.test.tsx`,
  `ui/src/components/_test/accessibility.test.tsx`); the viewer's shell over a real core with a structure loaded, every
  tab, the cards, the command menu, an open row menu, the inspectors, compact mode, a trajectory and a volume, and the
  empty viewer (`viewer/src/ui/_test/accessibility.test.tsx`). The helper checks itself (`ui/src/testing/_test`): a
  check that cannot fail checks nothing.
- **Off, always** (`AxeUnavailableInJsdom`): what needs layout or painted colors, which jsdom does not have: color
  contrast, target size, scrollable regions; and the whole-page rules (title, language, landmarks, bypass), because a
  view is mounted in a bare container. A view turns off more only with a reason in a comment where it asks.
- **Guards of what jsdom cannot see**: the touch targets are read from the CSS (`touch-targets.test.ts` in ui and
  viewer: the rail's tabs, the canvas chrome, buttons and the scene tree grow to `--msp-touch-target`); the long press
  is tested with a fake canvas (`viewer/src/sources/_test/long-press.test.ts`).
- **What axe found** and was fixed: an open select's list had no name; the empty snapshots list was a listbox with
  nothing in it, and a snapshot's option held buttons (it is a list of items now: a button applies it, its buttons are
  beside it); the compact tab bar's menu buttons sat inside the tab list; a volume's iso slider had a typographic minus
  in `aria-valuemin`.

## Audit by view

"Test" is what an automated test verifies; the rest of each row is manual (screen readers, real devices) and stays so. A
cell with a dash has nothing to verify there.

| View                                      | Keyboard                                                                                                                                                                                                                  | Names and structure (axe)                                         | Announced                                                    | Touch target in compact mode                           |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------ |
| Rail, tab strip, panels                   | Test: roving tab list; a closed panel gives its tab the focus (`focus.test.tsx`, `dom/_test`)                                                                                                                             | Test: floating and compact                                        | —                                                            | Test (CSS): `.tab`, More and menu buttons 44px         |
| Canvas chrome (tools, status pill)        | Test: the canvas stack's submenus (`submenu.test.tsx`); the collapse button, the folded toolbar's menu, arrows, Escape, focus handover and bounded, scrolling menus (`collapse.test.tsx`); menu button pattern is by code | Test: with and without data, the toolbar folded and its menu open | Test: the one polite status line (`focus.test.tsx`)          | Test (CSS): tools, pill and its buttons grow           |
| Command menu                              | Test: modal combobox, focus kept and returned (`command-menu-a11y.test.tsx`)                                                                                                                                              | Test: open                                                        | Test: what a search found is said once typing pauses         | By `--msp-row-height`; manual                          |
| Home, Help, Settings, Log                 | Test: the keys recorder and remap (`help-panel.test.tsx`); settings are native fields                                                                                                                                     | Test: loaded and empty                                            | Test: a remap is said politely; the Log is `aria-live="off"` | By `--msp-control-height`; manual                      |
| Structure: scene tree, inspector, menus   | Test: ARIA tree, arrows, Home and End, Enter (`scene-tree-virtual.test.tsx`, `tree-keys.test`)                                                                                                                            | Test: node picked, row menu open                                  | —                                                            | Test (CSS): rows, expander, row buttons                |
| Selection tab and bar                     | Native buttons                                                                                                                                                                                                            | Test: with a selection                                            | The summary is a status region (see the gaps)                | Test (CSS): bar and its buttons                        |
| State tree                                | Test: ARIA tree (`state-tree-virtual.test.tsx`)                                                                                                                                                                           | Test: cell picked (form, actions)                                 | —                                                            | As the scene tree; manual                              |
| Snapshots tab, Story card                 | Test: Alt+Up and Down reorder (`snapshots-panel.test.tsx`); Enter and Space apply, by code                                                                                                                                | Test: with snapshots                                              | —                                                            | Test (CSS): row and its buttons, 24px (44px compact)   |
| History                                   | Native buttons                                                                                                                                                                                                            | Test: with steps                                                  | A turn's reply is said once when finished, by code           | By `--msp-row-height`; manual                          |
| Export, Canvas, Animation cards (dialogs) | Test: a card takes, traps and returns the focus (`focus.test.tsx`)                                                                                                                                                        | Test: open                                                        | —                                                            | By `--msp-control-height`; manual                      |
| Params forms, select lists, color picker  | Test: select arrows, Enter, Escape (`select.test.tsx`)                                                                                                                                                                    | Test: every control kind, lists open                              | —                                                            | Test (CSS): `.small` and `.tiny` buttons; see the gaps |
| Context menu                              | Test: Escape and a press outside (`context-menu.test.tsx`); the arrows are by code                                                                                                                                        | Test                                                              | —                                                            | By `--msp-row-height`; manual                          |
| Volume levels and curve                   | Test: markers are sliders stepped by the arrows (`volume-levels.test.tsx`)                                                                                                                                                | Test: iso slider                                                  | `aria-valuetext` in σ                                        | See the gaps                                           |
| Trajectory bar                            | Native controls                                                                                                                                                                                                           | Test: with an ensemble                                            | —                                                            | Manual                                                 |
| Sequence band                             | Seq\* owns its keyboard (`packages/sequences`); not verified here                                                                                                                                                         | Test: band mounted over a structure                               | The cursor is said politely, by its package                  | Seq\*'s own; manual                                    |
| The canvas itself                         | Targets by text: select, measure and focus take targets, no pointer needed                                                                                                                                                | Image with a text alternative (the shell)                         | Results through the one polite region                        | Test: a long press gives the hover label (fake canvas) |

**Manual, not covered**: a screen reader pass of each row (NVDA with Firefox, VoiceOver with Safari, TalkBack), because
jsdom has no accessibility tree; contrast of every theme and of the canvas ink over the molecule (axe's `color-contrast`
needs a browser); the rendered size of every target (axe's `target-size` needs layout); the Prism assistant's panels
(not mounted in these tests); reading order with real layout; a real touch device (long press timing against the
browser's own context menu, pinch and drag on the canvas, the bottom sheet with the on-screen keyboard).

**Gaps found, not yet fixed**:

- Controls smaller than 24px that are not growing with `--msp-touch-target`: a CSS scan finds fixed 16 to 23px sizes in
  the color picker, params form rows (20px and 22px buttons), the playback deck, the audio bar, quick controls, the
  selection lists, measurements, the sequence strip's buttons and the volume transfer. Which of them are buttons, and
  which are exempt (inline, spaced by 24px, a mouse-only drag handle), is decided with `target-size` run in a real
  browser (a Storybook test run or Playwright), not by reading the CSS.
- Reordering by dragging (Snapshots; other `dragToReorder` lists) has the keyboard (Alt+Up and Down) but no
  single-pointer alternative for a touch user (WCAG 2.5.7): a Move up and Move down in the row's menu.
- The Selection summary is a status region (implicitly polite) beside the one live region for results: a selection
  changing by dragging may speak twice; it is to be heard with a screen reader before deciding.
- A long press shows the hover label until the next touch; nothing dismisses it earlier, and a screen reader user gets
  nothing from it (the canvas text and the targets are their way).
