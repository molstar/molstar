# Mol\* UI stage: the specification

`@molstar/ui/stage/*` (`packages/ui/src/stage/`): the app's main area as a stage. The app draws its canvas (a WebGL
canvas, an SVG, a document); the stage lays the chrome over it (a band on top, a lead, the bottom bar, the hover label,
the empty mark) in the canvas's ink. Its sizes follow the foundations (`styles §8`); the shell places it as its `main`.
Code comments cite this spec's sections (`stage §2`). The code is the truth for how.

| File                                                 | What                                                                                             |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `canvas-ink.ts`                                      | `luminance`, `themeForColor`, `parseColorAlpha`, `inkForColor`: the chrome's theme by color (§2) |
| `stage-model.ts`                                     | `StageModel` and the hover label's types (§1)                                                    |
| `stage.tsx`                                          | `Stage`: the layers and slots, what it measures, the pointer rules (§2)                          |
| `hover-label.tsx`                                    | `HoverLabel`, `HoverFollower` (§3.1)                                                             |
| `toolbar.tsx`                                        | `StageToolbar`, `collapsedRows`, `toolRow` (§3.2)                                                |
| `empty-stage.tsx`, `../components/splash/splash.tsx` | `EmptyStage`, `GlassSplash` (§4)                                                                 |
| `glass-empty-mark.tsx`                               | `GlassEmptyMark`: the Mol\* glass as an empty mark, opt in (§4)                                  |
| `../components/logo/theme-change.ts`                 | `onThemeChange`: a theme switch around an element (§2 Ink, §4)                                   |
| `placeholders.tsx`                                   | `StagePlaceholder`, `BandPlaceholder`: stand-ins for stories and tests, unpublished (§6)         |

## 1. The model

- **Every state is the model's, DOM-free** (`StageModel`): the background and its ink, `empty`, `busy`, the hover lines
  and their placement, the band (`band`, `topInset`), the size (`compact`), the empty mark's phase and `quiet`. Every
  change is a method (`setBackground`, `setEmpty`, `setBusy`, `setHover`, `showBand`, …). Why: an app, a script or an
  agent drives the stage with no document (headless: `src/_test/headless.test.ts`).
- **What the view measures is an input with a default**: `setTopInset` and `setSize` are written by the view's
  observers; headless, `topInset` is 0 and the size unknown (0: not compact). Nothing waits on a measurement.
- **Its sizes are options**: `compactWidth` (`StageCompactWidth`, 480) and `shortHeight` (`StageShortHeight`, 360);
  `compact` and `hoverAt` read them, and the view's compact layout follows `compact` (`data-compact`), not a fixed
  container query. Why: what is narrow is the app's (a dense toolbar folds sooner).
- **One way to clear the label**: `setHover([])` (an equal list set anew changes nothing).
- **The mark's wait is the model's**, on injectable timers (`timers`, the globals by default, read at each call): empty
  and idle `emptyDelayMs` (700) brings the mark up (`mark: 'in'`) and quiets the canvas; any `busy` restarts the wait (a
  load never flashes it); work alone never takes a mark down; something drawn takes it out (`'out'`, `markFadeOutMs`
  250, `markFadeOutMs`: the view's transition takes it), only if it showed (else straight to `'gone'`). A fade out
  always ends in `'gone'`, whatever comes meanwhile (emptied again, busy): the view sees the mark go, so a mark that
  plays as it shows (the splash) plays again; the next wait starts then. Why in the model: headless, the same states
  follow (with fake timers in tests).
- **`quiet`**: the canvas is faded out while the mark is up, and from the start until something is drawn (nothing of an
  empty canvas shows before the mark). A stage made with `mark: false` is never quiet.

## 2. The stage

- **Layers**: the app's canvas box (`data-msui-stage-canvas`, `data-quiet` while quiet), then one chrome layer in the
  canvas's ink: the band (`top`), the empty mark, `lead` (top left, below the band), the hover label, the bottom bar
  (`status` | `chips` | `corner`), `overlays`. The stage's background is the canvas's color.
- **Ink** (`inkForColor`): the chrome's theme is the canvas's, by its color's WCAG luminance at 0.179 (where black and
  white text contrast equally): `dark` chrome on a dark canvas. Why: the UI theme says nothing about what the app draws;
  a light toolbar on a black scene is lost. The app's color is kept as given (`--msui-stage-color`); the model parses
  hex, `rgb()`/`hsl()` (numbers or percentages, an alpha), the CSS names and `transparent`. What it cannot read
  (`oklch()`, `var()`) or a color not opaque (alpha below 1) is `resolved: false` (light, headless); the view then asks
  the CSS engine: the stage's computed background (drawn on a 1 px `OffscreenCanvas` where it is in another syntax)
  composited over each background behind it up to the first opaque one (source over in sRGB, as CSS blends; past the
  last, white). Why: never a white stage for a color the app meant, and a nearly transparent black is the page's color,
  not black. It resolves once per color, and again only on a theme switch around the stage (`onThemeChange`:
  `data-msui-theme` on it or an ancestor, `prefers-color-scheme`; not its own ink's), while the color is unresolved: a
  `var()` follows the theme, and a few style reads at a switch cost nothing between them.
- **The chrome takes no pointer; its content does** (`:where(.slot) > *`, no specificity: a component's own
  `pointer-events`, the hover label's, the status line's, wins whatever the order the styles load in): the canvas under
  it stays the app's everywhere the chrome draws nothing.
- **The band** is a slot over the canvas (a docked band is the shell's: a pane above the stage, shell §2.1; which one,
  or both, is the app's choice (the shell's `top` option), and the shell hands the overlay's view to this slot): its
  backdrop is the canvas color, fading out below its content over `--msui-stage-band-fade` (24), the faded part letting
  the pointer through. The chrome below it starts at `--msui-stage-top-inset`, the band's content height measured by one
  `ResizeObserver` into `topInset` (the app frames its scene below it). Hidden (`showBand(false)`): paused and inert
  (`activity`), `topInset` 0.
- **Sizes** (`styles §8`): the chrome keeps `--msui-stage-edge` (8) from the stage's edges; the bottom bar's columns a
  unit pair apart; compact (`StageModel.compact`, `data-compact`) the status takes a row of its own above the chips and
  the corner.
- **Pointer**: a pinch (ctrl + wheel) is held (non-passive, `preventDefault`) only over the canvas box and the bottom
  bar, so the page does not zoom there and a wheel elsewhere (the band, a panel) stays passive; Safari's `gesture*` are
  held on the root.

- **`insets`** (a value or a state: px from the bottom-left and bottom-right corners): what of the stage's bottom edge
  something over it covers (a shell's shared chrome over this pane: its port's insets, shell §2.1); the bottom bar keeps
  clear. The stage does not know what covers it.
- **`toolbarShown`** (a value or a state; default shown): whether the corner's toolbar shows is the app's, from what it
  knows: a pane's status from its layout (shell §2.1: active, hovered, focused), so several stages show one toolbar. The
  stage reads no layout. Hidden, it fades out and takes no pointer; it stays in the tab order, and the focus in it shows
  it (`:focus-within`, whatever the app does: a keyboard user never works in an invisible toolbar).

## 3. Over the canvas

### 3.1 The hover label

- Lines of text or runs (`code`, `sub`, `sup`, `small`, `em`, `strong`, a link's look), written by `textContent` only: a
  label never parses markup. **The first line is the title**, plain or runs (one way: no role on a line); a line of runs
  may be an `item` (indented).
- **No `role="tooltip"`**: nothing it describes is focused (a tooltip is referenced by `aria-describedby`; the canvas
  under the pointer has no element to reference it). It is plain text while shown, out of the accessibility tree while
  hidden, and not live (it changes at pointer rate). What is hovered is said the app's way (a status `announce`).
- **Auto** (default, `hoverAt`): the corner on a roomy stage, beside the pointer on a narrow (`compactWidth`) or short
  (`shortHeight`) one: a small pane or a window, where a corner label would cover the middle of the canvas or be
  squeezed. The corner while the size is unknown (headless). **Corner**: the stage's bottom right, above the corner's
  toolbar (`--msui-hover-label-*` on the corner), its own width (the cell is as narrow as a folded toolbar); it does not
  move. **Cursor**: beside the pointer, flipped near the edges: one read (the layer's box, the label's size) and one
  write (two offsets as custom properties, from cached `px` strings) a frame, on `frames`; the layer and its box are
  cached, read again after the stage resizes (`observeResize`), anything scrolls (a capture listener on the document) or
  the pointer comes in; 0 B a move and its frame (`_test/hover-label.quick.ts`).
- It appears at once; an emptied label stays one frame (a fast move from one item to the next swaps its text in place,
  no blink), then fades out with its last text; a held button (a drag, a rotation) hides it at once; out of the
  accessibility tree while hidden. It shows only over the canvas and `data-msui-hover-source` (the band). Its listeners
  are attached once, scoped.

### 3.2 The toolbar

- `StageToolbar` is `components/toolbar` in the `corner`: its menus open **upward**, as tall as the room above (the
  stage's top to the button, less the menu's gap and the edge: `--msui-stage-menu-room`), measured once as a menu is
  about to open (a capture listener on the bar, no layout read per frame); a row wider than the corner scrolls across
  (`ScrollArea`). Only while no menu of it is open, and never for a target inside a menu (a submenu's row is a menu
  button too, and the popover the bar's descendant): measuring from a row would shrink the open menu.
- **Folded** (`collapsed`: `auto` by `StageModel.compact`, or fixed, or a state): one button whose menu lists the rows
  (`collapsedRows`): `leadRows`, then each block's shown items, a hairline between blocks; a toggle is checkable, an
  action a row, a menu a submenu led by the item's own action (`withOwnAction`, `components/toolbar`'s, also its tap
  menu's: unless a row of the menu is it by id, never by label: two actions may read the same). The rows are built as
  the menu opens (the items' current faces and states), never while it is open (rebuilt under the pointer, its submenus
  never open). The row stays mounted, hidden, inert and paused (`activity`). `auto` without `stage`: folded while the
  element it is placed in (the closest `[data-msui-stage]`, else its parent) is narrower than `StageCompactWidth` (that
  element's width, not its own: its own changes as it folds).
- **Text and icon are the app's**: `labels` (`StageToolbarLabels`, `DefaultStageToolbarLabels` English; the row's name
  is `labels.tools`, no `aria-label` beside it), `foldedIcon` (`wrench`).
- On touch a tap opens a button's menu with its own action first (`components/toolbar`'s rule).

## 4. The empty stage and the splash

- `EmptyStage` shows while `StageModel.mark` is not `gone`: the app's mark (`Stage`'s `emptyMark`, given whether the
  stage is compact, rendered again when that changes). **None by default**: the stage is a base and names no brand
  (nothing shows; the wait and `quiet` still run). The Mol\* glass is opt in, `GlassEmptyMark` in a module the stage
  does not import (an app without it does not bundle the glass). It fades in slowly and out over `markFadeOutMs`
  (`--msui-stage-mark-out`, `data-phase`). Only the mark takes the pointer.
- `GlassSplash` is the glass wordmark large (`components/logo/glass.ts`, travel 1, the splash's gathering;
  `data-msui-splash`): one canvas sized for its pieces; it assembles each time it shows (its scope becomes active),
  leans under the pointer, shatters on a press and gathers again. It asks for frames only while something moves, never
  at rest, never while hidden (scope paused) or off screen (`IntersectionObserver`: what was to come is asked for when
  it is back). Reduced motion (followed as it changes: on, the word whole and still at once) or no 2D canvas: one still
  frame. Glass only: other fields come with their own spec.
- **The glass reads the page only at changes**: its colors (CSS probes) and pixel ratio as it shows and again on a theme
  switch (`onThemeChange`) or a ratio change (a `resolution` query for the ratio now), one frame drawn then; the
  pointer's box as numbers, read as the pointer comes in and again after a scroll or a resize. A move or a frame
  allocates nothing. Its name is `labels.mark` (`WordmarkLabels`: the word it draws).

## 5. Several stages

A stage is one canvas; several stages are the shell's (shell §2.1): stage views in the main area's region, each with its
own `StageModel` and slots, so nothing here changes for a split. The pane a person works in is the shell's `activePane`;
whether a toolbar is per stage or shared over the main area is the app's.

## 6. Placeholders

- `StagePlaceholder` (a three-chain assembly drawn space-filling, SVG) and `BandPlaceholder` (residue boxes under a
  ruler, a selection, the hovered column written once a frame) stand in for a 3D canvas and a sequence band in stories
  and tests; they feed the model's hover lines (`structureLines`, `residueLines`). The stories never load a real canvas.
  Molecular stand-ins are not the library's: not published (`exports` null, left out of `files`); stories and tests
  import them by relative path.

## Tests

`_test/`: the model and the ink in plain Node (`stage-model.test.ts`, fake timers); the views in jsdom (`stage`,
`hover-label`, `toolbar`, `placeholders`; `components/splash/_test`, `components/logo/_test`), axe over the stage's
states (empty with the mark, loaded with the toolbar, a dark band, the label at the cursor); the CSS rules a test cannot
see in jsdom read from the generated module. `test:perf`: `hover-label.quick.ts` (a cursor move and its frame, with the
layer invalidated now and then; a corner move with holds: 0 B).
