# Mol\* UI stage: the specification

`@molstar/ui/stage/*` (`packages/ui/src/stage/`): the app's main area as a stage. The app draws its canvas (a WebGL
canvas, an SVG, a document); the stage lays the chrome over it (a band on top, a lead, the bottom bar, the hover label,
the empty mark) in the canvas's ink. Its sizes follow the foundations (`styles §8`); the shell places it as its `main`.
Code comments cite this spec's sections (`stage §2`). The code is the truth for how.

| File                                                 | What                                                                                        |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| `canvas-ink.ts`                                      | `luminance`, `themeForColor`, `parseColor`, `inkForColor`: the chrome's theme by color (§2) |
| `stage-model.ts`                                     | `StageModel` and the hover label's types (§1)                                               |
| `stage.tsx`                                          | `Stage`: the layers and slots, what it measures, the pointer rules (§2)                     |
| `hover-label.tsx`                                    | `HoverLabel`, `HoverFollower` (§3.1)                                                        |
| `toolbar.tsx`                                        | `StageToolbar`, `collapsedRows` (§3.2)                                                      |
| `empty-stage.tsx`, `../components/splash/splash.tsx` | `EmptyStage`, `Splash` (§4)                                                                 |
| `placeholders.tsx`                                   | `StagePlaceholder`, `BandPlaceholder`: stand-ins for stories and tests (§6)                 |

## 1. The model

- **Every state is the model's, DOM-free** (`StageModel`): the background and its ink, `empty`, `busy`, the hover lines
  and their placement, the band (`band`, `topInset`), the size (`compact`), the empty mark's phase and `quiet`. Every
  change is a method (`setBackground`, `setEmpty`, `setBusy`, `setHover`, `showBand`, …). Why: an app, a script or an
  agent drives the stage with no document (headless: `src/_test/headless.test.ts`).
- **What the view measures is an input with a default**: `setTopInset` and `setSize` are written by the view's
  observers; headless, `topInset` is 0 and the size unknown (0: not compact). Nothing waits on a measurement.
- **The mark's wait is the model's**, on injectable timers (`timers`, the globals by default, read at each call): empty
  and idle `emptyDelayMs` (700) brings the mark up (`mark: 'in'`) and quiets the canvas; any `busy` restarts the wait (a
  load never flashes it); work alone never takes a mark down; something drawn takes it out (`'out'`, `markFadeOutMs`
  250), only if it showed (else straight to `'gone'`). Why in the model: headless, the same states follow (with fake
  timers in tests).
- **`quiet`**: the canvas is faded out while the mark is up, and from the start until something is drawn (nothing of an
  empty canvas shows before the mark). A stage made with `mark: false` is never quiet.

## 2. The stage

- **Layers**: the app's canvas box (`data-msui-stage-canvas`, `data-quiet` while quiet), then one chrome layer in the
  canvas's ink: the band (`top`), the empty mark, `lead` (top left, below the band), the hover label, the bottom bar
  (`status` | `chips` | `corner`), `overlays`. The stage's background is the canvas's color.
- **Ink** (`inkForColor`): the chrome's theme is the canvas's, by its color's WCAG luminance at 0.179 (where black and
  white text contrast equally): `dark` chrome on a dark canvas. A color that does not parse counts as white. Why: the UI
  theme says nothing about what the app draws; a light toolbar on a black scene is lost.
- **The chrome takes no pointer; its content does** (`.slot > *`): the canvas under it stays the app's everywhere the
  chrome draws nothing.
- **The band** is a slot over the canvas (a docked band is the shell's: a pane above the stage, shell §2.1; which one,
  or both, is the app's choice (the shell's `top` option), and the shell hands the overlay's view to this slot): its
  backdrop is the canvas color, fading out below its content over `--msui-stage-band-fade` (24), the faded part letting
  the pointer through. The chrome below it starts at `--msui-stage-top-inset`, the band's content height measured by one
  `ResizeObserver` into `topInset` (the app frames its scene below it). Hidden (`showBand(false)`): paused and inert
  (`activity`), `topInset` 0.
- **Sizes** (`styles §8`): the chrome keeps `--msui-stage-edge` (8) from the stage's edges; the bottom bar's columns a
  unit pair apart; below `StageCompactWidth` (480, a container query on the stage and `StageModel.compact`) the status
  takes a row of its own above the chips and the corner.
- **Pointer**: a pinch (ctrl + wheel) is held (non-passive, `preventDefault`) only over the canvas box and the bottom
  bar, so the page does not zoom there and a wheel elsewhere (the band, a panel) stays passive; Safari's `gesture*` are
  held on the root.

- **`insets`** (a value or a state: px from the bottom-left and bottom-right corners): what of the stage's bottom edge
  something over it covers (a shell's shared chrome over this pane: its port's insets, shell §2.1); the bottom bar keeps
  clear. The stage does not know what covers it.
- **`toolbarShown`** (a value or a state; default shown): whether the corner's toolbar shows is the app's, from what it
  knows: a pane's status from its layout (shell §2.1: active, hovered, focused), so several stages show one toolbar. The
  stage reads no layout. Hidden, it fades out and takes no pointer; it stays in the tab order (the focus coming in makes
  its pane focused, and the app shows it).

## 3. Over the canvas

### 3.1 The hover label

- Lines of text or runs (`code`, `sub`, `sup`, `small`, `em`, `strong`, a link's look), written by `textContent` only: a
  label never parses markup.
- **Auto** (default, `hoverAt`): the corner on a roomy stage, beside the pointer on a narrow (`StageCompactWidth`, 480)
  or short (`StageShortHeight`, 360) one: a small pane or a window, where a corner label would cover the middle of the
  canvas or be squeezed. The corner while the size is unknown (headless). **Corner**: the stage's bottom right, above
  the corner's toolbar (`--msui-hover-label-*` on the corner), its own width (the cell is as narrow as a folded
  toolbar); it does not move. **Cursor**: beside the pointer, flipped near the edges: one read (the layer's box, the
  label's size) and one write (two offsets as custom properties, from cached `px` strings) a frame, on `frames`; 0 B a
  move (`_test/hover-label.quick.ts`).
- It appears at once; an emptied label stays one frame (a fast move from one item to the next swaps its text in place,
  no blink), then fades out with its last text; a held button (a drag, a rotation) hides it at once; out of the
  accessibility tree while hidden. It shows only over the canvas and `data-msui-hover-source` (the band). Its listeners
  are attached once, scoped.

### 3.2 The toolbar

- `StageToolbar` is `components/toolbar` in the `corner`: its menus open **upward**, as tall as the room above (the
  stage's top to the button, less the menu's gap and the edge: `--msui-stage-menu-room`), measured once as a menu is
  about to open (a capture listener on the bar, no layout read per frame); a row wider than the corner scrolls across
  (`ScrollArea`).
- **Folded** (`collapsed`: `auto` by `StageModel.compact`, or fixed, or a state): one button whose menu lists the rows
  (`collapsedRows`): `leadRows`, then each block's shown items, a hairline between blocks; a toggle is checkable, an
  action a row, a menu a submenu led by the item's own action. The rows are built as the menu opens (the items' current
  faces and states). The row stays mounted, hidden, inert and paused (`activity`).
- On touch a tap opens a button's menu with its own action first (`components/toolbar`'s rule).

## 4. The empty stage and the splash

- `EmptyStage` shows while `StageModel.mark` is not `gone`: the app's mark (`Stage`'s `emptyMark`, given whether the
  stage is compact; default `DefaultEmptyMark`, the glass splash; `() => null` for none: the stage is a base, the mark
  is the app's); it fades in slowly and out in 250 ms (`data-phase`). Only the word takes the pointer.
- `Splash` is the glass wordmark large (`components/logo/glass.ts`, travel 1, the splash's gathering): one canvas sized
  for its pieces; it assembles each time it shows (its scope becomes active), leans under the pointer, shatters on a
  press and gathers again. It asks for frames only while something moves, never at rest, never while hidden (scope
  paused) or off screen (`IntersectionObserver`: what was to come is asked for when it is back). Reduced motion or no 2D
  canvas: one still frame. Glass only: other fields come with their own spec.

## 5. Several stages

A stage is one canvas; several stages are the shell's (shell §2.1): stage views in the main area's region, each with its
own `StageModel` and slots, so nothing here changes for a split. The pane a person works in is the shell's `activePane`;
whether a toolbar is per stage or shared over the main area is the app's.

## 6. Placeholders

- `StagePlaceholder` (a three-chain assembly drawn space-filling, SVG) and `BandPlaceholder` (residue boxes under a
  ruler, a selection, the hovered column written once a frame) stand in for a 3D canvas and a sequence band in stories
  and tests; they feed the model's hover lines (`structureLines`, `residueLines`). The stories never load a real canvas.

## Tests

`_test/`: the model and the ink in plain Node (`stage-model.test.ts`, fake timers); the views in jsdom (`stage`,
`hover-label`, `toolbar`, `placeholders`; `components/splash/_test`), axe over the stage's states (empty with the mark,
loaded with the toolbar, a dark band, the label at the cursor). `test:perf`: `hover-label.quick.ts` (a cursor move and
its frame, a corner move with holds: 0 B).
