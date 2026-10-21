# Styling

How prism-ui is styled and the rules every Prism view follows (the look, the spacing, performance). The code is the
truth for how: `src/styles/` (`theme.css`, `reset.css`, `registry.ts`) and `scripts/css.mjs`.

## 1. Decisions

- **CSS modules**, one next to each component (`button.module.css`); no Sass. Native nesting, custom properties,
  `color-mix()` and `@layer` cover what Sass did.
- **Tokens are `--msp-*` custom properties** in `styles/theme.css`; the reset is `styles/reset.css`; `styles/base.ts`
  imports both.
- **Themes are token sets** chosen by `data-msp-theme` on any element: `light` (default, on `:root`), `dark`, `blue`,
  `contrast`, `auto` (follows `prefers-color-scheme`); `prefers-contrast: more` gets the contrast edges. A theme sets a
  few base colors (background, text, accent, states, brand, `--msp-shade`: the direction shades go); every other color
  derives from them with `color-mix(in oklab, …)` wherever a theme is set, so a subtree with its own theme gets
  consistent shades. A host overrides any token on any element; tokens are declared on `:root` and `:host`, so they
  cross a shadow root.
- **Everything in `@layer msp.*`** (lowest first: `msp.reset`, `msp.theme`, `msp.components`, `msp.overrides`), so host
  styles outside the layers win. Only the reset uses `!important` (`[hidden]`, reduced motion). `@property`
  registrations and the shared keyframes sit outside the layers.
- **The old `msp-*` global classes are not reused**: prism-ui and the classic UI coexist on one page.

## 2. Tokens

All in `styles/theme.css`, by group: surfaces (`--msp-color-bg`, then derived `-panel`, `-raised`, `-control`, `-hover`,
`-pressed`, `-border`, `-border-strong`; `-well` a step the other way), text (`-text`, `-text-muted`, `-text-subtle`),
states (`-accent`, `-accent-strong`, `-danger`, `-success`, `-warning`, `-info`, `-focus-ring`, `-scrim`), brand
(`--msp-color-brand-1…4`, the logo's colors), glass (`--msp-glass-bg`, `-panel-bg`, `-filter`, `-edge`, `-shadow`), type
(`--msp-font-*`, `--msp-line-height`), size (`--msp-row-height`, `--msp-control-height`, `--msp-label-width`,
`--msp-radius-sm/md/lg`), elevation (`--msp-shadow-*`, `--msp-z-*`), motion (`--msp-motion-fast`, `-base`, `--msp-ease`,
`--msp-ease-reveal`), the over-canvas ink (`--msp-ink*`), fades (`--msp-fade*`, registered in `reset.css`).

- **Layout values are set inline by the shell**, not the theme: `--msp-dock-size`, `--msp-canvas-*` insets, the rail's
  `--msp-rail-*` steps. In compact mode fixed-size targets grow to `--msp-touch-target`.
- **A module may define its own tokens** for what varies inside one component (`--msp-tree-indent`, `--msp-progress`)
  and override a shared one locally (Settings' `--msp-label-width`).

## 2.1 Spacing

Every padding, margin and gap is a step of the scale or a rule named after it:

| Token                  | Value                  | For                                                                               |
| ---------------------- | ---------------------- | --------------------------------------------------------------------------------- |
| `--msp-space-1` … `-6` | 2, 4, 8, 12, 16, 24 px | list rows; inside a row; groups in a form; a card's padding; blocks; empty states |
| `--msp-inset`          | 8px                    | where text starts inside a control or a row                                       |
| `--msp-panel-inset`    | 4px                    | a panel's content inset at its sides and bottom                                   |
| `--msp-canvas-edge`    | 8px                    | the margin of everything floating over the canvas                                 |

- **Use the rule when the value means it** (a row's text inset is `--msp-inset`, not `--msp-space-3`).
- **`styles:check` rejects** a padding, margin or gap in px other than a 1px hairline; geometry that is not spacing
  stays in px with an `off-scale: <why>` comment on its line (`scripts/css.mjs`). Stories' inline styles use tokens too.

## 3. Build

- **Each `.css` under `src/` compiles to `_generated/x.css.ts`** next to it (committed, excluded from Biome) with
  lightningcss: a `*.module.css` gets classes and keyframes renamed `msp-<component>-<name>` (duplicate prefixes are an
  error), wrapped in `@layer msp.components`, and exports the typed `classes` map; the module registers its CSS with
  `defineStyles(id, css)`. Grid areas, container names and names starting with `msp-` (the shared keyframes) stay
  global.
- **A component imports its generated module**, so a page carries exactly the CSS of what it imports. Edit the `.css`,
  never the generated module; regenerate after formatting.
- **Styles are adopted, never linked**: `adoptStyles(document | shadowRoot)` keeps one constructed sheet per document (a
  `<style>` where unsupported), shared by its shadow roots and updated on later registrations (a lazy extension's; a
  burst applies once, on a microtask). A top-level element gets `.msp-root`, which scopes the reset (`:where`), so host
  pages are untouched. `styleText()` gives the whole CSS for a static file.
- **Tooling**: `pnpm styles`, `styles:check` (missing, stale or orphaned modules; the spacing scale), `styles:watch`;
  Storybook runs the same generator as a Vite plugin (`cssModulesPlugin`). Packages built on prism-ui run its
  `scripts/css.mjs` with their own `--prefix`.
- **Icons are generated**: `scripts/icons.mjs` copies the used Lucide icons (the `ICONS` list, `lucide-static` pinned)
  into `components/icon/_generated/lucide.ts`; `pnpm icons`, `icons:check`. Molecular glyphs are hand-drawn in
  `components/icon/molecular.ts` on the same 24px, 2px round-stroke grid.

## 4. Rules

- **Flat, few lines.** Areas part by background tone and spacing; a hairline only where a block would merge with a
  surface of the same tone. Everything over the canvas is flat glass (`--msp-glass-*`, the hairline as a spread shadow,
  no drop shadow); what opens over content for a moment (the command menu, menus, floating options) adds a soft shadow.
- **Panel text is 12px** (`--msp-font-size-sm`); 11px for captions and hints; 14px only where something reads a step
  larger (what is typed in the command menu, the canvas chips).
- **One text column, 8px into its box**: panel content inset 4px, every box pads its text 8px, so text starts 12px from
  the panel's edge, the header title's column.
- **Wells and inline blocks** set a block apart: a well (`--msp-color-well`) for an open section or options; an inline
  block (the control color at 40%, radius md, inset hairline) for what a row opens under itself. No accent bars.
- **Icons only where they say what text does not**: icon-only controls, state, top-level tree rows. Not on child rows,
  menu items, headers or beside words (except + Add). `Icon({ name, size, label })` is decorative unless labelled.
- **Accent for actions and the picked item only** (a primary button, the picked option, the active tab, the focus ring);
  settings (a slider's fill, a switch that is on) are `--msp-color-text-muted`.
- **A row's actions show on its hover, always on touch** (David, 2026-10-21): icon actions beside a row's text (an eye,
  a ⋯) are `opacity: 0` until the row is hovered or holds the focus, or the action's menu is open; under
  `@media (hover: none)` always shown, since a touch screen has no hover to find them by. What tells the row's state or
  where it leads stays shown: a hidden item's eye, a chevron (at the row's edge, the actions left of it). Why: quiet
  rows on a mouse, nothing out of reach on a finger.
- **The arrow cursor**; the text cursor only where text is typed; special cursors only for special gestures.
- **Cut-off text fades out, no ellipsis**: `msp-fade-overflow` (`reset.css`, CSS only, on the element's inline scroll
  timeline); `data-tip-overflow` shows the full text in the tooltip only while cut off.
- **Motion**: state changes take `--msp-motion-fast`; opening content grows and fades in `--msp-motion-base` (`Show`'s
  `reveal`, `dom/motion.ts`); a virtual tree's opened node too (its run of rows comes in, view-models.md §3: the user
  expects it to open as a tree did before it was virtualized). Nothing animates on first render, on a tab switch, for a
  list's other changes or its scrolling, or under reduced motion (the reset finishes every animation at once).
- **Every scrolling area is a `ScrollArea`** (`components/scroll-area`): overlay scrollbars, fading edges where there is
  more to scroll (`msp-fade-scroll-*`).
- **Shared controls**: `Select` opens inline (or floats with `over`); `ParamsForm` renders `PD.Params` as they are
  (essentials, then Advanced; on/off params as a switch); one `ColorPicker`; one `ActionFooter` ends every block applied
  by a button, and Enter in a one-line field of that block applies it (as Mol* v5's actions: the field commits first;
  not in a text area, a list box, a check or a slider, nor while it cannot apply); `Choices` for a row of one-click
  actions.
- **Focus rings** are a 2px `--msp-color-focus-ring` outline; fields, rows and footer buttons draw theirs inside.
- **The canvas follows the theme when it is plain**: each theme has `--msp-canvas-bg`; choosing a theme sets it as the
  3D canvas's background only while that background is a plain color (no image, gradient or skybox style), and at the
  start unless the options give one (`viewer/src/ui/shell.tsx`).
- **Colors, spacing and radii from tokens**; literals only for what must read whatever the theme. No inline styles
  except runtime values (positions, sizes, progress, data colors), often as custom properties.
- **State is styled through attributes the DOM already carries** (`[aria-expanded]`, `[aria-selected]`,
  `[data-dragging]`); `classList` toggles classes by reactive flags.
- **Fields have no edge** in the default themes (`--msp-field-edge` transparent); the contrast theme and
  `prefers-contrast: more` draw a 3:1 edge.

## 5. Performance

The UI must feel instant. Budgets: an input answered in the frame it arrives in (a pointer handler under 2 ms, a click's
work under 16 ms or deferred behind visible feedback); no task over 50 ms from an interaction; no layout shift when
content loads.

- **No layout thrashing.** Never read layout (`getBoundingClientRect`, `offset*`, `client*`, `scroll*`,
  `getComputedStyle`) after writing in the same frame; read from `ResizeObserver` entries or once at a frame's start. A
  positioned popup measures once when it opens.
- **Pointer-rate paths** (hover, drag, scroll, wheel): one update per frame, the DOM touched only on a change,
  `transform` rather than `top`/`left`; no nodes, list rebuilds, closures or objects per event. Delegate events from a
  container.
- **Bind the smallest node to the smallest observable**; never rebuild a subtree for one value. `For` keyed over stable
  ids, `VirtualList` for lists that can outgrow a screen, dense data (residues, tracks) on a canvas.
- **Style recalculation**: toggle state on the element that changes, not an ancestor; selectors one class deep, no
  universal or long descendant selectors, no `:has()` on large containers.
- **Compositing**: animate only `transform` and `opacity`; `will-change` only while animating. `backdrop-filter` over
  the 3D canvas re-blurs every frame: keep glass there small (a frameless tab draws none), likewise large blurred
  shadows.
- **Hover is instant**: no transition on a hover's colors (it reads as lag and makes neighbours blink on a sweep).
- **Containment**: panels and cards set `contain` (`layout paint`, or `content` where the size is known), so a change
  inside does not relayout the shell; off-screen rows use `content-visibility: auto` where a list is not virtualized.
- **Hidden UI is inert and does no UI updates.** What is hidden but kept mounted (a folded section, a panel behind an
  overlay card, a closing animation) is `inert` (no pointer, tooltip, focus or screen reader) and **paused**: its scope
  is inactive (`activity(flag, …)` in `dom/control-flow`, `Scope.setActive`), which drops every `bind` subscription
  under it. Showing again resubscribes once and each state replays its value: nothing rebuilt (the DOM, scroll position
  and drafts stay). Views do nothing for this; one that runs its own work (an observer, a rAF loop, a timer, a store
  subscription) pauses it with `onActive(active => …)`. Prefer removing content that is cheap to build (`Show`); keep it
  mounted and paused where state or cost matters. Tooltips never show for an inert element.
  `views/shell/_test/hidden.test.tsx`.
  - **Checked on the real shell** (`viewer/src/ui/_test/hidden-quiet.test.tsx`, jsdom, 1CBS): with the Structure tab
    shown and the toolbar folded, the Selection tab with its sections folded, and the compact rail closed over the
    docks, a load, picks, a selection, hide and show, renames and trajectory frames write nothing inside `[inert]` but
    the `inert` toggles (a MutationObserver on the shadow root), and the runtime's probe (`RuntimeProbe`) counts no
    binding write and no `watch` render in a paused scope, against 156, 44 and 39 in what shows. A raw subscription
    planted in a hidden part is caught. Tabs that are not shown are not mounted at all; the hidden-but-mounted parts are
    the folded toolbar, folded sections and docks behind the compact rail.
  - **New base primitives and hot views come with a guard** (packages/prism/AGENTS.md, Tests): a quick case (bytes a
    step, its positive control) and a heap guard; the base primitives (a state set, derive, bind, a `watch` render and
    flush, the keyed selector and store) read 0 bytes a step.
  - **A view subscribes to the narrowest source that determines its output** (view-models.md §3): hidden or shown, a
    change that writes nothing is budgeted waste (`waste.test.tsx`, `waste-budget.json`).
  - **Imperative views**: the sequence band's track gutter rebuilt its icons on every Seq\* layout version with raw
    subscriptions; paused, it now marks itself stale and rebuilds once on resume. Left as they are: Seq\*'s section
    surface keeps scheduling its frames while hidden (it skips drawing at no size; pausing it needs a Seq\* API); the
    PAE plot (model-archive) follows its colors by `bind` now (it pauses with its card); menus' and popovers' document
    listeners and observers live only while they are open; one-shot frames (a scroll into view, a font read) do not
    repeat.
  - **What the UI costs while a structure animates** (Internals › Performance › Animation, `watch-bench/anim.mjs`:
    1GRM's 5 NMR models looping at 5 a second, 6 s a setup, Chromium on an M4 with its GPU, a CPU profile split by
    source): every setup holds 60 fps (16.7 ms frames, p95 18.7 ms); the main thread works 0.37–0.53 ms a frame, the UI
    (Prism, Seq\*) 0.008–0.024 ms of it:

    | setup                 | busy ms/frame | UI ms/frame (share) | Mol\* ms/frame | writes | renders | DOM mutations/frame              |
    | --------------------- | ------------- | ------------------- | -------------- | ------ | ------- | -------------------------------- |
    | all panels closed     | 0.51          | 0.010 (2%)          | 0.37           | 0.3    | 0       | trajectory bar 0.2               |
    | the Structure tab     | 0.45          | 0.010 (2–3%)        | 0.32           | 0.3    | 0       | trajectory bar 0.2               |
    | the state tree        | 0.43          | 0.015–0.024 (4–6%)  | 0.28           | 1      | 0.8     | state tree 0.5, bar 0.2          |
    | the sequence band     | 0.38–0.42     | 0.008 (2%)          | 0.28           | 0.2    | 0       | bar 0.2 (Seq\* draws its canvas) |
    | the command menu open | 0.37–0.40     | 0.009–0.013 (2–3%)  | 0.27           | 0.3    | 0       | bar 0.2                          |

    The rest of a frame's work is the engine's (deps, native, GC: 0.08–0.15 ms). The state tree is the one view that
    works per model: a model's cell record is replaced, so its row is made anew (`For` by item). The scene tree and the
    inspector write nothing while the models play. Long animation frames: at most one in 6 s per setup.
- **Fonts and images** reserve their size; fonts have metric-compatible fallbacks.
- **Check it**: record the interaction in the browser's performance panel (CPU throttled 4x): forced reflows, long style
  recalculations, long tasks. A story with a heavy view is the place to measure.
