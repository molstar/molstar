# Styling

## 1. Decisions

- **CSS Modules** for component styles (`button.module.css` next to `button.tsx`): local class names, dead-class
  detection, no global naming scheme to police. Replaces SCSS; native CSS nesting, custom properties, `color-mix()`, and
  `@layer` cover what Sass mixins and colour functions do today.
- **Design tokens as custom properties with the `--msp-` prefix**, defined once in `styles/theme.css`; the reset is
  `styles/reset.css`. `styles/base.ts` imports both (every prism-ui root needs them).
- **Themes are token sets** switched at runtime by `data-msp-theme` on any element (`light`, the default on `:root`;
  `dark`; `blue`; `auto` following `prefers-color-scheme`). A theme sets a few base colors (background, text, accent,
  the states, the brand colors, and `--msp-shade`, the direction shades go: black on light themes, white on dark ones);
  every other color is derived from them with `color-mix(in oklab, …)` wherever a theme is set, so a subtree with its
  own theme gets consistent shades. A host can override any token on any element. The shell switches the UI theme from
  the rail's foot, and themes the canvas chrome separately from the canvas color ([shell.md](shell.md)).
- **Everything in `@layer msp.*`** (reset, theme, components), so host styles outside the layers win when they target
  prism-ui deliberately. Only the reset uses `!important` (`[hidden]`, reduced motion).
- **Light DOM** ([rendering.md](rendering.md) §5): one stylesheet styles all elements; no per-shadow-root injection.

## 2. Tokens

All in `styles/theme.css` unless noted. "Derived" means computed from the theme's base colors (the percentage is the
step toward `--msp-shade`, or toward the background for text).

| Group       | Tokens                                                                                                                                                                                                                                                                                                                              |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Surfaces    | `--msp-color-bg` (the theme's base); derived, lowest to strongest: `-panel` (3%), `-raised` (5%), `-control` (6.5%), `-hover` (9%), `-pressed` (13%) (the dark themes step further: 10%, 15%, 20%, `--msp-step-*`), `-border` (15%), `-border-strong` (28%); `-well` is the base color itself (a step from the panel the other way) |
| Text        | `--msp-color-text` (the theme's); derived `-text-muted` (labels, secondary text) and `-text-subtle` (hints, placeholders, disabled)                                                                                                                                                                                                 |
| States      | `--msp-color-accent`, `-accent-strong` (accent text on an accent wash), `-danger`, `-success`, `-warning`, `-info` (the theme's); derived `-accent-muted`, `-focus-ring` (the accent); `-scrim`. `-info`, `-accent-muted` and `-scrim` are defined, not used by the mocks yet                                                       |
| Brand       | `--msp-color-brand-1…4`: the logo's colors left to right (the wordmark and its flash), lighter on the dark and blue themes                                                                                                                                                                                                          |
| Glass       | `--msp-glass-bg` (the panel color, 70% opaque), `--msp-glass-panel-bg` (95%, for larger surfaces read over the canvas), `--msp-glass-filter` (blur and saturate), `--msp-glass-edge` (a hairline color), `--msp-glass-shadow` (that hairline as a 1px spread shadow: no layout)                                                     |
| Node kinds  | `--msp-accent-cyan`, `-red`, `-gray`, `-green`, `-purple`, `-blue`, `-orange` (the classic card accents for extension and node-kind brands; defined, not used by the mocks yet)                                                                                                                                                     |
| Type        | `--msp-font-sans`, `--msp-font-mono`, `--msp-font-size-xs` (11px), `-sm` (12px), `-md` (14px), `-lg` (16px, unused), `--msp-line-height` (1.4)                                                                                                                                                                                      |
| Size        | `--msp-row-height` (32px), `--msp-control-height` (28px), `--msp-label-width` (110px; Settings overrides it), `--msp-space-1…6` (2, 4, 8, 12, 16, 24px), `--msp-radius-sm` (2px, small inline marks), `-md` (4px, controls and rows in a surface), `-lg` (8px, every surface floating over the canvas)                              |
| Elevation   | `--msp-shadow-overlay` (floating select options), `--msp-z-overlay`, `-dialog`, `-tooltip`; `--msp-shadow-float` and `--msp-z-toast` are defined, not used yet (the other drop shadows are literals in their modules)                                                                                                               |
| Motion      | `--msp-motion-fast` (100ms), `--msp-motion-base` (220ms), `--msp-ease`, `--msp-ease-reveal` (ease in and out, for content that opens and closes); under `prefers-reduced-motion` the reset finishes every animation and transition at once                                                                                          |
| Layout      | Set inline by the shell, not the theme: `--msp-dock-size` per dock, `--msp-canvas-edge` (8px), `--msp-canvas-top-inset`, `-right-inset`, `-rail-inset`, `-panel-right`, `-pinned`, `--msp-right-panel-width` ([shell.md](shell.md))                                                                                                 |
| Over-canvas | `--msp-ink`, `--msp-ink-muted`, `--msp-ink-shadow` under `[data-msp-ink="light"]` and `[data-msp-ink="dark"]` (text over the canvas following its background). Defined, not set anywhere yet: the canvas chrome follows the canvas through the canvas layer's own theme and base color ([shell.md](shell.md))                       |
| Fades       | `--msp-fade`, `--msp-fade-top`, `--msp-fade-bottom`: registered (`@property`) lengths in `styles/reset.css`, animated by the fades below                                                                                                                                                                                            |

Modules define a few tokens of their own for what varies within one component (`--msp-tree-indent`, `--msp-range-fill`,
`--msp-progress`, `--msp-stack-row`), and may override a shared one locally (Settings' `--msp-label-width`).

The light, dark and blue backgrounds, text and accent colors are those of today's Sass skins (the old UI's
`skin/{light,dark,blue}.scss` over `skin/base/_colors.scss`; light is the inverse of dark), so the default look stays
recognisable; the shades derive from them rather than copying each Sass variable. The state colors start from the skins'
log colors, adjusted where they would not read on a theme's background (darker on light; `-info` lighter on dark and
blue; blue's `-danger`). Not yet: state-tree type-class colors and sequence colors as tokens.

## 3. Build

Implemented in `scripts/css.mjs` and `src/styles/registry.ts`. `tsc` leaves CSS imports untouched in `lib/`, which
neither Node nor bundler-free browsers can load, so CSS is compiled into TypeScript:

1. **Each `.css` file under `src/` gets a generated module `_generated/x.css.ts` in its folder** (committed, excluded
   from Biome). lightningcss compiles it; for `*.module.css` it renames classes and keyframes deterministically to
   `msp-<component>-<name>` (the file name, with its folder when they differ: `select`, `structure-scene-tree`;
   duplicate prefixes are an error) and wraps the result in `@layer msp.components`. Grid area and container names stay
   global (layout contracts between modules: the shell grid, the docks, the canvas container). Names that already start
   with `msp-` are global too: a module refers to the shared keyframes in `styles/reset.css` (`msp-fade-overflow`,
   `msp-fade-scroll-*`) as they are, so a module's own names never start with `msp-` (they are left out of `classes`).
   The module exports the typed class map (`classes`) and registers the CSS with `defineStyles(id, css)`.
2. **Components import the generated module** (`import { classes } from './_generated/button.module.css.js'`), so a page
   carries exactly the CSS of the components it imports.
3. **Elements adopt the styles.** On connect, every `MspElement` calls `adoptStyles(document)`: one constructed
   stylesheet per document (a `<style>` element where constructed sheets are unavailable), updated when more modules
   register (lazily loaded views; a burst of registrations applies once, on a microtask). A top-level element (not
   inside another prism-ui root) gets `.msp-root`, which scopes the reset (`:where(.msp-root)`, so host pages are
   untouched). The tokens are on `:root` and on every `[data-msp-theme]`, not on `.msp-root`. No stylesheet to link, in
   any build setup; `styleText()` gives the whole CSS for a static file if a distribution wants one.
4. **Tooling**: `pnpm styles` regenerates, `pnpm styles:check` fails on missing, stale or orphaned modules (for CI),
   `pnpm styles:watch` regenerates on change, and Storybook runs the same generator as a Vite plugin
   (`cssModulesPlugin`). Jest maps `./_generated/x.css.js` to `x.css.ts`. Edit the `.css` file, never the generated
   module; Biome formats `.css`, so regenerate after formatting.
5. **Icons** are generated too: `scripts/icons.mjs` copies the inner markup of the Lucide icons the UI uses (the `ICONS`
   list, from `lucide-static`, a pinned devDependency) into `components/icon/_generated/lucide.ts` (committed).
   `pnpm icons` regenerates, `pnpm icons:check` fails when it is missing or stale.

Layers, lowest first: `msp.reset`, `msp.theme`, `msp.components`, `msp.overrides` (declared by `styles/registry.ts`;
`msp.overrides` is empty so far). The `@property` registrations and the shared keyframes in `styles/reset.css` sit
outside the layers.

## 4. Rules

- **Minimal borders, flat surfaces.** Areas are separated by background tones (`--msp-color-bg`, `-panel`, `-raised`,
  `-control`, `-well`) and spacing, not lines. A hairline only where a block would otherwise merge with a surface of the
  same tone: inline blocks, the inspector's line, glass over a canvas of the same tint. Everything over the canvas is
  flat glass: `--msp-glass-bg` (or `--msp-glass-panel-bg` for a floating panel), `--msp-glass-filter` and the
  `--msp-glass-shadow` hairline, no drop shadow. What opens over content for a moment (the command menu, the rail's
  menus, floating select options, the busy card) adds a soft drop shadow to lift it.
- **Panel text is 12px.** The root's size is 14px (`--msp-font-size-md`); panels, canvas cards and the command menu set
  `--msp-font-size-sm` on their content (section lists, Home, Settings, the Structure tab, Measure and Superpose,
  Snapshots, the playback card), so one size runs through every panel; `Button`'s default size is 14px, and every view
  passes `size="sm"` (12px) or `"xs"` (11px). 11px (`-xs`) for captions, hints, descriptions and secondary lines. 14px
  only where something should read a step larger: what is typed in the command menu, the canvas's selection and focus
  chips, the file drop's title.
- **One text column, 8px into its box.** Panel content is inset 4px (`--msp-space-2`) at the sides and bottom; every box
  inside (section headers, Home's action rows and captions, list rows, fields, dropdown buttons) pads its text 8px
  (`--msp-space-3`). Text therefore starts 12px from the panel's edge, the column of the panel header's title (padded
  12px), and card headers on the canvas keep the same 12px. Controls keep the column too: a dropdown's and a text
  field's text, a slider's track and a switch all start 8px into the control column; nested params step their labels in
  by 8px per level while their controls stay in the form's column. Heights: section headers `--msp-row-height`, buttons
  and dropdowns `--msp-control-height`, form controls 24px (rows 2px apart, the label column `--msp-label-width`).
- **Wells and inline blocks.** Two ways to set a block apart inside a panel:
  - a **well** (`--msp-color-well`): an open section, an open Home action, the Structure tab's options, presets and
    styles, Snapshots' blocks, the image preview stage. A step from the panel the other way, so the controls in it stand
    out; a list inside a well (Measure's points, Superpose's entries) steps back to the panel color, so it reads as one
    block;
  - an **inline block**: what a row opens right under itself (a dropdown's options, the color picker, an on/off
    setting's settings, a list param's items, a param's help, a measurement's style, one's own endpoint in the Assistant
    tab), and an agent turn in History (dimmed once undone). Its surface is the control color mixed 40% toward
    transparent (a measurement's style sits on the well color instead), radius md, with an inset hairline
    (`box-shadow: inset 0 0 0 1px var(--msp-glass-edge)`) that keeps it apart inside a well of the same tone. No accent
    bar (the mocks dropped it).

  On dark themes the well is darker than the panel (shades go toward white), so a small form in a well of its own reads
  as a dark box: a small form sits on the panel like any control (the Assistant tab's key forms), and what needs a frame
  is an inline block.

- **Icons where they say something text does not.** Icons are the label of icon-only controls (the rail, the canvas
  stack, row buttons: eye, ⋯, settings, drag handle, remove), show state (visibility, a check, progress, a message's
  kind), anchor top-level objects in a tree (a structure, a volume, the measurements), and help scan the command menu.
  They are not repeated next to words: no icons on child tree rows (indentation and the description say what a row is),
  menu items, section and group headers, footer actions (a spinner only while running), or labelled header buttons
  (their icons only in the narrow layout, where the words go); the one kept beside a word is the + before Add
  (Snapshots' header, an object list's last row). Set defaults (`list-restart`) and refresh (`refresh-cw`) never share
  an icon. Sliders (`sliders-horizontal`) mean a row's own params unfolding under it (an on/off setting's, a list
  item's, a measurement's style) and nothing else; a panel's options are a gear (`settings`). The set is Lucide (the
  used icons only, generated, §3) plus `MolecularIcons` (`components/icon/molecular`: representations, picking levels,
  structure, volume, superpose, the Mol\* mark) drawn on the same 24px grid in the same 2px round stroke.
  `Icon({ name, size, label })` returns an `<svg>` in `currentColor`, 18px by default (rows mostly use 14px), decorative
  (`aria-hidden`) unless it has a `label`. A new icon is a name added to `ICONS` in `scripts/icons.mjs` and
  `pnpm icons`.
- **Accent for actions and the picked item only**: a primary button (a soft accent wash with `--msp-color-accent-strong`
  text), the picked option, the active tab, the focus ring. Settings are not actions: a slider's filled part and knob, a
  switch that is on, and the ring of a picked swatch are all `--msp-color-text-muted`. Inline blocks are set apart by
  their surface and hairline, not an accent bar.
- **The arrow cursor; the text cursor only where text is typed.** An app, not a page: the reset sets `cursor: default`
  on the root and `text` on text inputs, textareas and editable content; buttons get the pointer (the arrow again when
  disabled). Rows that select (tree rows, command menu rows, snapshots) keep the arrow. A special cursor only for a
  special gesture: `grab` on drag handles (and a line graph's points), `crosshair` on the color field and a line graph's
  plot, `col-resize` on splitters, `help` on a label with help (also underlined, dotted).
- **Cut-off text fades out, no ellipsis.** A label that does not fit fades over its last 20px; the fade shows only when
  it overflows (a scroll-driven animation on the element's own inline scroll timeline, `msp-fade-overflow` in
  `styles/reset.css`, CSS only; where scroll timelines are unsupported the text is clipped). An element opts in with
  `overflow: hidden`, `white-space: nowrap`, the `mask-image` gradient over `--msp-fade` and the animation as longhands
  (the recipe is in `reset.css`: the `animation` shorthand cannot carry a timeline). `data-tip-overflow` shows the full
  text in the app's tooltip, only while it is cut off. Nothing uses `text-overflow: ellipsis`.
- **Motion** (Overview › Foundations › Motion): state changes (hover, press, a chevron turning) take
  `--msp-motion-fast`; content that opens grows to its height and fades in, and closing plays it back, in
  `--msp-motion-base` (`Show`'s `reveal`, `dom/motion.ts` for content hidden by hand); panels slide in from the side;
  option lists drop in. Nothing animates on the first render, on a tab switch, for the rows of a list, or with reduced
  motion.
- **Every scrolling area is a `ScrollArea`** (`components/scroll-area`, adapted from Base UI's): overlay scrollbars that
  look the same on every platform (a thin thumb in the text color, wider under the pointer, shown on hover, while
  scrolling and while dragged, none when nothing overflows; a press on the track jumps there and drags on) and edges
  that fade where there is more to scroll (`msp-fade-scroll-*` in `styles/reset.css`, CSS only: no top fade until
  scrolled, no bottom fade at the end; areas that scroll only across do not fade). `class` goes to the root, which
  carries the area's size and surface (background, radius, hairlines), so those never fade; `contentClass` goes to the
  content, which carries the layout; `axis` is `y` (default), `x` or `both`; `viewportRef` gives the scrolling element.
  Nothing else sets `overflow: auto`.
- **Selects open inline, or float when asked.** `Select` expands like Mol\*'s action menu: its options (an inline block,
  the current one checked and in the accent) push the content below down, right under the button, or under the whole row
  in forms (`SelectButton` and `OptionList`). With `over` the options float instead: a manual popover in the top layer,
  under the button and at least as wide (above it when there is no room below), following it as things scroll, as wide
  as its options need, opaque (`--msp-color-raised`, hairline, `--msp-shadow-overlay`), so no panel's clipping or fade
  cuts it. `over` is for short settings rows where moving the layout would be jarring (the Assistant tab).
- **The params form** (`ParamsForm`, `components/params`) renders `PD.Params` as they are, one control per type,
  essentials first and the rest under Advanced. Group and Advanced headers are quiet rows (regular weight, muted, a
  small chevron after the name), not titles; a group inside a group steps its header in like its labels; a flat group
  (PD's `isFlat`) has no header, its rows at the parent's level. What a row opens unfolds under it as an inline block in
  the form's control column (a select's options, the color picker, a multi-select's checks, a mat4's grid, a line
  graph's editor, a color list's presets). Mol\*'s on/off pattern (a mapped param whose options are exactly `on` and
  `off`: occlusion, shadows, outline) is a switch and, while on, its settings behind the sliders button. An object list
  is its count in the row and one quiet list of its items under it: a click anywhere on an item opens its form inside
  the list, remove on hover, + Add as the last row. A param's description and legend (a table or a scale; `help(value)`
  when it has one, so a theme choice shows the chosen theme's) are behind a `?` after its label (on hover, or while
  open), opening under the row; rows carry no native title. Multiline text is a growing field; files a field like a
  dropdown's that also takes drops; plain values (`PD.Value`) have no row.
- **One color picker** (`ColorPicker`, `components/color`; no system picker), opened under the row: a smooth palette
  (hues across, lighter up, darker down, greys under it), the palette's hues, then the picked color's eight weights
  (itself among them, placed by its lightness); a field that reads most ways colors are written (`parseColor`: hex,
  0–255 or normalized channels, `rgb()`, `hsl()`, CSS names) and shows decimal channels. An empty color field says what
  it takes.
- **One action footer** (`ActionFooter`, `components/action-footer`) ends every block applied by a button (Home's
  actions, the state tree's updates, the Selection tab's sections, Animation, the remote upload, every Export section):
  set defaults on the left (only with a reset), optionally a quiet info line (the image's size) and a lesser action in
  the text color (Copy, Clear), then the action at the right in the accent, only as wide as its label (the empty row
  between is not a button, so a stray click runs nothing). No lines or fills; hover darkens and underlines. `disabled`
  (not possible yet) is apart from `busy` (running: a spinner).
- **Choices** (`Choices`, `components/choices`): a labelled row of one-click actions (the quick styles' Representation
  and Style rows) at equal widths, wrapping into rows of equal widths; `active` marks the one in effect.
- **Focus rings** are a 2px outline in `--msp-color-focus-ring`, 1px outside by default (the reset); fields, dropdowns,
  footer buttons and rows draw theirs inside (`outline-offset: -2px`), so a control at a scrolling panel's edge never
  has its ring cut off.

- Colors come from tokens; literal colors only for what must read whatever the theme (the color field, a marker's ring,
  white text on a danger fill, drop shadows). Spacing and radii come from tokens; a component's own sizes (a 24px form
  control, a 26px tree row) are literals in its module.
- No inline styles except values computed at runtime: positions, sizes, progress and data colors (swatches, gradients),
  often as custom properties (`--msp-dock-size`, the canvas insets, `--msp-progress`). The few static ones left: the
  wordmark's gradient stops, a view host's `display: contents`, the mock canvas placeholder.
- State is styled through the attributes the DOM already carries (`[aria-expanded]`, `[aria-selected]`,
  `[aria-checked]`, `[data-dragging]`, `[data-axis]`); classes compose with template strings, and the `classList` prop
  toggles classes by (reactive) flags.
- The old `msp-*` global classes are not reused; prism-ui and the old UI can coexist on one page.
