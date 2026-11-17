# Mol\* UI foundations: the specification

`@molstar/ui/styles/foundations.css`: the design tokens (`--msui-*`) and the alignment primitives (`.msui-*`). This spec
is the rules and their reasons; the stylesheet's comments and the foundation stories (`src/_stories/foundations/`,
Storybook › UI › Foundations) cite its sections (`styles §5`).

## 1. Color

- A theme (`data-msui-theme`: `light` (also the default on `:root`), `dark`, `blue`, `contrast`, `auto`) sets the base
  colors only: `--msui-bg`, `--msui-fg`, `--msui-accent`, the states, `--msui-shade` (the direction shades go: black on
  a light base, white on a dark one), the step sizes and the text mixes.
- Everything else is derived where a theme applies (`color-mix` in oklab): surfaces (`sunken`, the surface, `raised`),
  controls (`control` is raised; `hover` and `active` a step each above), text (`text`, `muted`, `subtle`: mixed toward
  the base as far as they keep 4.5:1 on a surface), lines, the muted accent, focus. Why: a subtree with its own theme
  gets consistent shades, and a host restyles by setting a few tokens.
- A dark base takes larger steps (`--msui-step-*`): a mix toward white reads as a smaller change than the same mix
  toward black.
- The accent is for actions and the picked item only; separation is by surface first, hairlines second.
- Contrast: the `contrast` theme and `prefers-contrast: more` keep muted text and lines at their contrast, and give a
  surface an edge (`--msui-surface-edge`): its step from a black base does not show.

## 2. Type

- Four sizes, each with a line height on the 4px grid: 11/16, 12/16 (the UI), 14/20, 16/24. Why: text of different sizes
  in neighboring columns lands on the same grid.
- Two weights: regular, and strong (600) for headers and the one action that matters.
- Tabular figures wherever numbers stack (`.msui-num`).
- **Overflow**: text that is cut off fades at its end (`msui-fade-overflow`, no ellipsis), but only variable content (a
  name, data in a row, a file name). The UI's own labels (buttons, choices, options, menu items, tabs) are never cut:
  they fit, their control grows, or their group wraps. Why: a cut label hides what an action does; a cut name is still
  recognizable and has its full text in a tooltip.

## 3. Spacing

- One unit (`--msui-unit`, 4px); the scale is its multiples (4, 8, 12, 16, 24, 32). Every padding, margin, gap and size
  is on the grid; only hairlines (1px) are off it. Why: anything placed by the scale lands where something else placed
  by it does.

## 4. Sizes

- A row is 24px; a form row is a list row; a control is a unit shorter (20px) and centered in its row, so stacked
  controls keep a gap and their text sits where a row's text does; a header is 32px; an icon 16px. Why: a list beside a
  form lines up row by row.

### 4.1 Icon sizing

- Three sizes on the grid: `--msui-icon-sm` (12px), `--msui-icon` (16px), `--msui-icon-lg` (20px). An icon is never
  taller than its text's line height.
- **By what the icon does:**
  - 16px where an icon identifies an object: in a list or tree row (24px), a header, a default (24px) button.
  - 12px where an icon marks a control or a command: in a 20px control (a small button, a select's chevron, a field's
    clear), in a menu or an option list (an item's icon, a check, a submenu chevron), next to small text.
  - 20px with large text and in empty states. Why: 12px text beside a 16px command icon reads heavy; an object's icon is
    what a row is recognized by, so it keeps the larger size.
- **An icon takes its size from where it sits.** The container sets `--msui-icon-size` (a button by its size, a slot, a
  header, a menu); an icon reads it (else `--msui-icon`). `Icon`'s `size` prop overrides it, for the rare icon that
  stands alone. Why: a button and its icon cannot drift apart, and one change to a container's size resizes its icons.
- An icon's slot keeps its width (`--msui-icon-slot`, 16) whatever the icon's size, so text stays in its column; a 12px
  icon is centered in it. An icon-only button is square at its height, the icon centered.
- Stroke icons scale their stroke with their size (a 24 viewBox, stroke 2): don't resize an icon by transform.
- **Toolbars** (`components/toolbar`):
  - A button is square on the grid, a header's height (`--msui-header`, 32px), with a 16px icon (`--msui-icon`): in a
    toolbar the icon is the label, so it takes the object size. A captioned button keeps the 16px icon, with xs text
    under it (a unit taller: 40px) or beside it (the icon where a square button has it).
  - A block of buttons is a surface: `--msui-radius-lg`, `--msui-surface-pad` between its edge and the buttons, whose
    radius is concentric (`--msui-radius-md`, §6); the gaps between blocks (`--msui-space-2`) are the grouping. Or
    (`grouping="separators"`) the bar is one block and a hairline (`--msui-line-strong`, an icon long) separates its
    groups, a unit of space on each side.
  - A menu is a surface of 24px rows: a 12px icon centered in the 16px slot at the inset, the text at inset +
    `--msui-indent`, a 12px check or submenu chevron on the inner right edge, hints in `--msui-text-subtle`, a group's
    heading xs muted on the text column. Labels are never cut (§2). The chosen one of a choice takes the accent text,
    not bold. A menu below or above its button is shifted by its padding, so its icon slot lies under the button's icon.
    Why: the bar's icons and the menu's share a column.

## 5. Alignment

Things that belong together visually share a column; each column is a token, so a change moves everything on it.

- **One inset** (`--msui-inset`, 12px): a header's, a row's and a form label's text start that far from their surface's
  edge. A surface pads its sides by `--msui-surface-pad` (4px, so a hovered row keeps off its edge); what is in it pads
  the rest (`--msui-inner`).
- **Icon slots** (`--msui-icon-slot`): a row's icon has a fixed slot, kept when empty (`.msui-slot` without content), so
  the text after it starts in one column with or without an icon. A list without icons uses `.msui-row-plain`.
- **A level is a slot** (`--msui-indent` = slot + gap, 16 + 4 = 20px): a child's icon sits under its parent's text
  (`--msui-level`). Why so small: a tree goes several levels deep in a narrow panel.
- **One label column** for every form (`.msui-form`): `--msui-label` (112px), or wider when a form's longest label needs
  it, up to about half the form's width; one width for the whole form, so values line up across its sections. Labels
  never wrap: a label past the cap stays on one line and fades at its end, with its full text in its tooltip. Why: a
  wrapped label makes a row two or three rows tall and breaks the form's rhythm; a field's label is the name of a
  parameter, recognizable from its start.
- **Nested form levels keep one value column.** A group, a choice's settings, an advanced fold or a list item inside a
  form is a level deeper: its labels indent by a form level (`--msui-form-indent`, 12px: form labels have no icon to
  align under, so a level is a step of the scale, not the tree's slot) inside the same label column, so every field at
  every depth starts on the same value column. Why: values that share a column can be compared at a glance; a form that
  steps its fields right at each level reads as a tree, not a form.
- **Values start where controls write** (`--msui-control-inset`): a plain value (`.msui-value`) starts where a control's
  text does, so read-only and editable rows line up.
- **A control sits in a row** (§4): a form beside a list lines up row by row.
- **Numbers by place** (`.msui-num`): tabular, right-aligned; units in their own column.
- **On the grid** (§3, §2): sizes, gaps and line heights.
- The Alignment stories draw the columns over a panel built with the primitives, and over the same content laid out ad
  hoc.

- **What loads is a placeholder in place** (`components/skeleton`): its shapes on the columns and lines of the content
  that replaces it (a row's slot and text column, a label and a control, text lines), shown after a moment (300ms), so
  nothing moves when the content arrives and a quick load shows nothing. Its shapes pulse together (opacity only); its
  words are text in a polite status region (said as it appears), and the area that loads is `aria-busy` meanwhile
  (`Loaded`), not the status region (busy would hold its words back).

## 6. Shape

- Radii by role: `sm` (2px) inline marks, `lg` (8px) surfaces, `md` controls and rows inside a surface. Elevation
  (`--msui-shadow-float`, `--msui-shadow-overlay`) only for what floats; `--msui-edge` for a surface on its own color.
  Layers: `--msui-z-*`.
- **Squircles** where the browser draws them (`@supports (corner-shape: squircle)`): `--msui-corner-shape` is `squircle`
  there, `round` elsewhere, and the reset sets it on every element of a root, so one token switches a subtree. A
  superellipse reads tighter than a circle at the same radius, so with squircles the radii grow by half (`sm` 3px, `lg`
  12px). Why: a squircle's curvature starts gradually, so a corner joins its sides without the visible kink of a
  circular arc, and nested corners stay calm.
- **Nested corners are concentric.** An element inset by `p` inside a parent with radius `R` takes radius `R − p` (never
  below 0): both curves share their centers, so the gap between them is even all around the corner. A larger inner
  radius bulges toward the parent's corner, a smaller one leaves a wedge of extra space there. `--msui-radius-md` is
  derived this way from the surface (`--msui-radius-lg` minus `--msui-surface-pad`), so a row or control inside a
  surface is concentric by construction, with round corners and with squircles. Compose further levels the same way:
  `max(0px, calc(<parent radius> - <inset>))`; where the inset is the radius or more, the inner element is square or
  takes `sm`.
- **A popover covers its column**: opened from a panel (a surface), it takes the panel's left edge and width, never
  wider; where there is no panel, it is as wide as its content, at most 340px. Why: a floating panel wider than the
  column it opens from overhangs the layout and reads as detached from it.
- A fill that touches its parent's edge (no inset) takes the parent's radius on the touching corners and none inside.

## 7. Motion

- `--msui-motion-fast` (100ms: hover, press, a chevron), `--msui-motion-base` (200ms: what opens or closes),
  `--msui-motion-slow` (360ms: what only fades); `--msui-ease` (out) and `--msui-ease-in-out` (opening and closing).
- Move what the eye follows; pick the move by what must not change (grow in a flow, slide a panel, fade what keeps its
  size); transform, opacity and grid tracks only; nothing on a first render or for rows of a list. Every duration is 0
  under `prefers-reduced-motion`.

### 7.1 Transitions

What enters and leaves the UI, one vocabulary (tokens `--msui-enter`, `--msui-exit`, `--msui-rise`, `--msui-ease-in`):

- **Enter slower than leave.** What enters takes `--msui-enter` (base, eased out: it arrives and settles); what leaves
  `--msui-exit` (fast, eased in: it gets out of the way). Why: arrival is what the eye follows; a departure that lingers
  delays what was asked for.
- **A float rises from its anchor** (`msui-float-in`: fade, and `--msui-rise` from the side it opened on; `Floating`
  marks the side it took, a flip included). Every floating layer (a select's list, a popover, a menu, a context menu)
  enters this way and leaves at once: a pick or Esc is done. Over its field (a select) it only fades, so the current
  option stays on the field's text. A list to pick from (options, a menu) enters fast (`--msui-motion-fast`): it is
  waited on; a panel (a popover) at `--msui-enter`. Each layer sets both offsets itself: they are custom properties,
  which inherit, and a layer inside another (a select in a popover, a submenu) must not rise by its parent's.
- **A fold grows and shrinks** (`Show`'s `reveal`): to its height on opening, back to nothing on closing, both over
  base, eased in and out. What slides or fades in (`reveal` `side`, `fade`) enters at `--msui-enter` and fades out at
  the exit. The closing plays through the runtime's `Motion.close`, installed with the styles (`adoptStyles`, so any
  root that takes them, `styles/transitions.ts`); without Web Animations or under reduced motion the content goes at
  once. The durations are registered `<time>` properties, so the closing reads them resolved whatever a theme wrote (a
  `calc()`, a `var()`).
- **A swap crossfades in place**: what replaces content in the same place (a placeholder by its content: `Loaded`)
  shares its cell; the old fades out (exit) as the new fades in. Nothing moves when the two have the same shape (§5,
  placeholders).
- **Not animated**: a first render, the rows of a list (scrolling, filtering), what follows the pointer, values that
  change in place (a number, a label). Why: motion there is noise, and costs on hot paths.
