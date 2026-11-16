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

## 3. Spacing

- One unit (`--msui-unit`, 4px); the scale is its multiples (4, 8, 12, 16, 24, 32). Every padding, margin, gap and size
  is on the grid; only hairlines (1px) are off it. Why: anything placed by the scale lands where something else placed
  by it does.

## 4. Sizes

- A row is 24px; a form row is a list row; a control is a unit shorter (20px) and centered in its row, so stacked
  controls keep a gap and their text sits where a row's text does; a header is 32px; an icon 16px. Why: a list beside a
  form lines up row by row.

## 5. Alignment

Things that belong together visually share a column; each column is a token, so a change moves everything on it.

- **One inset** (`--msui-inset`, 12px): a header's, a row's and a form label's text start that far from their surface's
  edge. A surface pads its sides by `--msui-surface-pad` (4px, so a hovered row keeps off its edge); what is in it pads
  the rest (`--msui-inner`).
- **Icon slots** (`--msui-icon-slot`): a row's icon has a fixed slot, kept when empty (`.msui-slot` without content), so
  the text after it starts in one column with or without an icon. A list without icons uses `.msui-row-plain`.
- **A level is a slot** (`--msui-indent` = slot + gap, 16 + 4 = 20px): a child's icon sits under its parent's text
  (`--msui-level`). Why so small: a tree goes several levels deep in a narrow panel.
- **One label column** (`--msui-label`, 112px) for every form (`.msui-form`): values line up across forms and sections.
- **Values start where controls write** (`--msui-control-inset`): a plain value (`.msui-value`) starts where a control's
  text does, so read-only and editable rows line up.
- **A control sits in a row** (§4): a form beside a list lines up row by row.
- **Numbers by place** (`.msui-num`): tabular, right-aligned; units in their own column.
- **On the grid** (§3, §2): sizes, gaps and line heights.
- The Alignment stories draw the columns over a panel built with the primitives, and over the same content laid out ad
  hoc.

## 6. Shape

- Radii by role: `sm` (2px) inline marks, `md` (4px) controls and rows inside a surface, `lg` (8px) surfaces. Nested,
  the inner radius is smaller, a padding apart, so inner corners follow the outer curve.
- Elevation (`--msui-shadow-float`, `--msui-shadow-overlay`) only for what floats; `--msui-edge` for a surface on its
  own color. Layers: `--msui-z-*`.

## 7. Motion

- `--msui-motion-fast` (100ms: hover, press, a chevron), `--msui-motion-base` (200ms: what opens or closes),
  `--msui-motion-slow` (360ms: what only fades); `--msui-ease` (out) and `--msui-ease-in-out` (opening and closing).
- Move what the eye follows; pick the move by what must not change (grow in a flow, slide a panel, fade what keeps its
  size); transform, opacity and grid tracks only; nothing on a first render or for rows of a list. Every duration is 0
  under `prefers-reduced-motion`.
