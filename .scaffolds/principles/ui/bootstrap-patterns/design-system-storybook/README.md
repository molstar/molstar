# Design system and Storybook {#design-system-storybook}

- **Implementations:** [TypeScript](typescript.md)

Rules §design-system, §thin-view-layer: one token-driven design system of presentational components; Storybook is the workbench and the documentation hub; stories run on the real architecture and double as tests.

## Problem

Without one token-driven component library, each product invents its own colors, spacings and buttons, and themes
need per-component code. Components that reach into models cannot be shown or tested on their own. Without one
workbench, UI is judged only inside the running app, docs are copied and drift, and design mocks leak into product code.

## Shape

- **Tokens**, defined once in the common layer: groups `color`, `space`, `radius`, `size`, `font`, `motion`. They are
  emitted once as CSS custom properties (`--<group>-<key>`, numbers in px) on the root, with theme and density
  overrides as further token sets. Breakpoints are exported constants beside the tokens (a container query cannot read
  a custom property); layouts follow their container, and the shell is the container; a story of a component used
  outside a shell wraps it in one.
- **Styling**: scoped styles that use only tokens; variants as data attributes (`data-variant`, `data-state`); dynamic
  values through CSS variables set inline.
- **Presentational components** (the design system): props in, events out, no models; may keep transient local UI
  state (hover, focus, open/closed, a draft).
- **Connected views**, in the product's UI package at feature boundaries: bind a model (§reactive-models) to
  presentational components, reading only the slice they use and passing plain values and callbacks down.
- **Icons** (with Lucide as the example set; another open set works the same way): data in one registry in the
  common layer, in Lucide's icon-node form (a list of `[tag, attributes]` on a
  24×24 grid), Lucide's icons plus custom ones. One frame applies to all (`viewBox 0 0 24 24`, no fill, stroke
  `currentColor`, width 2, round caps and joins), and one `Icon` component renders them, sized by tokens, with an
  accessible label when it has one. An app bundles only the icon names its UI draws.
- **One Storybook per repository**, sections in order: docs, then `Reference/<design>` for each reference package,
  then one per package that declares a title in its manifest (`"storybook": { "title": "Design system/Components" }`),
  sorted by title. Story titles come from file names; folders such as `src/` and `stories/` never show.
- **Stories**: a presentational component's story is its props; a connected view's story is a model on the in-memory
  composition (`local()`, §ports-and-adapters), seeded through commands. One story per vocabulary sample covers every
  node kind (§conformance).
- **Checks at the gate**: the Storybook build; interaction and accessibility checks on every real story; visual
  regression baselines per story, reviewed on change. Reference stories are excluded.
- **Documentation hub**: pages render module specs, the command catalog (§commands), conformance status, error-code
  catalogs (§error-codes) and the rules from their sources (§rendered-docs); hand-written pages only for narrative
  with no other home.
- **Reference section** (§reference-mocks): throwaway mocks in `prototypes/<design>-reference/`, tagged `reference`,
  each naming its owning `.design/` file.

## Design flow (ties to §design-plan-implement)

1. `.design/<topic>.md` states intent and links the stories that illustrate it by story id (designs may reference Storybook, never `.execution/`).
2. Explore first in the reference package (fast, throwaway), then as real stories on the `local()` composition; review on the Storybook preview (a per-PR deploy helps).
3. Wire into the app only after the real stories settle, then delete the reference entries.

## Invariants

- Shared components use only tokens: no raw colors, spacing or font values (a lint), and no runtime CSS-in-JS.
- Themes and density are token sets, not per-component code.
- Shared components are restyled from outside only through documented CSS variables.
- Presentational components never touch models; only connected views do.
- A new package appears in Storybook by declaring its title, with no change to the Storybook config.
- No story needs a backend.
- Every shared component has stories for its states (default, empty, loading, error, long content, dense, disabled,
  dark/light, narrow) (§storybook-workbench).
- Every reference entry names an existing `.design/` file; nothing in `lib`, `packages`, `apps` or `services` imports
  a reference package; reference stories are never tested or baselined.
- Custom icons use only `path`, `circle`, `rect`, `line`, `polyline`, `polygon`, `ellipse`; carry no `fill`, `stroke`,
  `stroke-width`, `transform` or `style`; keep coordinates within 0–24 (a gate check).
- Product code names icons, never pastes SVG; a UI source that names an icon missing from the app's bundled set fails a quick test.
- Docs are generated from the registries and sources, never copied.

## Conformance cases

- Tokens `{space: {1: 4}, color: {bg: "#0b0d10"}}` → `--space-1:4px;--color-bg:#0b0d10;`.
- A package manifest with title `Design system/Components` → a section with that prefix; a package without one → no section.
- `prototypes/checkout-redesign-reference` → section prefix `Reference/checkout-redesign`.
- A reference story naming a `.design/` file that does not exist → the reference gate fails.
- A file under `lib/` importing from a `-reference` package → the reference gate fails.
- A custom icon with a `fill` attribute, a `text` element, or a coordinate of 30 → the icon check fails.
- Button story `Press`: click the button → `onPress` called once.
- The Storybook build succeeds, and every real story passes its interaction and accessibility checks.

## Open

- The icon check says coordinates stay within 0–24 but not how coordinates inside a path's `d` data are checked.
- §storybook-workbench lists the states every shared component's stories cover; no check enforces it.
