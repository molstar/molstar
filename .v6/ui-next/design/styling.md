# Styling

## 1. Decisions

- **CSS Modules** for component styles (`button.module.css` next to `button.tsx`): local class names, dead-class
  detection, no global naming scheme to police. Replaces SCSS; native CSS nesting, custom properties, `color-mix()`, and
  `@layer` cover what Sass mixins and colour functions do today.
- **Design tokens as custom properties with the `--msp-` prefix**, defined once in `styles/tokens.css`.
- **Themes are token sets** switched at runtime by `data-msp-theme` on the shell root (or any subtree) (`light`, `dark`,
  `blue` to keep the current variants, `auto` following `prefers-color-scheme`). A host can override any token on the
  element.
- **Everything in `@layer msp`** (reset, tokens, components), so host styles outside the layer win when they target
  ui-next deliberately, and ui-next never needs `!important`.
- **Light DOM** ([rendering.md](rendering.md) §5): one stylesheet styles all elements; no per-shadow-root injection.

## 2. Tokens

| Group       | Tokens (examples)                                                                                                                                                                       |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Colour      | `--msp-color-bg`, `-panel`, `-raised`, `-hover`, `-pressed`, `-border`, `-border-strong`, `-text`, `-text-muted`, `-accent`, `-danger`, `-warning`, `-success`, `-focus-ring`, `-scrim` |
| Accents     | `--msp-accent-cyan`, `-red`, `-gray`, `-green`, `-purple`, `-blue`, `-orange` (today's brand accents for extension views)                                                               |
| Type        | `--msp-font-sans`, `--msp-font-mono`, `--msp-font-size-xs…lg`, `--msp-line-height`                                                                                                      |
| Size        | `--msp-row-height` (32px today), `--msp-space-1…6`, `--msp-radius-sm…lg`, `--msp-label-width` (110px today)                                                                             |
| Layout      | `--msp-dock-left-size`, `-right-size`, `-top-size`, `-bottom-size`, `--msp-activity-bar-size` (the shell writes sizes here)                                                             |
| Elevation   | `--msp-shadow-float`, `--msp-shadow-overlay`, `--msp-z-overlay`, `-dialog`, `-toast`, `-tooltip`                                                                                        |
| Motion      | `--msp-motion-fast`, `--msp-motion-base`, `--msp-ease`                                                                                                                                  |
| Over-canvas | `--msp-ink`, `--msp-ink-muted`, `--msp-ink-shadow`, set from `data-ink="light                                                                                                           | dark"`on viewport overlays (derived from the canvas background in`ViewportModel.ink`) |

The light, dark, and blue themes are derived from today's Sass variables (`skin/base/_colors.scss`, `_vars.scss`), so
the default look stays recognisable. State-tree type-class colours and sequence colours become tokens too.

## 3. Build

Implemented in `scripts/css.mjs` and `src/styles/registry.ts`. `tsc` leaves CSS imports untouched in `lib/`, which
neither Node nor bundler-free browsers can load, so CSS is compiled into TypeScript:

1. **Each `.css` file under `src/` gets a generated module `_generated/x.css.ts` in its folder** (committed, excluded
   from Biome). lightningcss compiles it; for `*.module.css` it renames classes deterministically to
   `msp-<component>-<class>` (the file name, with its folder when they differ; duplicates are an error) and wraps the
   result in `@layer msp.components`. The module exports the typed class map (`classes`) and registers the CSS with
   `defineStyles(id, css)`.
2. **Components import the generated module** (`import { classes } from './_generated/button.module.css.js'`), so a page
   carries exactly the CSS of the components it imports.
3. **Elements adopt the styles.** On connect, every `MspElement` calls `adoptStyles(document)`: one constructed
   stylesheet per document (a `<style>` element where constructed sheets are unavailable), updated when more modules
   register (lazily loaded views). Top-level elements get `.msp-root`, which scopes the reset and the tokens. No
   stylesheet to link, in any build setup; `styleText()` gives the whole CSS for a static file if a distribution wants
   one.
4. **Tooling**: `pnpm styles` regenerates, `pnpm styles:check` fails on missing, stale or orphaned modules (for CI),
   `pnpm styles:watch` regenerates on change, and Storybook runs the same generator as a Vite plugin. Jest maps
   `./_generated/x.css.js` to `x.css.ts`. Edit the `.css` file, never the generated module; Biome formats `.css`, so
   regenerate after formatting.

Layers, lowest first: `msp.reset`, `msp.theme`, `msp.components`, `msp.overrides`. Themes are set with `data-msp-theme`
(`light`, `dark`, `blue`, `auto`) on any element; derived tokens are recomputed wherever a theme is set, so a subtree
can use another theme.

## 4. Rules

- Components use tokens, not literal colours or sizes.
- No inline styles except values computed at runtime (positions, sizes from the shell).
- Class composition through a small `cx()` helper; state classes bound to `ReadonlyState`s via `classList`.
- Icons stay SVG components (the existing ~86 Material-style icons port as data), coloured with `currentColor`.
- The old `msp-*` global classes are not reused; ui-next and the old UI can coexist on one page.
