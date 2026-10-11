# Styling

## 1. Decisions

- **CSS Modules** for component styles (`button.module.css` next to `button.tsx`): local class names, dead-class
  detection, no global naming scheme to police. Replaces SCSS; native CSS nesting, custom properties, `color-mix()`, and
  `@layer` cover what Sass mixins and colour functions do today.
- **Design tokens as custom properties with the `--msp-` prefix**, defined once in `styles/tokens.css`.
- **Themes are token sets** switched at runtime by `data-theme` on the shell root (`light`, `dark`, `blue` to keep the
  current variants, `auto` following `prefers-color-scheme`). A host can override any token on the element.
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

`tsc` does not process CSS, so library output needs one extra step:

1. **Compile CSS Modules** with deterministic, readable names (`msp-<file>-<class>`, e.g. `msp-button-primary`).
   lightningcss supports CSS Modules with a configurable naming `pattern`, nesting lowering, and minification in one
   native dependency. esbuild's `local-css` loader also handles `.module.css` but its naming is not configurable (to be
   verified during setup); deterministic names matter because hosts and tests may read them and because library and app
   builds must agree.
2. **Emit per module** a `*.module.css.js` (the class map) and a `*.module.css.d.ts` (typed keys), next to the compiled
   component in `lib/`; TS sources import `./button.module.css` and resolve to these at build and type-check time.
3. **Concatenate** all modules plus tokens into `lib/ui-next.css`. Apps import that one file
   (`@molstar/plugin-ui-next/ui-next.css`); app builds can instead let esbuild bundle the modules from source.

This is added to the shared workspace scripts (`scripts/workspace/assets.mjs` and `scripts/esbuild/app.mjs`) and to the
browser distribution. Biome formats and lints CSS.

## 4. Rules

- Components use tokens, not literal colours or sizes.
- No inline styles except values computed at runtime (positions, sizes from the shell).
- Class composition through a small `cx()` helper; state classes bound to `ReadonlyState`s via `classList`.
- Icons stay SVG components (the existing ~86 Material-style icons port as data), coloured with `currentColor`.
- The old `msp-*` global classes are not reused; ui-next and the old UI can coexist on one page.
