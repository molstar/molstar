# Prism

The new Mol\* UI. The repository's rules ([../../AGENTS.md](../../AGENTS.md)) apply; these are Prism's own. Designs:
[.v6/prism-ui](../../.v6/prism-ui/README.md) (start with its README; `design/contexts.md` holds the contracts).

- Packages: `ui` (the common base: its own DOM JSX runtime, components, the shell, the command menu), `viewer` (the
  viewer mode, Mol\*'s own UI), `stories` (the story mode), `assistant` (the AI assistant's UI). The modes and the
  assistant are opt-in: an app gets one only by depending on it, and `ui` never imports them
  (`ui/src/_test/opt-in-boundary.test.ts`).
- Not React: views are `.tsx` on the runtime in `ui/src/dom` (`Show`, `For`, `observableState`); state lives in models,
  views bind to their observables.
- Styling: CSS modules compiled to `_generated/` by `styles`; every color and spacing from `--msp-*` tokens (the scale
  is `--msp-space-1..6`; `styles:check` rejects other values). The look is flat and quiet: few borders and shadows, the
  accent only for actions and the picked item, text that does not fit fades rather than ends in an ellipsis. Check light
  and dark, wide and narrow, in Storybook.
- **Snappy is a requirement** (David, 2026-10-14; the rules in [styling.md §5](../../.v6/prism-ui/design/styling.md)):
  no layout reads after writes in the same frame, pointer-rate paths that only update what changed (once per frame,
  nothing rebuilt or allocated per event), animation of `transform` and `opacity` only, flat single-class selectors,
  contained and virtualized subtrees. Measure an interaction in the browser's performance panel before calling it done.
- **Hidden UI is inert and does no UI updates** (David, 2026-10-14;
  [styling.md §5](../../.v6/prism-ui/design/styling.md)): what the shell hides but keeps mounted (a folded section, a
  panel behind an overlay card) is `inert` (no hover, tooltip, focus or screen reader) and its scope is paused, so its
  bindings hold no subscriptions until it shows again (`activity()`, `Scope.setActive`; the runtime does it once, views
  do nothing). A view that subscribes, observes or loops by itself pauses that with `onActive(active => …)`; content
  that is hidden by removing it needs none of this.

## Commands

`pnpm -C storybook storybook` or `mise run storybook` (port 6006; the top-level `storybook/` shows the stories of the
packages registered in `storybook/packages.ts`, every Prism package among them; stories go in `_stories/` folders),
`pnpm -C packages/prism/<pkg> styles` after editing a CSS module (`styles:check` verifies the generated files and the
spacing scale), `typecheck:storybook` for stories and tests.
