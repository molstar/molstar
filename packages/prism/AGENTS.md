# Prism

The new Mol\* UI. The repository's rules ([../../AGENTS.md](../../AGENTS.md)) apply; these are Prism's own. Designs:
[.v6/prism-ui](../../.v6/prism-ui/README.md) (start with its README; `design/contexts.md` holds the contracts).

- Packages: `ui` (the common base: its own DOM JSX runtime, components, the shell, the command menu, Storybook),
  `viewer` (the viewer mode, Mol\*'s own UI), `stories` (the story mode), `assistant` (the AI assistant's UI). The modes
  and the assistant are opt-in: an app gets one only by depending on it, and `ui` never imports them
  (`ui/src/_test/opt-in-boundary.test.ts`).
- Not React: views are `.tsx` on the runtime in `ui/src/dom` (`Show`, `For`, `observableState`); state lives in models,
  views bind to their observables.
- Styling: CSS modules compiled to `_generated/` by `styles`; every color and spacing from `--msp-*` tokens (the scale
  is `--msp-space-1..6`; `styles:check` rejects other values). The look is flat and quiet: few borders and shadows, the
  accent only for actions and the picked item, text that does not fit fades rather than ends in an ellipsis. Check light
  and dark, wide and narrow, in Storybook.

## Commands

`pnpm -C packages/prism/ui storybook` (port 6006; it shows the stories of the other Prism packages too),
`pnpm -C packages/prism/<pkg> styles` after editing a CSS module (`styles:check` verifies the generated files and the
spacing scale), `typecheck:storybook` for stories and tests.
