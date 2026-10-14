# Copied: Mol\* Prism UI's tokens and reset

A local copy, so the Seq\* Storybook has Prism's look before Seq\* moves into the Mol\* v6 repository. At the move
this folder is deleted and the Storybook imports Prism's own styles (`@molstar/prism-ui`'s `styles/base`).

- **Origin**: Mol\* v6, `packages/prism/ui/src/styles/theme.css` and `reset.css`, at commit `10cde25a3` (2026-10-11).
- **License**: MIT, Mol\* contributors; the files keep their header.
- **Unchanged.** `prism.css` (ours) declares Prism's layer order (as `styles/registry.ts` does) and imports both. The
  tokens are `--msp-*` on `:root` and on every `[data-msp-theme]` (`light`, `dark`, `blue`, `auto`); the reset applies
  only under `.msp-root`, so pages without it are untouched.
