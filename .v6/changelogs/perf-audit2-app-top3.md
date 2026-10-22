## Performance: perf audit 2, app-level (startup, Prism start code, per-change walks)

- Startup: `ElementSymbol` memoizes the symbols it reads (at most 4096) instead of a table of 242k strings built when
  `@molstar/model` loads: every app, test worker and future worker starts with 22 MB less heap (35.3 → 13.3 MB after
  importing `PluginContext` + `DefaultPluginSpec` in Node), and that import is ~2× faster
- Prism: AlphaFold DB's PAE reader is its own module (`@molstar/prism-viewer/sequence/molstar/pae-matrix`), so the
  model-archive extension no longer pulls all of Seq\* into the page's start code: 6.31 → 5.16 MB raw, 1.96 → 1.55 MB
  gzip loaded before the viewer runs
- Build: an ESM app's entry budget (`molstar.budgets`) counts what loads at start (the entry and every chunk it imports
  statically), not the 7 KB entry stub; Prism's is 5.5 MB, so a static import of a lazy part fails the build
- Plugin: `StructureInfo` visits the root models only when a root model is created or removed, not for every
  component; the measurement manager keeps its group from cell events instead of searching the whole tree on every
  state change. 500 structures: showing them all ~6 % faster; changing the theme, overpaint or removing one structure
  58 → 46, 57 → 37, 52 → 37 ms
