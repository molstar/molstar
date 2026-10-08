# v6 changelog draft

Entries for the v6 release notes. `CHANGELOG.md` keeps fixes that also apply to 5.x; changes that only exist in v6 go
here, and the v6 `CHANGELOG.md` section will be written from this list. Detailed API changes and migration notes are in
the [API change ledger](breaking-v6-changes.md).

## Workspace and tooling

- [Breaking] Mol\* is published as versioned ESM workspace packages (`@molstar/core`, `@molstar/model`,
  `@molstar/plugin`, ...) instead of `molstar/lib/...`; there is no CommonJS build (#1951)
- TypeScript 7 builds, Biome lint and formatting (#1958)

## Plugin composition (#1961)

- [Breaking] Plugin registries start empty and a plugin registers only what `PluginSpec.registry` lists.
  `DefaultRegistry` (`@molstar/plugin/default-registry`), `DefaultPluginSpec()` and `DefaultPluginUISpec()` provide the
  full built-in set; `Viewer.create` and its options are unchanged
- [Breaking] `PluginSpec.actions`, `animations` and `customFormats` are removed; they are registry entries now
- [Breaking] `StateTransforms` is removed from the library; transformers, formats, presets and selection queries are
  split into modules by functionality. The classic Viewer global keeps `molstar.lib.plugin.StateTransforms`
- [Breaking] `DefaultPluginSpec` and `DefaultPluginUISpec` move to `@molstar/plugin/default-spec` and
  `@molstar/plugin-ui/default-spec`; `createPluginUI` and the view models require a spec
- `plugin.register(entry)` registers an entry at run time and returns an undo; registries count registrations by
  provider and reject a different provider under a registered name
- Presets resolve by id or an optional `alias`; the built-in short names (`'default'`, `'auto'`, ...) are aliases
- Script languages other than MolScript are enabled by importing `@molstar/model/script/transpilers/<lang>` (or `all`)
- Snapshot loading checks all transformer ids before changing state; unregistered representation and theme names fall
  back to the registry default with a warning
- `@molstar/plugin-ui/default-ui` has the registry-free UI defaults that `DefaultPluginUISpec()` is composed from
- An empty theme or format registry points to the catalog or `DefaultRegistry` in its warning or error
- The `external-structure` and `external-volume` color themes are registered again, as the `ExternalColorThemes` entry
  (`@molstar/plugin/themes/external`)
- mesoscale-explorer no longer bundles the full built-in registry (4.88 MB to 4.21 MB)
- New `examples/slim-plugin` (SDF and ball-and-stick only)
- Example data files move from `examples/` to `data/examples/`
