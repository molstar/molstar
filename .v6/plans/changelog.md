# v6 changelog draft

Entries for the v6 release notes. `CHANGELOG.md` keeps fixes that also apply to 5.x; changes that only exist in v6 go
here, and the v6 `CHANGELOG.md` section will be written from this list. Detailed API changes and migration notes are in
the [API change ledger](breaking-v6-changes.md).

## Workspace and tooling

- [Breaking] Mol\* is published as versioned ESM workspace packages (`@molstar/core`, `@molstar/model`,
  `@molstar/plugin`, ...) instead of `molstar/lib/...`; there is no CommonJS build (#1951)
- TypeScript 7 builds, Biome lint and formatting (#1958)
- Node 22.13 or later
- The documentation is an Astro Starlight site (`apps/docs`) that keeps the v5 docs (`/docs/v5/`) next to the v6 docs
  (`/docs/v6/`); `/docs/` redirects to the default version and the old page URLs redirect to v5

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
- Script languages other than MolScript are enabled by importing `@molstar/query-language/transpilers/<lang>` (or `all`)
- Snapshot loading checks all transformer ids before changing state; unregistered representation and theme names fall
  back to the registry default with a warning
- `@molstar/plugin-ui/default-ui` has the registry-free UI defaults that `DefaultPluginUISpec()` is composed from
- An empty theme or format registry points to the catalog or `DefaultRegistry` in its warning or error
- The `external-structure` and `external-volume` color themes are registered again, as the `ExternalColorThemes` entry
  (`@molstar/plugin/themes/external`)
- mesoscale-explorer no longer bundles the full built-in registry (4.88 MB to 4.21 MB)
- New `examples/slim-plugin` (SDF and ball-and-stick only)
- Example data files move from `examples/` to `data/examples/`

## Standalone MolQL and MVS authoring

- [Breaking] Language builders, symbol tables, text parsers, and transpilers move to `@molstar/query-language`;
  executable query compilation and molecular `Script` operations stay in model
- `@molstar/mvs-builder/molql` exposes `MolScriptBuilder` and `compileScript` for MolScript, PyMOL, VMD, and Jmol,
  without loading the molecular query runtime or requiring plugin language registration
- MVS standalone validation and `mvs-validate` reject malformed MolQL, unknown callable names, and invalid argument
  names using symbol-definition tables; custom vocabularies use an explicit lookup
- `mvs-validate` reports file/expression paths and continues checking later inputs after malformed JSON or file errors
- [Breaking] Remove the duplicate `@molstar/mvs-builder/expression` module; use
  `@molstar/query-language/language/expression`

## MVS selectors through model indices

- New MolQL generator `element-schema` (`structure.generator.element-schema`, `sel.atom.element-schema` in MolScript)
  with the data symbols `schema-item` (`structure.type.schema-item`) and `empty` (`core.type.empty`, `empty-value`): it
  resolves tables of `StructureElement.SchemaItem` rows (`label_asym_id`, `auth_seq_id`, `beg_label_seq_id`, ...)
  through the indices of the model instead of testing every atom for every row, as a single set. 1000 rows on a
  720k-atom structure take about 17 ms (seconds before), 10,000 rows about 130 ms. A table can be given as rows
  (`schema-item` for each) or by columns (one `schema-item` whose fields are lists); the compiler passes a list through
  as it is, so the cost of the expression does not grow with its size, and `Schema.toExpression` switches to columns
  from 10 rows
- `null` and lists of values (of `string`, `number`, `boolean` and `null`) are literals in the arguments of an
  expression (`null`, which is what JSON makes of `undefined` in a list, is a missing value)
- `StructureElement.Schema.toExpression`, `toLoci` and `toBundle` (`Loci.fromSchema`, `Bundle.fromSchema`) use it, and
  MVS does for the selectors given as rows of `component` nodes (the plugin's own components, which the classic UI can
  modify), inline tooltips and color layers. Annotation components, labels, annotation colors and tooltips, and
  primitive positions resolve rows through the same indices directly
- Code which reuses buffers between calls registers them with `ScratchMemory` (`@molstar/core/util/int-buffer`) under a
  context (`SchemaSelectionScratch` for the rows), which keeps what is small after each call, and
  `ScratchMemory.release(...contexts)` frees it (every context when none is named); MVS calls it when a load is done
- Rows are resolved the same way everywhere: a residue without `label_seq_id` (ligands, waters) no longer matches a
  `label_seq_id` range or `end_label_seq_id` as MolScript did (it read a missing value as 0), `atom_id` rows select
  every atom with that id, `operator_name` is supported in rows, JSON annotation values of the wrong type (`"10"` for a
  number) are converted (see the API change ledger), and a range of `auth_seq_id` with
  `beg_auth_ins_code`/`end_auth_ins_code` is a range of residues in the order of the input file (new in `SchemaItem`,
  not in MolViewSpec yet)
- Components, layers and labels made from rows are named by what they select in a compact form (`A`, `A432`, `A:10-50`,
  `A, B, C`, `120 residues`) instead of a dump of the rows
- Annotation rows use reusable typed-array buffers, a sort-free range union, and paint each element once (later rows
  win, and the pass stops when every element has a row), so overlapping rows no longer cost the sum of their sizes: the
  union of 20,000 rows allocates 26 MB instead of 1.4 GB, and 20,000 overlapping rows take about 3 ms instead of 1.5 s
- The row indices are shared by the frames of a trajectory (10 frames of 720k atoms, 100 rows: about 800 ms to 150 ms),
  and rows without chain criteria (`residue_index`, insertion code, `comp_id`) use whole-model lookups instead of a
  lookup in every chain

## Mol\* UI v-next, phase 1

- Prism, the new UI (`packages/prism/*`), and the AI agent packages (`packages/ai/*`): private, not published, mocked in
  Storybook (`storybook/`, the repository's Storybook for the registered packages)
- State actions and transformers take `display.triggerLabel`, the label of the button that runs them ("Download",
  "Open"); a UI falls back to a generic one without it

## Video export

- [Breaking] Rename the MP4 export extension/package, public APIs and viewer option to video export; retain its
  serialized behavior ID
- Add a video renderer backend interface and make Mediabunny the default, with runtime-generated codec/container
  parameters
- Support all Mediabunny video codecs when an encoder is available, with quality/bitrate, FPS and advanced codec
  settings
- Native PD group presets provide Balanced (default), Fast preview, High quality and Smooth motion video settings;
  presets preserve the selected format and codec, with Automatic hardware selection; unsupported FPS/resolution disables
  rendering instead of silently replacing the codec
- Reject invalid API settings, accept fractional Mbps bitrates, and discover codecs supported outside the High preset
- Reuse codec discovery across formats and regenerate UI schemas only for resolution, FPS or format changes
- Keep H.264 WASM available through an explicit backend import for headless use; the MVS CLI continues to use it
- [Breaking] Encoded video results include data, MIME type and extension; backend settings replace top-level
  quantization
- Fix encoder cleanup on failure/cancellation, the extra terminal frame and swallowed headless file-write errors
