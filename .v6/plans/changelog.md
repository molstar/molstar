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
- The Prism viewer over a real plugin (`@molstar/prism-viewer`, `apps/prism-viewer`): one API whose operations every
  surface runs (commands, views, the URL, the agent), one undo history with checkpoints, built-in features as
  extensions (canvas settings, snapshots, export, measurements, superposition, state tree, loaders, remote states,
  animation), lazy where heavy; the viewer's agent on `@molstar/ai-session` (`agent/`), with live runs opt in
  (`pnpm test:live`); the session's record shown in History, kept in the browser on request; the assistant on the real
  agent (`apps/prism-viewer` passes it: the person's own keys, kept scrambled)
- Model Archive's quality in the Prism viewer (`@molstar/model-archive-extension/prism-viewer`): color by pLDDT or
  QMEAN (commands, `annotations.apply`), the pLDDT track in the sequence; sequence behaviors for extensions
- A colorblind-safe data palette in the Prism viewer (`canvas.palette`, `?palette=colorblind`): the categorical themes
  take Okabe–Ito, the core's new `okabe-ito` color list; the light restore after a reload; touch targets of 44px in
  compact mode
- Prism's high-contrast UI theme, and stronger borders and muted text in every theme when the system asks for more
  contrast
- Prism accessibility (gate G8): every generated form control is named by its label; the command menu is a modal
  combobox (the highlighted option as the active descendant, the result count said, focus kept and given back); cards,
  menus and closed panels give the focus back; menus, the canvas toolbar's menus, the picking level and the rail's tab
  lists take the arrows; muted, subtle, accent and status text at 4.5:1 or more in every theme (light's accent,
  success and warning darker, blue's accent, danger and info lighter); field and pressed edges in the high-contrast
  theme and under `prefers-contrast: more`; selections and foci by target, and the assistant's finished reply, said
  in the one live region; the sequence's keyboard cursor said and ringed; compact panel targets at 44px
- The Prism viewer's Storage view (Settings): what it keeps in the browser and where (kept records with size and age,
  preferences, the light restore, the assistant's keys by provider and scope), each with Forget, and Forget everything
- Mol*'s extensions in the Prism viewer, commands first: model and geometry export (model export now writes PDB too),
  video export and Zenodo (lazy), assembly symmetry, membrane orientation (ANVIL), Kinemage, and the data-only ones
  (validation reports and partial charges through `annotations.apply`); the single-file `prism-viewer.js`
- Prism viewer: the unit cell of a crystal structure as a row under it (show, hide, color, opacity; `structure.unit-cell`)
- Prism viewer: move and rotate a structure by hand (v5's Conformation Transform): the structure inspector's Transform and `structure.transform`, `structure.transform-reset`
- Prism viewer: `load.data` loads a file's contents (text, base64 or bytes) with a format; `load.url` and `load.data` load every model or an assembly, and the URL takes `all-models=1` and `assembly=<id>`
- Prism viewer: with remote states on (`ui.showRemoteStates`), Home lists the state server's states; a click opens one
- Prism viewer: after a load, Home offers what the structure carries beyond B-factors and secondary structure ("AF-… has pLDDT": Show applies it, one Undo)
- Prism viewer: the PAE plot of a ModelCIF model (model-archive extension): a heat map card, a click selects the two residues; `model-archive.pae` and "Show the PAE plot"
- Prism viewer: Help's basics (load, select, focus, the command menu, undo, the assistant) above the bindings
- Prism viewer: an Illumination toggle and, where WebXR is, an XR button on the canvas toolbar (`canvas.xr`).
- Prism viewer: an expanded viewport (`view.expand`: the canvas fills the viewer, Escape brings the panels back) and `view.fullscreen`.
- Prism viewer: the Export card's Camera section: the camera in numbers (as v5's), Copy MVS camera; `export.viewport`
- Prism viewer: the Selection tab's summary (what is selected, the structures and chains it spans, Focus, Clear; with nothing selected, a way into Select mode); the style layers applied (color, transparency, glow, material, clipping) listed under its style form, each edited or removed, Clear all (`representation.layers`, `.layer-update`, `.layer-remove`, `representation.reset` with `all`); the superpositions made listed under Superpose and removed (`superpose.list`, `superpose.remove`)
  and `load.viewport` (a viewport state as JSON); `camera.set` and a richer `camera.get`. Session & state is v5's State
  panel: State, Session, Open, and the Save options (`session.save`'s `options`). The Image section chooses Image or
  Animation (the `export.picture` slot, whose views get the shared frame, `ExportFrame`); `video.export` takes `axes`.
- Prism viewer: the plugin's toasts show in the status column.
- Prism viewer: the focus history (`focus.history`, `focus.again`, the "Focus again…" command), as v5's focus picker
  lists recent foci.
- Prism viewer: volume streaming's view, detail and channels (`volume.streaming-view`), with a "Volume streaming" row
  and inspector under the streamed structure.
- Prism viewer: a map's levels in the inspector of a volume, an isosurface and streamed density: the map's histogram
  with a marker per iso value to drag (one undo step a drag), the value in σ or map units, each level's color and
  opacity; an isosurface's tree row shows its color.
- Prism viewer: a direct volume's transfer function in the inspector: the opacity curve over the map's histogram (a
  log scale, Mol*'s presets), the area under it in the colors painted there, the volume-value theme's stops on the same
  axis (a stop on a point moves with it), Colors at peaks (a popover: the theme's colors or a palette; a preset moves
  the stops onto its peaks with their colors); previews while dragging
  are throttled, the final value applies on release, one undo step a drag.
- Prism viewer: a structure representation's quick settings above its whole form: its type and color theme as
  choices, one color from a row of swatches, the opacity and the size as sliders (a drag is one undo step); a map
  level's color from the same row.
- Prism viewer: a story's audio (Mol*'s markdown audio commands, MVS narration) has its own bar at the top of the
  playback card (play or pause, the position to click or drag, the time, stop, close) and a canvas button while there
  is a player, so it can be paused or stopped with the card shut; the playback card opens for a story loaded before
  its views came.
- Prism viewer: the quick type and color of a representation drawn as what they are (the types' icons in a segmented
  strip, the themes as tiles in their colors on the structure); pLDDT offered only for a prediction (ModelCIF quality
  metrics, or no experimental method); a component's inspector lists its representations to open; Contents shows its
  counts as figures and its entities as rows (kind, chains, atoms), each opened in the tree; Focus is one row for the
  scene; the Structure panel's header has no Component button (the command menu adds components).
- Prism UI: a select's options are grouped by Mol*'s categories, each group a fold (the current option's open); a
  form's select opens as a popover under its field; `CommandManager.open` takes values to open a command with.
- Prism viewer: a component's Add representation and Set coloring open their commands in the command menu for it,
  the component (or a preset component's keyword) filled in, and one Esc closes a command opened that way; Modify by
  selection acts at once (include, subtract or intersect the selection).
- Prism viewer: the playback card is the Story (a book icon), one frame, its text wrapping inside it; a description's
  first heading that only repeats the snapshot's title is not shown; ← and → step through a story (quiet where there is
  none); the Snapshots tab's Play is an icon beside its options.
- Prism viewer: a component's Modify by selection is three icons on its row (include, subtract, intersect), always
  there, unavailable with nothing selected; the inspector's eye is the size of the tree's; Export's first section says
  "Image & animation" once an extension adds the animation; a quick panel's sliders sit together.
- Prism UI: a number field without bounds is text with a decimal keyboard (a native number field showed the locale's
  decimal comma, "0,0"), and one without a step shows its value as it is instead of rounded to a whole number; a typed
  comma reads as the decimal point. A section can have a `heading` that follows what it offers.
- `PluginAnimationManager.stop` waits for a frame under way: applying a snapshot during another's transition no longer
  lands that transition's frame over the new snapshot (stepping through an MVS story showed the wrong scene).
- Prism viewer: the playback card is quieter: the snapshots row has no title (the position says it on hover and lists
  them all), the snapshot's transition scrubber is behind a toggle in it, and a story's audio has a mute there, the
  Story card shut or not; stepping through a story posts no status message (said to screen readers; an operation's
  `quiet`).
- Prism viewer: SeqViewSpec archives (`.svsx`) load by URL, as a file or dropped: the document, the structure the
  archive holds (opened as a file), and its relative URIs resolved inside it.
- Prism viewer: saving a session can include its record (the form's Record, off by default: prompts and
  replies are in it), as `record.jsonl` in a `.molx` or a `record` field in a `.molj`; opening such a file shows the
  record in History, read-only.
- `PluginStateSnapshotManager.serialize` takes `extra` (files in a `.molx` beside `state.json`, top-level fields in a
  `.molj`), and `open` returns what it did not use, so a host can keep its own data in a session file.
- Prism viewer: Home lists the records kept in earlier sessions; one opens in History, read-only, beside what is
  loaded ("This session" goes back); each can be forgotten there, and this session's kept copy from History. At most
  `ui.keptRecordsMax` (20) are kept: keeping one more names the oldest and asks before forgetting it.
- Prism viewer: the Story card stays open on a click outside it (Esc or its ✕ close it) and its description fills the
  card; a first heading is left out when the title already says it (a part counter like "[1/2]" aside), whatever order
  the title and description change in. Overlays can be `persistent` and `fill`.
- The camera: a transition under way when the animation loop's clock starts again (after an image is drawn, e.g. a
  snapshot's thumbnail) keeps going instead of flinging the camera far out (a negative or huge radius, an empty view);
  a state snapshot keeps where an animating camera is going (`Camera.getSettledSnapshot`) and the view a pending camera
  reset frames (after a load), not the camera before it.
- Prism viewer: an inspector's header has the item's eye at its end (Show or Hide; Turn on or off for the unit cell,
  the streaming and an extension's row), in place of a Turn on row while it is off.
- Prism viewer: the canvas's playback card has a row each for the snapshots (the current one, its position and the
  list, previous, play, next), the current snapshot's transition (a scrubber, play from the start or stop:
  `snapshot.seek`, `snapshot.transition`) and the trajectory's models; extensions add rows in the `playback` slot.
- Prism viewer: the sequence's No colors is plain letters (a neutral strip; plain text in the text style) instead of
  Seq*'s default residue colors.
- Prism UI: a form's presets open as a popover under their button (`PresetsList`), and a list inline in a panel lets
  the wheel scroll the panel at its ends (`ScrollArea`'s `chain`).
- `ColorTheme.Palette` (the volume-value theme): a color list with offsets interpolates between its stops (as
  `ColorScale` does) instead of switching halfway between them, and offsets with many decimals no longer make a huge
  palette.
- Prism viewer: load a SeqViewSpec into the sequence (`load.svs`, `.svsj` files, `svs-url`): it shows over the loaded
  structure it names (loaded if missing) until another structure is focused or loaded, and Undo and sessions keep it.
- Prism viewer: `load.mvs` takes v5's options: `append` (keep the Snapshots list), `keepCameraOrientation`,
  `snapshot`, `extensions`.
- Prism viewer: snapshot descriptions run Mol*'s markdown commands (focus, highlight, camera, snapshots, audio) and
  show MVSX images, color swatches and tables (`snapshot.link`; prism-ui's `Markdown` takes a host).
- Prism viewer: the playback card opens by itself when a loaded story has snapshots, and shows a single snapshot's
  description.
- Prism viewer: the image export's crop box can be dragged and resized (and moved with the arrow keys).
- Prism viewer: Help lists the canvas keys (view, moving, snapshots).
- Prism viewer app: the `map-provider`, `show-toggle-fullscreen`, `debug-mode` and `timing-mode` URL parameters are read.
- Prism viewer: the structure inspector sets the model of a trajectory and the structure type (the model, an assembly,
  symmetry mates or unit cells); the structure row's menu has Preset, Assembly, Crystal symmetry and Stream density.
- Prism viewer: extensions' rows in the scene tree (`ViewerExtensionUI.rows`, drawn with `SceneRow`); assembly symmetry
  and ANVIL's membrane each have a row under the structures they apply to.
- Prism viewer: the `<prism-viewer>` element (`@molstar/prism-viewer/element`; `spec`, `pdb`, `theme`,
  `hide-controls`), defined by the single file.
- The Prism viewer loads its assistant on first use (the Assistant tab, the Ask row, or `?assistant`), not at start; the
  command menu renders its options as they scroll in; History reads only what is new
- `@molstar/ai-tools`: a dynamic toolkit (`toolkit.add`, lazy skills disclosed during a turn); `ai-session` reads it as
  it is
- Prism viewer: Seq* is its own chunk (`sequence.runtime` in the spec), loaded when a structure is first shown: the
  app's first load is 1.25 MB gzipped (from 1.68)
- Prism viewer: the assistant's Inspector (Trace, Tools, Limits) over the real agent, opened from the Assistant tab's
  Inspect rows, the turn card's trace and the Assistant trace/tools/limits commands; extensions' views may add `pages`
- Prism viewer: the params form lays params out as v5 does (a choice's own params behind its settings button, an
  Advanced fold, a fold per category); the scene tree's focus group (always under a structure, holding the focus
  target and surroundings); Preset and Styles open the command menu in their groups; cards stay open while the canvas
  is turned; a lazy extension that fails to load says so and offers Retry
- The snapshot manager's `getStateSnapshot` and `serialize` take `keepList` (the current scene without changing the
  list)
- State actions and transformers take `display.triggerLabel`, the label of the button that runs them ("Download",
  "Open"); a UI falls back to a generic one without it
- `PluginState.registerCustomState(key, { get, set })` lets a component keep its own state in plugin snapshots
  (`snapshot.custom[key]`): saved by `getSnapshot`, restored by `setSnapshot` after the data tree (`undefined` when the
  snapshot has no value), carried by snapshot manager entries and `.molx` / `.molj` sessions
- Every plugin path that changes the canvas props emits `events.canvas3d.settingsUpdated` (initial props, `clear(true)`,
  snapshot, transition-frame and snapshot-animation restore, canvas context props); new
  `PluginContext.setCanvas3DProps` and `setCanvas3DContextProps` set the props and emit. A direct `canvas3d.setProps`
  still emits nothing

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

## MessagePack, UTF-8 and zip hardening (Seq\* port)

- MessagePack (`@molstar/io/common/msgpack/encode`, `decode`): extension types. `encodeMsgPack(value, ext?)` takes a
  `MsgPackExtEncoder`, `decodeMsgPack(buffer, ext?)` a `MsgPackExtDecoder`; without a decoder an extension value reads
  as `{ type, data }`, and `skipMsgPackValue` skips them
- MessagePack encoding follows JSON's rules: a number that is not finite is written as nil (it threw), an integer
  outside the int32 and uint32 ranges as a float 64 (it threw), and the keys of an object whose value is `undefined`
  are left out (they were written as nil). BinaryCIF files no longer carry a nil `mask` key
- MessagePack decoding of untrusted input: a cut buffer throws (it read as data), a length header that claims more
  items than the bytes left throws, arrays grow as items are read rather than at the claimed length, and a
  `__proto__` key is an own key, as `JSON.parse` makes it. On a 340 KB `.bcif` the checks cost about 3% of a decode
- Fix UTF-8 encoding of characters outside the Basic Multilingual Plane (`utf8Write`, `utf8ByteCount` in
  `@molstar/core/util/utf8`): a surrogate pair was written as two 3-byte sequences instead of one 4-byte sequence;
  ASCII and BMP output is unchanged
- Zip hardening (`@molstar/core/util/zip/zip`, a security fix): `unzip` looks for the end record only in the last
  64 KB (a truncated archive looped for ever), reads the central directory first and takes sizes, name and CRC-32 from
  it (archives written with a data descriptor, as macOS Finder writes them, have zeros in the local header), checks that
  no local entry is listed twice, and that every entry has its declared size and CRC-32; `inflateRaw` and `inflate` with
  a buffer stop at its end and return what they inflated. New optional `UnzipLimits` (`maxEntries`, `maxBytes`) as the
  fourth argument of `unzip`. Invalid or damaged archives that were read before now throw

## Seq\* input streams

- The Seq\* viewer's input for consumers is `ctx.interaction`: rxjs streams (`hover`, `press`, `click`, `drag`, `menu`,
  `wheel`, `key`, `leave`, `activate`) of fresh payloads in Mol\*'s button and modifier terms, scoped by `on`, built only
  while observed; the behaviors' base events are `ctx.input.events`. The track menu is `interaction.menu`, the browser's
  context-menu request (the `contextMenu` behavior, `events.menu` and `events.activate` are removed). Behaviors list
  their bindings as Mol\*'s `Binding` (`ctx.behaviors.bindings()`)

## Seq\* line styles

- A `ranges` layer's `line_style` channel (format 1.0): `uniform` with a `line_style` or `field_line_style` over a column
  of the rows draws each item `solid`, `dashed`, `dotted`, `wavy`, `double`, `box` or `highlight`, in tracks, text and
  SVG exports; `bookmark.add` takes an optional `line_style`. `PathBatch` strokes may be dashed (`dash`, `gap`,
  `phase`), and `Paint2D` includes `setLineDash` and `lineDashOffset`.

## Seq\* selection

- [Breaking] The `select` behavior's `shiftPositions` defaults to `toggle`: Shift with a click or a drag on positions
  adds or removes them (it extended a range from the anchor); `extend` stays a param

## Asset manager byte cache

- `AssetManager.resolve` reads URL assets through a URL-keyed, ref-counted byte cache: a URL (and its headers) is
  downloaded once whoever asks first, concurrent requests share the download, the bytes are held while any resolved
  asset is held and dropped with the last `dispose`, and a failed download is not cached. Requests with a body and
  `store: false` bypass it
