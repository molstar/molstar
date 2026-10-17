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
- `pnpm deploy:local` and `pnpm deploy:remote` are back: they build the workspace and prepare (and, remotely, push) the
  GitHub Pages folder, which now also carries the Prism viewer (`prism/`) and the UI Storybook (`ui-storybook/`);
  `--<app>-prefix <folder>` moves an app (`pnpm deploy:remote -- --prism-prefix prism-preview`)

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

- Prism, the new UI (`packages/prism/*`), and the AI agent (`packages/ai`): private, not published, mocked in Storybook
  (`storybook/`, the repository's Storybook for the registered packages)
- The AI packages are one package, `@molstar/ai` (`packages/ai`), with explicit entries (`.`, `./harness`,
  `./providers`, `./providers/{anthropic,openai,google}`, `./models`, `./trace`, `./testing`, `./narration`): the nine
  `@molstar/ai-*` packages are folders of it, and Biome enforces the layering between them
- The AI model adapters are providers, a thin layer over the AI SDK (`createAnthropicProvider`, `createOpenAIProvider`,
  `createGoogleProvider`; the port `Provider`, `ProviderBindings`; a session's `models.providerFor(id)`): a request is
  `streamText`'s own arguments, retained as JSON and sent as it is; an inline attachment's `digest` is
  `murmur3-128:<hex>`
- The AI transcripts are called threads, not fibers (`createThreadHost`, `ThreadId`, `@molstar/ai/testing/threads`);
  records a developer kept in the browser under the old stored names (table `fiber`, scope `fiber:`) are unreadable
- New `examples/minimal-prism-viewer`: the smallest Prism viewer composition (an mmCIF structure as cartoon or
  ball-and-stick, no features, extensions, sequence or assistant), a Vite app; its README lists what it leaves out and
  its bundle size
- The Prism viewer over a real plugin (`@molstar/prism-viewer`, `apps/prism-viewer`): one API whose operations every
  surface runs (commands, views, the URL, the agent), one undo history with checkpoints, built-in features as extensions
  (canvas settings, snapshots, export, measurements, superposition, state tree, loaders, remote states, animation), lazy
  where heavy; the viewer's agent on `@molstar/ai-session` (`agent/`), with live runs opt in (`pnpm test:live`); the
  session's record shown in History, kept in the browser on request; the assistant on the real agent
  (`apps/prism-viewer` passes it: the person's own keys, kept scrambled)
- Model Archive's quality in the Prism viewer (`@molstar/model-archive-extension/prism-viewer`): color by pLDDT or QMEAN
  (commands, `annotations.apply`), the pLDDT track in the sequence; sequence behaviors for extensions
- Prism's own data palette, the default (`representation.defaultPalette`, the option and the operation (was
  `canvas.palette`): `prism`, `classic` for Mol*'s own Many-Distinct, `colorblind`; `?palette=`): the categorical themes
  take the core's new `many-distinct-prism` list (Many-Distinct redrawn in OKLCH, Mol*'s green and orange first); a spec
  gives the lists (`spec.palettes`: Prism's and the colorblind-safe one, hex strings or numbers); a session without a
  palette opens on the classic one
- A colorblind-safe data palette in the Prism viewer (`canvas.palette`, `?palette=colorblind`): the categorical themes
  take Okabe–Ito, the core's new `okabe-ito` color list; the light restore after a reload; touch targets of 44px in
  compact mode
- Prism's high-contrast UI theme, and stronger borders and muted text in every theme when the system asks for more
  contrast
- Prism accessibility (gate G8): every generated form control is named by its label; the command menu is a modal
  combobox (the highlighted option as the active descendant, the result count said, focus kept and given back); cards,
  menus and closed panels give the focus back; menus, the canvas toolbar's menus, the picking level and the rail's tab
  lists take the arrows; muted, subtle, accent and status text at 4.5:1 or more in every theme (light's accent, success
  and warning darker, blue's accent, danger and info lighter); field and pressed edges in the high-contrast theme and
  under `prefers-contrast: more`; selections and foci by target, and the assistant's finished reply, said in the one
  live region; the sequence's keyboard cursor said and ringed; compact panel targets at 44px
- Prism UI runtime: a plain-JS reactive base, rxjs only at its boundary: a state is its latest value plus listeners
  (`observableState`, `State`), `derive` is a lazy shared derived value, `toState` turns a stream into a state at a
  model's edge; components render as a function of state: `watch(sources, render)` writes what changed, once per task
  (one microtask scheduler, parents first; `flush()` renders now), `dom/writes.ts` writes only what differs,
  `keyedSelector(...).key(k)` is a keyed source; the trees' rows, the command menu's rows and the params form's controls
  use it (a scene-tree row 16 → 3.4 KB), views no longer pipe states (308 sites became `derive`), a guard test keeps new
  rxjs out of views, and a virtual list no longer lays out its new rows inside the scroll event
- Prism UI runtime: work wanted by key (`workQueue`: the last want of a key wins, run once in its lane, task, frame or
  idle, by priority, an optional budget carrying the rest over) and `latest` (each call aborts the one before). The
  scene tree flattens once a burst (expanding 1,000 components at once: 401 → 19.5 ms in jsdom); the command menu
  filters once a frame, for the latest text; the virtual lists keep rows made ahead of a fast wheel or flick, by its
  speed (1,200 px wheel steps: 2 scrolls of 40 onto the skeleton, from 40 of 40)
- Prism viewer: the scene tree follows the scene's structure as one change stream (`SceneModel.changes`, once per Mol*
  transaction) and flattens again only what changed: with 1,000 components open, 454 live listeners instead of 8,457
  (the same for any number), a mount 16 ms instead of 25 and 5.7 MB kept instead of 23, a toggle 60 µs instead of
  765; the state tree's toggle no longer walks the whole tree. Guards keep every view's subscriptions to what is shown
- Seq\*: a section surface pauses on its own (`Surface.setActive`, active by default): paused, it asks for no frame
  (pans, hovers and resizes only gather; a pending frame is cancelled), and one frame on resume draws all that waited.
  The Prism viewer's sequence band pauses its surface with the band's scope while hidden
- The Prism viewer's Storage view (Settings): what it keeps in the browser and where (kept records with size and age,
  preferences, the light restore, the assistant's keys by provider and scope), each with Forget, and Forget everything
- Mol*'s extensions in the Prism viewer, commands first: model and geometry export (model export now writes PDB too),
  video export and Zenodo (lazy), assembly symmetry, membrane orientation (ANVIL), Kinemage, and the data-only ones
  (validation reports and partial charges through `annotations.apply`); the single-file `prism-viewer.js`
- Prism viewer: the unit cell of a crystal structure as a row under it (show, hide, color, opacity;
  `structure.unit-cell`)
- Prism viewer: move and rotate a structure by hand (v5's Conformation Transform): the structure inspector's Transform
  and `structure.transform`, `structure.transform-reset`
- Prism viewer: `load.data` loads a file's contents (text, base64 or bytes) with a format; `load.url` and `load.data`
  load every model or an assembly, and the URL takes `all-models=1` and `assembly=<id>`
- Prism viewer: with remote states on (`ui.showRemoteStates`), Home lists the state server's states; a click opens one
- Prism viewer agent: `agent.skills` and `agent.core` add to or remove from the defaults by id (`{ add, remove }`),
  declared in the dependency-free `agent/skills` (a spec that extends them no longer brings zod);
  `assistant({ provider, model })` from `@molstar/prism-assistant/spec` is the lazy loader;
  `ViewerAgent.resume({ decisions, hooks })` and `AgentSession.resume({ decisions, hooks })` take options, and
  `ViewerAgent` exposes `busy` and a read-only `inspect` instead of its session
- Prism viewer: after a load, Home offers what the structure carries beyond B-factors and secondary structure ("AF-… has
  pLDDT": Show applies it, one Undo)
- Prism viewer: the PAE plot of a ModelCIF model (model-archive extension), as AlphaFold DB draws it: residue ticks,
  chain bands for a complex, the cell under the pointer lit in 3D and the sequence, a click selects its two residues, a
  dragged box (moved, resized by its edges and corners, cleared by its ✕) marks its aligned and scored ranges in 3D and
  the sequence (the rest greyed) with its mean error; a pLDDT strip along both axes (shown; its toggle hides it); the
  card from the toolbar and from a "Predicted aligned error" row in the Structure panel (an extension's row can open a
  card: `SceneRow` `onOpen`); `model-archive.pae` and "Show the PAE plot"
- Prism viewer app: remote states (Save and Open on the state server, Home's section and the card) are on (no rail
  button: Home's section and the command menu open them; uploading the current state is a section of the Export card)
- Prism viewer: volume streaming starts in the Auto view (the density around the whole structure, around the selection
  while there is one), X-ray maps too; a row's eye that turns something on (Volume streaming, the unit cell, an
  extension's row with controls) picks it, so its controls are in the inspector at once (the inspector's eye on the
  "off" row no longer leaves it on a row that is gone)
- Prism viewer on phones: the splash runs to the canvas's right edge (it was cut short by a rail that is not there), and
  the status messages and the selection stay on the canvas as on a desktop (the selection beside the toolbar, the
  messages floating above), so neither coming or going resizes the canvas or moves its toolbar; the floating sequence
  band spans the width
- Storybook: Viewer / Dynamic PDB: dynamicpdb.com's entry dpdb_gijqn2no (3K7F, deposited and qFit) through an extension
  written as a host would (`packages/prism/viewer/src/_stories/dynamic-pdb`, the API's responses as fixtures: its API
  answers only its own origin): `dynamic-pdb.load`, `.entry`, `.track`; the waters removed; a card with the entry, its
  mutations (hover lights them, a click selects), each model's metrics, eye and density (the volume server); sequence
  tracks for the mutations, B-factor, RSCC, qFit's alternate conformers and their Cα spread, each with its eye in the
  card (some hidden at first); the side panel pinned
- Prism viewer: `@molstar/prism-viewer` is the embedder's entry (`createPrismViewer`, `DefaultPrismViewerSpec`, the
  spec, viewer and API types); its subpaths are a curated list (no `./*`), `/api` no longer exports the runner,
  `viewer.sources` is read-only (marks are set through `highlight.mark`), and the unused `project`, `selectionItems` and
  `selectionBounds` of the targets port are gone
- Prism viewer: `viewer.extensions.mount(ext)` adds an extension to a running viewer and returns its removal; an
  extension's API area named like one the viewer or another extension has is refused, naming who has it
- Prism viewer: Home's Open Files is a folding block like its actions, holding a drop area that browses on a click (the
  formats it reads on hover)
- Prism viewer: Home's heading is the wordmark as glass (the splash's, without the light or Prism's triangle): it
  assembles as Home opens, leans under the pointer and shatters on a press; nothing is drawn while it rests
- Prism viewer: the empty canvas's splash: with the wordmark (now without Prism's triangle), once the canvas has been
  empty and the viewer idle for 0.7 s (work restarts the wait, so long loads show nothing), light rushes in behind the
  word and leaves the Mol\* star along its arms as bubbles in small spectra of the logo's colors; it answers the pointer
  and touch (pushed aside, dragged and stirred, more light toward it, a burst on a press). The word is glass: it gathers
  smoothly from drifting pieces as the splash starts, its pieces lean a little under the pointer, and a press on it
  shatters it before it snaps back (again while in pieces: further, back from where they are), the snap back drawing the
  light in for a moment like gravity before it is thrown out, unsplit until then. The light goes through air:
  Navier–Stokes on a coarse grid round the star (stable fluids: advection, drag, vorticity confinement, a pressure
  projection), stirred by slow eddies of its own, by the pointer's strokes (the light swirls on behind them) and by a
  shatter's gust, each speck moving at its own velocity plus the air's (a Runge–Kutta step). Beside it, an experiment in
  Storybook (Components / Splash: Quantum, Side by side): the light as a wave, ψ on a 128 × 128 grid moved on by the
  Schrödinger equation (split-step Fourier), emitters radiating the packets (the flow as strong as the glass is
  assembled), the wordmark's glass a potential well its pieces carry (a shatter scatters it; the star splits packets by
  wavenumber like a prism), the specks born where |ψ|² says and riding the probability current (Bohmian trajectories) as
  they fade; the pointer a barrier, a press a kick; on a dark canvas the wave is added light; the wave in the logo's
  colors by its wavenumber (brand-1 the longest waves to brand-4 the shortest). The light comes in as a broad front
  going straight right (no one spot holds it), and starts already flowing, faint (packets along the beam going right),
  in a round vignette. Faster: the FFT with no permutations (decimation in frequency forward, in time back, the kinetic
  table bit-reversed), Strang steps merged, separable emitters, a fast atan2, work only inside the vignette (about 2.2×
  less a frame); the glass's moving pieces drawn from a sprite atlas instead of a clip and the word's strokes each; no
  allocations in the frame paths. The empty canvas picks one at random each time it shows (the wave two times in three;
  `EmptyCanvas`'s `field`; the viewer's `ui.splash`: `random`, `quantum`, `light` or `off`, the app's `?splash=`). While
  the empty canvas's mark is up the canvas is quiet (transparent, its axes too); what loads then comes in gradually (0.9
  s, over the canvas's own color); a load the mark never showed for appears at once. The sequence band fades in and out
  (opacity only, 0.36 s) instead of growing to its height and collapsing (the canvas-drawn band was laid out again every
  frame); `Show`'s `reveal="fade"`, `--msp-motion-fade`. A right panel replacing another swaps as a scoped view
  transition; a lazy panel's first open (State tree, Assistant) swaps like any other (its stub's tab stays until the
  real one is in, instead of closing and opening again). Choosing a UI theme sets the canvas background to the theme's
  (`--msp-canvas-bg`) when the canvas has a plain color (no image, gradient or skybox style); at the start too, unless
  the options give one. The preferences (the theme among them) are kept in localStorage, else sessionStorage, else
  nowhere; every storage access is guarded. Export has its own section in the rail (a hairline above it); the picking
  levels' rows stay 28px with a mouse (44px only on touch). The embedded layout (`ui.layout.compactRail`, the app's
  `?hide-rail`): the rail hidden behind a button in the canvas's corner; Storybook: Viewer / With extensions, Embedded;
  the Seq* introduction uses it. The canvas toolbar fades in and out with the canvas's content (opacity, like the
  sequence band); hiding a compact rail closes its panel first (no flash when it comes back); History in the house
  style. `options.ui.splash: 'off'` turns it off; a still frame under reduced motion; Storybook: Components / Splash,
  with a control for each option
- Prism viewer: a load started from Home (its blocks, Open files, the database and URL forms, a drop while it is open)
  shows the Structure tab once something came in
- Seq\*: `createSeqViewer` takes its spec (no default), the default composition is its own module (`default-spec.ts`)
  and its behaviors and entries are pure, so a host's partial spec bundles only what it lists (an empty spec ~720 kB
  minified, the default ~1,040 kB); the Prism viewer's sequence is composed in one place, `PrismSequenceSpec()`
- Prism viewer API: one operation per state: `trajectory.frame` shows a structure's model (`structure.model` is its
  alias until the release); an empty `representation.remove` or `superpose.remove` is the only one there is, else
  ambiguous, and `all: true` clears
- Prism viewer API: one set of words: `list` for collections (`structure.get`, `volume.get`, `particles.get` and
  `measure.get` are now `.list`, the old ids aliases until the release; the `structure.list` state is
  `structure.loaded`), `add`/`remove`/`intersect` for `component.modify` (was `union`/`subtract`), `preset` for
  `load.entry`'s representation preset (was `representation`), a structure as `{ ref, label }` in `focus.get`,
  `focus.history`, `contents.get` and `component.add`; with no `structure`, an operation uses the focused one, else the
  one holding the selection, and `structure.remove` needs a name when several are loaded
- Prism viewer API: cells of the state tree are addressed by a node selector (`node`), CSS-like text compiled to Mol*'s
  `StateSelection` (`component[type=ligand] > representation`, `#ref`, a ref alone): `state.*`,
  `representation.update`/`visibility`/`remove`, `volume.*` and `particles.*` act on every match as one step (the
  raw-ref `ref` and `representation` fields are gone; a representation type stays a scope); new `state.select` and
  `state.apply` (a transformer under each match, its new refs returned); `node()` writes selectors in scripts
- Prism viewer API: a value left out is kept and a value given is set, with no sentinel values; switches take the state
  wanted (commands flip them): `canvas.lighting`'s occlusion, shadow, outline, fog and illumination are
  `{ name: 'on', params } | { name: 'off' }` (`canvas.illumination` is now its illumination switch), `canvas.xr`,
  `view.expand`, `view.fullscreen` and `volume.streaming` take `on: boolean`; `export.image`'s and `video.export`'s
  `size` is `{ name: 'hd' } | { name: 'custom', params: { width, height } }`; `camera.set` sets the view,
  `canvas.camera` the projection; `canvas.get` reports the lighting, camera, clipping and XR
- Prism viewer API: `animation.play` takes an animation by its name without the namespace (`animate-camera-spin`), its
  title (case-insensitive) or its full id, fails `ambiguous` for a name two animations share and `not-found` (listing
  them) for an unknown one, and with no name plays the one playing, else asks which (`invalid`) rather than the first;
  `animation.get` gives these names (and the full `id`); the Play command is one command whose form lists the animations
  registered when it is shown (the per-animation preset commands are gone); new `sequence.rail` sets the tracks' rail
  (navigator, height; a step), which the `sequence.rail.*` commands now run, and `sequence.get` reports it
- Prism viewer API: one loading vocabulary: `format` is `'auto'` by default (`load.url`, `load.data`, `load.volume`,
  `load.mvs`, `load.trajectory`), told by Mol*'s format registry from the extension (else the data), and `binary` is
  optional (told from the extension or the format); `load.data`, `load.volume`, `load.mvs`, `load.svs` and
  `session.open` take exactly one of `url`, `data` (text, base64 or bytes; `load.data`'s `bytes` is gone) and `file`;
  `session.open` also takes a state server's `id` and a `type` (auto: a zip is a session), and `load.snapshot` is the
  same operation
- Prism viewer: review fixes (2026-10-16, the branch review): `api.group` is all or nothing for every call in it; a
  Structure panel row removes only its own representation; targets, hover and selection work on coarse (PDB-IHM) models;
  the sequence band draws a loaded SeqViewSpec and keeps its band when the followed structure changes; a file that does
  not parse fails instead of reporting success, a session with one unreachable part opens the rest; History shows a
  redone step as done; less work at pointer rate (hover, sliders, the transfer function, snapshots, the Image preview);
  the eager built-in features' operations (`load.entry`, `canvas.lighting`, `session.save`, …) are typed as present; the
  assistant downloads only the model SDK it uses; PDB export keeps modified residues in their chain
- Prism viewer: Home, the Structure panel, the Selection tab and Help are built-in features (`home`, `structure-panel`,
  `selection-panel`, `help`; their tabs keep their ids), so a spec without them bundles none of their views
- Prism viewer: the scene's API is built-in features (`structure`, `volume`, `trajectory`, `annotations`), listed in the
  default spec; a spec without one leaves its areas and commands out of the viewer and the bundle, and reading a missing
  area says which feature to add
- Prism viewer: the query languages target text is read in come from the spec (`languages`; the default spec lists
  PyMOL, VMD, Jmol and MolScript), so a smaller viewer bundles no transpilers
- Prism viewer: marks, a view's residue ranges shown in the 3D view and the sequence until cleared (`highlight.mark`,
  `highlight.unmark`; `highlight.set` with `sequence` lights the sequence too)
- Prism viewer: Help's basics (load, select, focus, the command menu, undo, the assistant) above the bindings
- Prism viewer: an Illumination toggle and, where WebXR is, an XR button on the canvas toolbar (`canvas.xr`).
- Prism viewer: an expanded viewport (`view.expand`: the canvas fills the viewer, Escape brings the panels back) and
  `view.fullscreen`.
- Prism viewer: `window.molstar.lib` in the app (and so in the `molstar` distribution's `build/prism-viewer/`): v5's
  `molstar.lib` without the classic UI, plus `lib.prism` (`createPrismViewer`, `DefaultPrismViewerSpec`,
  `FullPrismViewerSpec`, `node`, `MS`)
- Prism viewer: the Export card's Camera section: the camera in numbers (as v5's), Copy MVS camera; `export.viewport`
- Prism viewer: the Selection tab's summary (what is selected, the structures and chains it spans, Focus, Clear; with
  nothing selected, a way into Select mode); the style layers applied (color, transparency, glow, material, clipping)
  listed under its style form, each edited or removed, Clear all (`representation.layers`, `.layer-update`,
  `.layer-remove`, `representation.reset` with `all`); the superpositions made listed under Superpose and removed
  (`superpose.list`, `superpose.remove`) and `load.viewport` (a viewport state as JSON); `camera.set` and a richer
  `camera.get`. Session & state is v5's State panel: State, Session, Open, and the Save options (`session.save`'s
  `options`). The Image section chooses Image or Animation (the `export.picture` slot, whose views get the shared frame,
  `ExportFrame`); `video.export` takes `axes`.
- Prism viewer: the plugin's toasts show in the status column.
- Prism viewer: the focus history (`focus.history`, `focus.again`, the "Focus again…" command), as v5's focus picker
  lists recent foci.
- Prism viewer: volume streaming's view, detail and channels (`volume.streaming-view`), with a "Volume streaming" row
  and inspector under the streamed structure.
- Prism viewer: a map's levels in the inspector of a volume, an isosurface and streamed density: the map's histogram
  with a marker per iso value to drag (one undo step a drag), the value in σ or map units, each level's color and
  opacity; an isosurface's tree row shows its color.
- Prism viewer: a direct volume's transfer function in the inspector: the opacity curve over the map's histogram (a log
  scale, Mol*'s presets), the area under it in the colors painted there, the volume-value theme's stops on the same axis
  (a stop on a point moves with it), Colors at peaks (a popover: the theme's colors or a palette; a preset moves the
  stops onto its peaks with their colors); previews while dragging are throttled, the final value applies on release,
  one undo step a drag.
- Prism viewer: a structure representation's quick settings above its whole form: its type and color theme as choices,
  one color from a row of swatches, the opacity and the size as sliders (a drag is one undo step); a map level's color
  from the same row.
- Prism viewer: a story's audio (Mol*'s markdown audio commands, MVS narration) has its own bar at the top of the
  playback card (play or pause, the position to click or drag, the time, stop, close) and a canvas button while there is
  a player, so it can be paused or stopped with the card shut; the playback card opens for a story loaded before its
  views came.
- Prism viewer: the quick type and color of a representation drawn as what they are (the types' icons in a segmented
  strip, the themes as tiles in their colors on the structure); pLDDT offered only for a prediction (ModelCIF quality
  metrics, or no experimental method); a component's inspector lists its representations to open; Contents shows its
  counts as figures and its entities as rows (kind, chains, atoms), each opened in the tree; Focus is one row for the
  scene; the Structure panel's header has no Component button (the command menu adds components).
- Prism UI: a select's options are grouped by Mol*'s categories, each group a fold (the current option's open); a form's
  select opens as a popover under its field; `CommandManager.open` takes values to open a command with.
- Prism viewer: a component's Add representation and Set coloring open their commands in the command menu for it, the
  component (or a preset component's keyword) filled in, and one Esc closes a command opened that way; Modify by
  selection acts at once (include, subtract or intersect the selection).
- Prism viewer: the playback card is the Story (a book icon), one frame, its text wrapping inside it; a description's
  first heading that only repeats the snapshot's title is not shown; ← and → step through a story (quiet where there is
  none); the Snapshots tab's Play is an icon beside its options.
- Download Structure from ModelArchive reads `modelarchive.org` (`www.modelarchive.org` does not answer a browser's
  cross-origin request, so the download failed).
- Prism viewer: a snapshot's ⋯ menu is a popover with Remove (its order is the drag handle's: no Move up and down);
  prism-ui's `MenuPopover` serves it and the tree's row menus.
- `StructureComponent` (a static one: the presets' Polymer, Water, …) takes a new label on update (it returned Unchanged
  before reaching the label, so a rename did nothing).
- Prism viewer: a component is renamed from its inspector's title (a pencil); Make a component can make one with no
  representation (v5's "< Create Later >") and opens the Structure panel at it; Export is last in the rail's scene
  group; Home's Zenodo starts folded, among the actions.
- Prism UI: switching from one card to another (Animation, Canvas) swaps its title and body in the one card, as a tab
  swaps a panel's content, instead of closing and opening a card again.
- Prism assistant: your own OpenAI-compatible endpoint runs turns: its optional context window and max output fields
  (32,768 and 4,096 tokens when left empty) become the provider's `modelLimits`
- Prism assistant: a call waiting to reach an address offers Allow once and Allow <host> (Remember keeps it in this
  browser, else for the session); the next calls to that host go without asking; the Storage view lists the allowed
  domains, each with Forget.
- Prism viewer: a tree row's ⋯ menu opens as a popover from its button (as the selects and presets do), not inline under
  the row; Escape or a press outside closes it. Home's extension sections (Zenodo) line up with its actions.
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
- Prism viewer: saving a session can include its record (the form's Record, off by default: prompts and replies are in
  it), as `record.jsonl` in a `.molx` or a `record` field in a `.molj`; opening such a file shows the record in History,
  read-only.
- `PluginStateSnapshotManager.serialize` takes `extra` (files in a `.molx` beside `state.json`, top-level fields in a
  `.molj`), and `open` returns what it did not use, so a host can keep its own data in a session file.
- Prism viewer: Home lists the records kept in earlier sessions; one opens in History, read-only, beside what is loaded
  ("This session" goes back); each can be forgotten there, and this session's kept copy from History. At most
  `ui.keptRecordsMax` (20) are kept: keeping one more names the oldest and asks before forgetting it.
- Prism viewer: the Story card stays open on a click outside it (Esc or its ✕ close it) and its description fills the
  card; a first heading is left out when the title already says it (a part counter like "[1/2]" aside), whatever order
  the title and description change in. Overlays can be `persistent` and `fill`.
- The camera: a transition under way when the animation loop's clock starts again (after an image is drawn, e.g. a
  snapshot's thumbnail) keeps going instead of flinging the camera far out (a negative or huge radius, an empty view); a
  state snapshot keeps where an animating camera is going (`Camera.getSettledSnapshot`) and the view a pending camera
  reset frames (after a load), not the camera before it.
- Prism viewer: an inspector's header has the item's eye at its end (Show or Hide; Turn on or off for the unit cell, the
  streaming and an extension's row), in place of a Turn on row while it is off.
- Prism viewer: the canvas's playback card has a row each for the snapshots (the current one, its position and the list,
  previous, play, next), the current snapshot's transition (a scrubber, play from the start or stop: `snapshot.seek`,
  `snapshot.transition`) and the trajectory's models; extensions add rows in the `playback` slot.
- Prism viewer: the sequence's No colors is plain letters (a neutral strip; plain text in the text style) instead of
  Seq*'s default residue colors.
- Prism UI: a form's presets open as a popover under their button (`PresetsList`), and a list inline in a panel lets the
  wheel scroll the panel at its ends (`ScrollArea`'s `chain`).
- `ColorTheme.Palette` (the volume-value theme): a color list with offsets interpolates between its stops (as
  `ColorScale` does) instead of switching halfway between them, and offsets with many decimals no longer make a huge
  palette.
- Prism viewer: load a SeqViewSpec into the sequence (`load.svs`, `.svsj` files, `svs-url`): it shows over the loaded
  structure it names (loaded if missing) until another structure is focused or loaded, and Undo and sessions keep it.
- Prism viewer: `load.mvs` takes v5's options: `append` (keep the Snapshots list), `keepCameraOrientation`, `snapshot`,
  `extensions`.
- Prism viewer: snapshot descriptions run Mol*'s markdown commands (focus, highlight, camera, snapshots, audio) and show
  MVSX images, color swatches and tables (`snapshot.link`; prism-ui's `Markdown` takes a host).
- Prism viewer: the playback card opens by itself when a loaded story has snapshots, and shows a single snapshot's
  description.
- Prism viewer: the image export's crop box can be dragged and resized (and moved with the arrow keys).
- Prism viewer: Help lists the canvas keys (view, moving, snapshots).
- Prism viewer app: the `map-provider`, `show-toggle-fullscreen`, `debug-mode` and `timing-mode` URL parameters are
  read.
- Prism viewer: the structure inspector sets the model of a trajectory and the structure type (the model, an assembly,
  symmetry mates or unit cells); the structure row's menu has Preset, Assembly, Crystal symmetry and Stream density.
- Prism viewer: extensions' rows in the scene tree (`ViewerExtensionUI.rows`, drawn with `SceneRow`); assembly symmetry
  and ANVIL's membrane each have a row under the structures they apply to.
- Prism viewer: the `<prism-viewer>` element (`@molstar/prism-viewer/element`; `spec`, `pdb`, `theme`, `hide-controls`),
  defined by the single file.
- The Prism viewer loads its assistant on first use (the Assistant tab, the Ask row, or `?assistant`), not at start; the
  command menu renders its options as they scroll in; History reads only what is new
- `@molstar/ai-tools`: a dynamic toolkit (`toolkit.add`, lazy skills disclosed during a turn); `ai-session` reads it as
  it is
- Prism viewer: Seq* is its own chunk (`sequence.runtime` in the spec), loaded when a structure is first shown: the
  app's first load is 1.25 MB gzipped (from 1.68)
- Prism viewer: the assistant's Inspector (Trace, Tools, Limits) over the real agent, opened from the Assistant tab's
  Inspect rows, the turn card's trace and the Assistant trace/tools/limits commands; extensions' views may add `pages`
- Prism viewer: the params form lays params out as v5 does (a choice's own params behind its settings button, an
  Advanced fold, a fold per category); the scene tree's focus group (always under a structure, holding the focus target
  and surroundings); Preset and Styles open the command menu in their groups; cards stay open while the canvas is
  turned; a lazy extension that fails to load says so and offers Retry
- The snapshot manager's `getStateSnapshot` and `serialize` take `keepList` (the current scene without changing the
  list)
- State actions and transformers take `display.triggerLabel`, the label of the button that runs them ("Download",
  "Open"); a UI falls back to a generic one without it
- `PluginState.registerCustomState(key, { get, set })` lets a component keep its own state in plugin snapshots
  (`snapshot.custom[key]`): saved by `getSnapshot`, restored by `setSnapshot` after the data tree (`undefined` when the
  snapshot has no value), carried by snapshot manager entries and `.molx` / `.molj` sessions
- Every plugin path that changes the canvas props emits `events.canvas3d.settingsUpdated` (initial props, `clear(true)`,
  snapshot, transition-frame and snapshot-animation restore, canvas context props); new `PluginContext.setCanvas3DProps`
  and `setCanvas3DContextProps` set the props and emit. A direct `canvas3d.setProps` still emits nothing

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
  outside the int32 and uint32 ranges as a float 64 (it threw), and the keys of an object whose value is `undefined` are
  left out (they were written as nil). BinaryCIF files no longer carry a nil `mask` key
- MessagePack decoding of untrusted input: a cut buffer throws (it read as data), a length header that claims more items
  than the bytes left throws, arrays grow as items are read rather than at the claimed length, and a `__proto__` key is
  an own key, as `JSON.parse` makes it. On a 340 KB `.bcif` the checks cost about 3% of a decode
- Fix UTF-8 encoding of characters outside the Basic Multilingual Plane (`utf8Write`, `utf8ByteCount` in
  `@molstar/core/util/utf8`): a surrogate pair was written as two 3-byte sequences instead of one 4-byte sequence; ASCII
  and BMP output is unchanged
- Zip hardening (`@molstar/core/util/zip/zip`, a security fix): `unzip` looks for the end record only in the last 64 KB
  (a truncated archive looped for ever), reads the central directory first and takes sizes, name and CRC-32 from it
  (archives written with a data descriptor, as macOS Finder writes them, have zeros in the local header), checks that no
  local entry is listed twice, and that every entry has its declared size and CRC-32; `inflateRaw` and `inflate` with a
  buffer stop at its end and return what they inflated. New optional `UnzipLimits` (`maxEntries`, `maxBytes`) as the
  fourth argument of `unzip`. Invalid or damaged archives that were read before now throw

- Prism viewer: a representation's row menu (⋯) offers Type and Color (the current ones checked) and Edit…; a load from
  an address that does not answer reports the failed fetch, not an internal error; the camera's field of view and
  clipping fields have readable labels.

## Seq\* input streams

- The Seq\* viewer's input for consumers is `ctx.interaction`: rxjs streams (`hover`, `press`, `click`, `drag`, `menu`,
  `wheel`, `key`, `leave`, `activate`) of fresh payloads in Mol\*'s button and modifier terms, scoped by `on`, built
  only while observed; the behaviors' base events are `ctx.input.events`. The track menu is `interaction.menu`, the
  browser's context-menu request (the `contextMenu` behavior, `events.menu` and `events.activate` are removed).
  Behaviors list their bindings as Mol\*'s `Binding` (`ctx.behaviors.bindings()`)

## Seq\* line styles

- A `ranges` layer's `line_style` channel (format 1.0): `uniform` with a `line_style` or `field_line_style` over a
  column of the rows draws each item `solid`, `dashed`, `dotted`, `wavy`, `double`, `box` or `highlight`, in tracks,
  text and SVG exports; `bookmark.add` takes an optional `line_style`. `PathBatch` strokes may be dashed (`dash`, `gap`,
  `phase`), and `Paint2D` includes `setLineDash` and `lineDashOffset`.

## Seq\* selection

- [Breaking] The `select` behavior's `shiftPositions` defaults to `toggle`: Shift with a click or a drag on positions
  adds or removes them (it extended a range from the anchor); `extend` stays a param

## Asset manager byte cache

- `AssetManager.resolve` reads URL assets through a URL-keyed, ref-counted byte cache: a URL (and its headers) is
  downloaded once whoever asks first, concurrent requests share the download, the bytes are held while any resolved
  asset is held and dropped with the last `dispose`, and a failed download is not cached. Requests with a body and
  `store: false` bypass it

## Prism UI cleanups

- Duplicate command ids throw when registered, naming the id (they replaced the earlier command silently).
- `ReactiveModel` has a `mountSignal`: an `AbortSignal` that aborts when the mount ends.

## Prism XR

- The Canvas card has an XR section (the headset view's resolution scale, scene size, post-processing and minimum target
  distance).
- `options.ui.xr` (`auto`, `always`, `never`; classic `ShowXR`) decides the toolbar's XR button: with `always` it shows
  where WebXR is not too, and clicking says so.

## Prism session record

- History shows what a run was pointed at and what it changed ("1CBS · A:10-20 · 3 nodes, view"), the agent turn's model
  (subtly), and the reply as the model says it, not only at the end.
- An Undo turn is one record entry (naming every step), and an undo or redo records who made it (the person, the
  assistant, a script): History says "Undone by the assistant".
- The kept records' index folds itself (one entry per record) once it has grown, so keeping and touching records for
  long does not grow the browser's storage without bound; no kept record is ever dropped by it.

## Prism layout

- The viewer keeps the layout the person leaves (docks, sizes, hidden tabs, which panel each dock shows) in its
  preferences and restores it at start (`persist.layout`, on by default; not for an embedded layout).
- Reset layout (`view.reset-layout`, and a rail tab's menu) puts the layout back to the spec's.
- A rail tab's menu (right click, Shift+F10) hides that tab or shows a hidden one; `view.show-tabs` shows them all, and
  opening a hidden tab by name brings it back.
- `ui.layout.mode: 'embedded'` is the layout for a viewer inside a page: the rail behind its corner button, the panel
  floating and never pinned, compact where narrow, no panel open at start.

## PAE from AlphaFold DB

- A structure that is an AlphaFold DB entry (an `AF-<UniProt>-F<n>` id in its entry id, name or address) can fetch its
  predicted aligned error from AlphaFold DB: "Fetch PAE from AlphaFold DB" in the structure row's menu, the
  `model-archive.fetch-pae` command and operation, and a button in the PAE card. The PAE is a download and a parse in
  the state tree, so it is undone with Undo and saved in a session.
- The Structure panel's PAE row says the right size for a partial ModelCIF metric (a domain's PAE of a bigger model).

## Prism commands: Select, keys, menus by command id

- The command menu has a Select group: every selection query of the plugin's registry by category (By query), and the
  loaded structures' entities, chains, ligands and other residues, and elements (In the loaded structures), with one
  mode (Replace, Add, Remove, Intersect) set in the group's header. The highlighted row lights its atoms in 3D and says
  how much it would select ("1 residue · 20 atoms").
- Keys can be remapped in Help (Change keys: record, clear, reset, give a key to a command without one, reset all); the
  person's keys are kept with the viewer's preferences. `ui.keybindings` is the app's layer under them. Help, the
  command menu and `?` show the keys in effect, and conflicts per scope.
- Commands can be scoped to a region (`scope: 'canvas'`, `'panel:<id>'`): their key works only there, and before a
  global command on the same key.
- prism-ui: `commandItem` / `commandItems` (menu entries by command id) and `CommandButton` (a button by command id);
  context menu entries show a key. The canvas's camera Reset menu is made of commands.

## Prism: virtualized scene and state trees

- The Structure panel's scene tree and the State tree tab render only the rows on screen: a structure of thousands of
  components or a state of 5,000 cells costs what a screenful does, and opening or closing a row only recomputes the
  list. The keyboard (arrows, Home, End, Right and Left to open, close, enter and leave, Enter) moves across rows that
  are not rendered, scrolling them into view; screen readers get each row's level and place (`aria-level`,
  `aria-posinset`, `aria-setsize`).
- The State tree's ⋯ menu is a popover (as the scene tree's), not a list under the row. An extension's row in the scene
  tree never opens its controls under it: they are in the inspector.
- prism-ui: `VirtualList` measures its rows' height, takes `edge`, `contentRef` and a handle (`scrollToIndex`);
  `flattenTree` gives each node its place among its siblings; `keepRows` keeps a flattened list's unchanged rows;
  `TreeNavigation` (`dom/tree-keys`) is the tree keyboard over a flattened, virtualized tree (it replaces `treeKeys` and
  `treeFocus`).

## Prism: several viewers on a page, and a camera to follow

- Viewers on one page no longer share what each owns: the query languages are a set per viewer (the page-wide ones only
  seed it), the kept records are in an IndexedDB named after the viewer's `persist.key` (`prism-viewer-records` for the
  default key), the assistant's model choice and allowed domains are kept per viewer, fullscreen follows only its own
  element, and the assistant's Esc and Inspector keys go to the viewer that has the focus. The assistant's keys stay the
  person's, shared by every viewer in the browser.
- `viewer.camera` (the view now, `changed`, `set(view, { silent })`) and `syncCameras(viewers)` keep viewers' cameras
  together with no echo.
- Stories: Viewer › Several viewers (three independent viewers; docking poses with a synced camera).

## Seq\* and Prism: coordinate systems with no viewer, the sequence's caches, hiding tracks

- Seq\*: `@molstar/seq-core/translate` has the `SystemRegistry` (a set's own systems by key, `seqstar:<ref>`, beside a
  host's, changes through `onChange`) and `Coordinates` with their conversions; `registry.map(target, to)` maps between
  a structure, a sequence set, an alignment's columns and a transcript with no viewer. The viewer's `ctx.systems` is one
  (`ViewerSystems`, plus its reactive `version`); the Mol\* adapter's `StructureLink` attaches to any registry.
- Prism: a SeqViewSpec loaded into the sequence reads a file the viewer already holds (a downloaded or opened file, by
  its address) and shares a structure it already parsed, instead of refusing it; files it does not hold are still
  refused, saying to load them first.
- Prism: a sequence track can be hidden from its right-click menu (Hide track) and brought back from the Sequence menu's
  Visible tracks submenu (a checkbox per track, Show all) or the faint "N hidden" under the tracks. Each change is one
  Undo step; hidden tracks are saved with snapshots and sessions and stay hidden across another structure, new colors
  and a reload. Commands `sequence.track.hide`, `sequence.track.show`, `sequence.track.show-all`; operations
  `sequence.track.hide` and `sequence.track.show`; `sequence.get` lists the tracks.
- prism-ui: a canvas menu row with `items` is a submenu, opened beside the menu (hover, click, →).

## Prism: automated accessibility checks and what they found

- prism-ui: `@molstar/prism-ui/testing/axe` (`expectNoAxeViolations`) runs axe-core in the unit tests; the shell, the
  components, every tab and card of the viewer, the command menu and the inspectors are checked
  ([accessibility.md](../../packages/prism/ui/_spec/accessibility.md), which also has the audit by view and the gaps).
- Fixed: an open select's options list has a name; the compact tab bar's menu buttons are beside the tab list, not in
  it; the Snapshots list is a list of items (a button applies a snapshot, its buttons are beside it) and its buttons are
  24px, 44px in compact mode; a volume's iso slider gives ARIA a plain number.
- CI: the Prism packages' `styles:check`, `icons:check` and `typecheck:storybook` run in a job of their own.

## Seq\*: alignment figures (ESPript-style), frames across tracks

- An alignment's `rows` draw in the text view: a block of text rows on every line, each row named on the left and
  numbered on the right in its own numbering, the group numbers above each line the reference row's. A streamed
  alignment's rows keep their fallback (the reference row's letters) in text.
- The style key `text_cells` (render prop `wrappedCells`): letters and rows drawn on their theme's cells in text, each
  letter in the ink that reads on its cell.
- The `espript` color theme and scheme: identical columns black, residues whose group holds 70% of the column grey (by
  groups, gaps not counted); its dark variant on a dark skin.
- A new layer type, `frame`: named regions outlined across the tracks `span` names (`from`, `to`), titled with each
  item's label, open on the side where a line break or the view cuts the item; over a rows track, the rows the item
  targets. In text and tracks, the SVG export alike; `svs.layer.span` for a span outside the section.
- `secondary`: helices are cylinders (rounded where the helix ends, square where it is cut); `labels: true` draws each
  item's label (α1, β2) over its shape.
- A range's `box` or `highlight` on an alignment row draws on that row in the text.
- The Annotations part: a place's `kind` (`ranges`, `points`, `frame`); a place of markers takes a glyph (the row's
  `glyph` column, `bookmark.add`'s `glyph`) instead of a line style.
- Stories: Text view › Alignment figure (the core's globins as an ESPript-like figure, with the panel's Frame, Box on a
  row and Markers places); Publication tracks gains a Frame and Markers.

## Prism: a compact canvas toolbar, bounded menus

- prism-ui: the canvas toolbar can be one button whose menu lists the tools (toggles as checkable rows, the mode
  choices, each tool's own menu as a submenu). It is the host's choice, set at creation by `ui.layout.toolbar`
  (`ShellOptions.toolbar`): `auto` (the default: the button in the compact mode, the bar otherwise, following the mode
  as the shell resizes), `expanded` or `collapsed`. The person has no control for it, and it is not part of the kept
  layout. `ShellModel.toolbarCollapsed` (read-only); `ShellView`'s `canvasLeadMenu` gives the menu the rows of the app's
  own lead buttons (the viewer's Interact, Select and picking level).
- prism-ui: the canvas toolbar's menus and submenus are as tall as the room between the canvas's top and their button
  and scroll inside it, instead of running off the top of a short viewer.
