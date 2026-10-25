# Prism: to do

What is missing, and what was considered and deferred. One entry each, with a pointer (a code path or a `_spec`
section); an entry leaves when it is done (its spec section says so then). The specs:
[viewer/\_spec](viewer/_spec/README.md) and [ui/\_spec](ui/_spec/README.md). The AI package (`@molstar/ai`) has its own:
[../ai/TODO.md](../ai/TODO.md); Seq\*'s is [../sequences/todo.md](../sequences/todo.md).

## Decisions open

- **One schema vocabulary; reconsider zod** (David). The workspace describes and checks data three ways: Mol\*'s
  `ParamDefinition` (operations, options, Settings: `viewer/src/options.ts`, the API areas), io-ts in MolViewSpec
  (`mvs/builder/src/tree/generic/field-schema.ts`, the builder and the runtime) and zod in `@molstar/ai` (the agent's
  tools: field shapes from each operation's `ParamDefinition`, `viewer/src/agent/params.ts`, turned into zod, then into
  the JSON Schema the model sees and the validator). The question is the workspace's, not the agent's: one vocabulary
  that gives types, validation with good messages and JSON Schema (for the model, for MolViewSpec's published schema,
  for SeqViewSpec), or each kept where it is. zod is never in the page's startup bundle (the tool declarations are plain
  data, `viewer/src/agent/skills.ts`; the assistant is a chunk, `assistant/src/spec.ts`), but it is ~192 kB minified in
  the assistant's chunk (English only, as `scripts/esbuild/app.mjs` builds it; ~453 kB with every locale, 95 kB
  gzipped), and the AI SDK's providers import it, so it leaves the page only with providers of our own. Candidates: JSON
  Schema straight from `ParamDefinition` with a small validator over that closed vocabulary; zod everywhere (v4 emits
  JSON Schema; MolViewSpec off io-ts); a smaller validator. The harness side, slice 12b: [../ai/TODO.md](../ai/TODO.md).
  Measured 2026-10-17 (esbuild, minified / gzip, English-only zod): a provider chunk on the AI SDK is 570–750 / 155–195
  kB (zod 193 / 54, `ai` 184, the unused gateway 111); the official SDKs (`@anthropic-ai/sdk`, `openai`,
  `@google/genai`) are 315–400 / 72–84 kB with no zod on the streaming path. A rewrite of `providers/` on them is ~900
  lines and 2–3 weeks with the fixtures; `tools/input.ts` then needs a ~120-line validator over the closed field
  vocabulary (or `zod/mini`, 18 / 6 kB).
- **Assistant keys are shared by every viewer in the browser** (README decision 29): Forget in one viewer is seen by
  another only after a reload; make the key store notify the other viewers if that matters.
- **Size**: MolViewSpec in the startup bundle (~197 kB, a lazy MVS feature with stub format providers); the Seq\* chunk
  (one default for now, everything but streaming; its floor of ~721 kB is aligners, SeqQL, streams and palettes:
  [../sequences/todo.md](../sequences/todo.md)); eager extensions in the app (~190 kB, lazy stubs owning their formats);
  panels shown on demand as static imports (command menu, History, Settings with the params form, volume inspectors);
  `sideEffects` declared almost nowhere (every module that registers on import must be listed).
- **Build budgets** (2026-10-16, app.md §1): the `index.js` budget no longer bounds the start (the start code is in
  shared chunks since `lib`'s lazy import): a budget over the start closure (65 files, 4.4 MB, 1.29 MB gzipped); lower
  the 8 MB `prism-viewer.js` budget (6.2 MB today).
- **`molstar.lib`** (`apps/prism-viewer/src/lib.ts`): also from the single file `prism-viewer.js` (+32 kB gzipped on
  every page that imports it)? `StateTransforms` for full v5 parity
  ([architecture.md](../../.v6/designs/architecture.md) §9.3 keeps it for the classic Viewer)?
- **Choices to keep or change**: built-in areas typed as present and throwing when left out; panels keeping their tab
  ids; History and Settings in the core; the runner re-checking `exists` only for loads.
- **Mol\*'s other apps on Prism**: docking-viewer (its preset buttons as toolbar items and commands), mesoscale-explorer
  (its panels as views; the test of the shell's modularity), mvs-stories (embeddable viewport and snapshot-stepper
  elements), and examples as API checks (`react`, `alpha-orbitals`, `ligand-editor`): port, or leave them on the classic
  UI (`apps/`, `examples/`).

- **Coloring a structure by data: two core color themes** (David, 2026-10-24; not decided, a sketch to think through).
  Aim: any sequence track (an extension's included) offers "Color structure by this", and MVS's coloring shares the same
  pieces. The palettes are already one implementation in core (`@molstar/core/util/color/value-palette.ts`, #2041).
  - `data-values` (graphics): values from a source through a palette. Params: `source` (hidden `PD.Value`,
    `{ name, params }`: a registered value source that recomputes the values, so Undo, snapshots and sessions keep only
    the reference), `palette` (`ValuePaletteParam`), `granularity` (residue or atom), `background`, `tag`. Applicable
    only where a registered source can give values; hidden from the pickers. Sources: computed from the structure
    (hydropathy, B-factor), read from the file (pLDDT, MVS annotations), held by a loaded SVS document. A sequence track
    whose layer is `values` or `ranges` with a `field_palette` would derive its menu entry from the layer (its data
    mapped to residues by `StructureLink`, its palette translated), no extension code needed.
  - `selector-colors` (graphics): few explicit selector → color items, inline. Params: `items` (hidden `PD.Value` of
    `{ selector, color }`: a static part, a schema item or several, painted in order), `background`, `tag`. For "chain A
    red", an SVS highlight, the agent's coloring; MVS's `color` nodes.
  - `tag` (this theme only for now): hidden `PD.Value<{ owner, id?, label? } | undefined>` on the theme's params, never
    read by the theme; provenance for the UI (which track's coloring is lit in the gutter, the legend's title, labels,
    the session record, a coloring whose track is gone).
  - Schema items and their matching to elements, now MVS's (`ComponentExpressionT`, `annotation-prop.ts`,
    `element-painter.ts`), would move to `@molstar/model` (a stacked PR, as #2041), for `selector-colors`, an
    MVS-annotation value source and MVS itself.
  - Open: `items` hidden or editable (`PD.ObjectList`); `source` one value or a name and its params apart; the "what you
    see" alternative (the band's resolved colors per residue) for layers whose color is not a palette over a field.

## Missing: the viewer

- **Scene defaults for maps, particles and shapes** (2026-10-25, [api.md §2 "Look"](viewer/_spec/api.md)): the scene
  defaults (the Structure panel's sliders) are the component manager's options, which only structure representations
  follow; a volume's, particles' or a shape's look is set by node and is its own. Following them too needs the volume
  and particles managers (and the shape loaders) to read the options and `setOptions` to reach their representations.

- **The Bonds field follows the registry live** (David, 2026-10-20): the structure's "Bonds" field (`sources/scene.ts`,
  `sourceForm`; [bond-providers plan §4.1](../../.v6/plans/bond-providers.md)) reads `plugin.bondProviders` when the
  form is built. A provider added or removed while the inspector is open shows after the next pick. Needs a change event
  on the plugin's registries (none has one today); the representation and theme selects would follow from the same
  event.

- **The Sequence node, open points** (2026-10-24, [sequence.md §3, §7.1, §8.1](viewer/_spec/sequence.md)): a provider
  whose track Seq\* refuses is put in again by every build (each refused, the document left as it was): drop its `in`
  after a refusal; hiding the band from the Sequence eye or No sequence first moves the band, then fades it (the top
  dock's hide, not the eye's); version 3's hidden tracks are kept only when the followed structure was saved (none under
  Auto): read them as every structure's instead?
- **More sequence tracks from the structure** (David, 2026-10-24): B-factor (a Seq\* computed provider over atom values
  in `viewer/src/sequence/molstar/`, none today), SIFTS coverage (`molstar:sifts_coverage`) as a provider.

- **Split viewports** (David, 2026-10-17): one Prism instance hosting several plugin contexts side by side (one shell,
  command menu, history and assistant over several canvases): comparing structures or states, linked cameras and
  selections, a 3D view beside a 2D one, a detail view of a pocket next to the whole. Today a viewer is one plugin
  (`viewer.ts`); the API, the sources and the journal assume it. Needs: which viewport an operation and the agent
  target, the shell's focus between viewports, per-viewport sequence bands, and the WebGL budget across them
  (`sources/canvas-presence.ts`).
- **A canvas coming back into view is dark until its context is restored** (~1.2 s after a far jump in the Fifty viewers
  story: Mol* rebuilds its resources). Keep a snapshot of the last frame when a context is released and show it until
  the restore draws (`sources/canvas-presence.ts`).
- A fast path in the runner for pointer-rate view operations; measure first: the hover is per element, not per frame
  (`viewer/src/api/run.ts`).
- `project()` and `selectionItems()` in `viewer/src/targets/evaluate.ts` read atomic-only properties (unused).
- The drop overlay and Home read `sources.formats` once, when the shell mounts (`viewer/src/ui/shell.tsx`,
  `viewer/src/features/home/shell.tsx`).
- Node selectors: `:has(...)`, positions (`:first`, `:nth(n)`); `node` on the structure-scoped operations
  ([nodes.md §5](viewer/_spec/nodes.md)).
- Themes and presets that fetch as they apply carry no addresses (only known services today;
  `viewer/src/agent/policies.ts`); the pLDDT track reads B-factors where the 3D theme reads ModelCIF scores
  (`packages/extensions/model-archive/src/prism-viewer.ts`).
- Mol*'s first framing of a scene (on load) frames it in the whole canvas, partly under the sequence band; resets frame
  below it (`viewer/src/sources/framing.ts`, [README](viewer/_spec/README.md) decision 28). Needs a reset hook in the
  graphics package's camera reset, after the WebGPU work.
- The generated form in a dialog ("open as form") for command params that are not flat; today they run at their defaults
  from the inline sentence (`ui/src/forms/flat-fields.ts`, [commands.md](ui/_spec/commands.md)).
- The app: hashed entry names; publishing from the release workflow (today `pnpm deploy:remote` publishes the app and
  its single file at `prism/` and Storybook at `ui-storybook/`) ([app.md §1–2](viewer/_spec/app.md)).
- Before Prism becomes the default viewer: accessibility and touch beyond the automated checks (a screen reader pass,
  `target-size` and contrast in a real browser, the sub-24px controls, a single-pointer alternative to drag-reordering;
  [accessibility.md](ui/_spec/accessibility.md) "Audit by view"), compact mode, several viewers on a page, performance
  budgets (a one-node change in a scene tree of 1,000 components or a state tree of 5,000 cells within a frame), and
  [parity.md](viewer/_spec/parity.md)'s missing rows; the classic `apps/viewer` stays beside it for a major version.
- Kept data: preferences keep the theme and the layout, not settings; kept records past 90 days are only marked; the
  Storage view does not refresh while open and does not list the assistant's own stores (`prism-assistant.*` in
  localStorage, `prism-assistant-keys` in IndexedDB) ([journal.md §6–7](viewer/_spec/journal.md)).
- Sources: no state of the latest registered presets, themes and actions (only Home's actions follow registrations,
  `viewer/src/sources/home.ts`); none for trajectories from coordinate files or segmentations
  ([state.md §1](viewer/_spec/state.md)).
- Parity gaps: SB-NCBR tunnels (their preset fetches from a service the policy does not know), Volumes & Segmentations
  and the PAE view (to be redesigned), XR checked on a headset, an extension's `reactRenderFn` in snapshot descriptions,
  a host's own bindings in Help ([parity.md §4–6](viewer/_spec/parity.md)).
- Targets: `src/targets/` beside `@molstar/query-language`, so MolViewSpec and stories share it; a schema field for
  alternate locations (query text only today, [targets.md §5](viewer/_spec/targets.md)); SeqQL's `molstar.select` reads
  with Mol\*'s transpilers, not the viewer's target compiler (`viewer/src/sequence/molstar/query.ts`).
- The sequence: a test of every built-in track provider and behavior together over a structure switch, a load and a
  clear; Seq\*'s own document history is on under the viewer's ([sequence.md §7.3](viewer/_spec/sequence.md)); track
  menus offer only coloring and Hide track (Show in 3D, results to `select.set`, hover-card actions: §5); a pLDDT tint
  of the letters in Text (§10); several sequence views (several Sequence nodes lit at once, a section each in the band:
  one `StructureLink` per complex, the band's height), focus-driven views (§3); Seq\*'s `trackLegends` for data-driven
  legends (`viewer/src/sequence/legend.ts`); a third overlap color for marks (`viewer/src/sequence/link.ts`,
  [marks.md](viewer/_spec/marks.md)).
- Embeddable parts: a viewport alone (canvas, hover label, toolbar) for mvs-stories and the docking viewer, and the
  island elements (`<msp-viewport>`, `<msp-sequence>`, `<msp-params-form>`, `<msp-snapshot-stepper>`;
  [rendering.md §3](ui/_spec/rendering.md)).
- Stories, its own app on prism-ui (there are no modes; `ModeChooser` is gone), is not built.
- Turning illumination on waits on its shaders' linking (core graphics, `packages/graphics`).
- The `dom` runtime's edge cases (`ui/src/dom/control-flow.ts`, `render.ts`): `For` shows its fallback while every entry
  waits on a time slice (`order`, not the items, decides); `Switch` with no cases never shows its fallback, and drops
  nested case arrays; a slice ends at a microtask, not a task boundary; a region rendered re-entrantly can insert
  content built in a disposed scope (no guard). `combineBehaviors` projects from combineLatest's raw values, so states
  set in one `batch` emit an intermediate pair (`ui/src/reactive/combine-behaviors.ts`).
- Core follow-ups, each additive (the viewer's sources work around them, `viewer/src/sources/`): a registry-changed
  event from `PluginContext.register`; what changed in `cell.stateUpdated`; `layout.setProps` and every
  `canvas3d.setProps` path emitting; running tasks as a behavior with one id space
  (`plugin/core/src/util/task-manager.ts`); stable string ids beside `PluginCommands`.
- **Clipping** ([clipping.md](viewer/_spec/clipping.md)), still open:
  - tests of the gizmo controller (the active cut following renumbering, the cancel paths: Esc, a change under a drag, a
    failed commit), and a session saved and opened with drawings;
  - the component manager taking global clip objects away shifts the later ones (`replaceOwned`): their keep-whole
    layers are not remapped then;
  - stories for a volume cut and a clip set from code; a mesoscale structure's long Cuts through list measured;
  - an untagged cut from code is tagged when drawn, so code that set it no longer finds it by its values.
- **Gizmo handles behind what is drawn cannot be grabbed** (`@molstar/graphics/util/gizmo`): picking takes the nearest
  object, so a ring or an arrow inside a structure is hidden from the pointer by it; 3D tools draw handles on top. Do it
  as render layers, not an "on top" flag: a representation's layer (its state), the layers drawn in order, each with its
  own depth (cleared, or kept), picked from the top layer down. Handles are one layer above the scene; the same takes
  later effects (an overlay, a highlight layer, a layer with its own post-processing). Canvas3D's helper scenes (the
  corner axes, `helper.handle`) are the model. A first layer can grow out of `HandleHelper`'s scene
  (`canvas3d/helper/handle-helper.ts`: drawn after the scene's depth is detached, picked after the scene, its own loci):
  any number of render objects in it instead of its one built-in handle, `VisualSet` adding gizmos there, and the pick
  pass's depth cleared before it if a handle behind a structure still loses the pick. A renderer change: Alex, after the
  WebGPU backend lands.

## Missing: prism-ui

- Commands: Select's one-run mode by a modifier on Enter (Shift+Enter removes, ⌥↵ intersects) and mode words in tails
  ([commands.md §5](ui/_spec/commands.md)); trees marked as a key scope (`data-msp-scope="tree"`) and commands using
  scopes (§6); more menus and buttons by command id: the scene tree's row menus, the selection bar, `CommandButton` in
  views (§7).
- No lint rule against `.value` inside JSX ([rendering.md §2.6](ui/_spec/rendering.md)); no `usePrism()` for React views
  (§5).
- States handed to library components as props (`SelectButton`'s `label`, `Button`'s `pressed`, `Select`'s and
  `TextInput`'s `value`, `ColorPicker`'s): typed `MaybeReactive`, so each is still a binding of its own
  ([rendering.md §2](ui/_spec/rendering.md)); a `Child` label or a watch inside would end the casts. A scene-tree row's
  layout (~0.1 ms a row at the frame) is now its largest cost: its CSS wants a look.
- A visual check of the virtualized scene and state trees (light, dark, compact) in Storybook; no scene-tree case in the
  performance story.
- The subscription guard covers prism-ui only (`ui/src/_test/subscriptions.test.ts`); the other Prism packages want it.
- The `--msp-ink` tokens are unused outside the canvas layer.

## Missing: the agent and the assistant

- **Tool and skill definitions need a redesign** (David, 2026-10-17; do separately). A skill's `guidance` restates its
  tools' ids and descriptions (`packages/extensions/kinemage/src/prism-viewer.ts`: `KinemageSkill` names `get_kinemage`,
  `show_kinemage_part`… in prose and again in `tools`), so the model reads them twice. `OperationToolSpec`
  (`viewer/src/agent/skills.ts`) is a second vocabulary beside the operation's own definition: a model-facing id, a
  description that repeats the operation's title, `fields`, `fixed` and `formats` per tool. Wanted: tools derived from
  operations (their title, params and descriptions; one name), a skill whose guidance says only what the tools cannot
  (when and why, how they combine), and every extension's skills moved over ([agent.md §2–3](viewer/_spec/agent.md)).
- Skills: SeqQL selections and the sequence skill, `viewer-code` (TypeScript against the API, api.md §7), opening panels
  as tools ([agent.md §3](viewer/_spec/agent.md)).
- `ask` parking the turn with the question in the turn card; today it ends the turn
  ([agent.md §3](viewer/_spec/agent.md)).
- Policies: steps made unsafe by a change outside the API when undo is off; the domain rules listed in Settings; data
  quoted in tool results ("Loaded 1CBS.") as it is in the state block ([agent.md §4](viewer/_spec/agent.md)).
- What the agent sees: colors, snapshots and the canvas in the state block; a state version for repeated reads
  ([agent.md §5](viewer/_spec/agent.md)).
- Without a page: the eval's run log replayed in the app; an MCP server over the headless viewer
  ([agent.md §6](viewer/_spec/agent.md)).
- The setup: models asked from each provider's `/models`, an own endpoint's models listed
  (`assistant/src/viewer/settings-model.ts`).
- The Assistant tab always says stored keys are scrambled; `KeyStore.protection()` (protected by a per-browser key, or
  only scrambled) is not shown (`assistant/src/views/ai-setup.tsx`).
- The command menu's assist: every request is a turn; reading words into one command (`assistant/src/viewer/ui.tsx`).
- Web and data tools for the agent in the viewer (search, fetch a page or a file; `viewer/src/agent/skills.ts`).

## Quadratic paths still open (the 2026-10-17 scan)

The scan's reachable hits fixed and guarded are in `test:scaling` (AGENTS.md, bulk operations). Open, by reach:

- **Mol\*'s particle hierarchy** (`managers.particles.hierarchy.current`, `packages/plugin`) is rebuilt over the whole
  state tree on its first read after a change (about 82 ms over 64 loads); a Simularium frame step reads it (1.4 ms with
  20 structures, 2026-10-21).
- **Particles in `PluginContext`**: `representation.particles` and `managers.particles` are built whatever the registry
  composes (their code imported); they move into `ParticlesRegistry` with the PluginContext rewrite.
- **`representation.set` across many components** (an `api.group` of one call per component) is one commit per
  component; only the updates inside one call are batched.
- **FIFO retirement in a projection** (`ui/src/reactive/transaction.ts`): retiring K in order refolds what is left each
  time (K²); inherent to the fold as it is.

## Considered and deferred

- **One file for the scene and its sequences** (David, 2026-10-24): MolViewSpec and SeqViewSpec together, each
  structure's sequence (its tracks and chains, [sequence.md §7.1](viewer/_spec/sequence.md)) bound to its MVS node by
  ref, not matched by address as `load.svs` does (§6). Later; until then the viewer loads one or the other.
- `agent.skills` as a function of the defaults: `{ add, remove }` was chosen (data, composes by id, needs no import of
  the defaults; [viewer-spec.md §2.4](viewer/_spec/viewer-spec.md)).
- Keys in `assistant(options)`: a page never passes a key; keys are the person's, typed into the Assistant tab
  ([assistant/\_spec](assistant/_spec/README.md), "Keys"). Dev keys from `.env.local` serve only the live test runs
  (`pnpm test:live`); a host with its own models passes `provided` (a server-side adapter).
- Field edges in the default themes: kept as they are (2026-10-16; the fill marks a field, the high-contrast theme and
  `prefers-contrast: more` draw an edge: [styling.md](ui/_spec/styling.md)).
- Compiled row templates for the `dom` runtime (rows cloned rather than created; what creating 10,000 rows still costs
  is element creation, ~50 ms): it needs a build step the runtime does not have
  (`ui/src/_stories/runtime-bench.stories.tsx`).
- Named parts for host styling (`data-msp-part`) and a theming page in Storybook: hosts theme through the `--msp-*`
  tokens only (`ui/src/styles/theme.css`).
