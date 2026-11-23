# Viewer composition: a spec, as Mol\*'s plugins compose

## Intent

The author, 2026-11-17: "let's be explicit in designs about how the "viewer" will compose? it should be similar to how
existing "plugins" already composes, and it should configure all the ports, asset mangers, etc".

A Seq\* viewer is composed the way a Mol\* v6 plugin is: from an explicit spec that lists what it contains, created by
one call, with nothing registered on import. One composition root names every registry entry, behavior, setting and
port adapter, so what a viewer can do and reach is visible in one place.

## Precedents

- Mol\* v6's plugin composition (`.v6/designs/plugin-composition.md`; `PluginSpec` in `@molstar/plugin/spec`):
  registries start empty; `spec.registry` lists declarative entries (formats, representations, themes, presets, …)
  registered in order, without dependency resolution; registration is atomic (every conflict found before anything
  changes, one error listing them), counted (the same object again increments, a different object under a key is an
  error) and undoable (`register` returns its undo); `behaviors` with params; `config` as `[item, value]` pairs;
  `assetManager` as a factory with ownership (`disposeWithPlugin`); `DefaultPluginSpec` for full defaults; base entry
  points take an explicit spec.
- The previous implementation: a viewer spec shaped like the plugin spec (registry, behaviors, config, render) and
  `createSeqViewer(spec, ports)` ([sources/seqstar](../sources/seqstar/README.md)).

## Decision

### The spec

```ts
interface SeqViewerSpec {
  /** Registry entries, registered in order: what the viewer can read, compute and draw. */
  registry?: readonly SeqRegistryEntry[];
  /** Behaviors with their params: the viewer's policies (hover, selection, navigation, follow, …). */
  behaviors: readonly SeqViewerSpec.Behavior[];
  /** The interaction grammar: a preset and binding overrides (ux-foundations.md, section 8). */
  bindings?: SeqViewerSpec.Bindings;
  /** Settings as [item, value] pairs, as the plugin's config: foundations' preset and numbers, limits, defaults. */
  config?: readonly [SeqConfigItem, unknown][];
  /** Port adapters, each a factory with ownership; unset ports take their default adapter. */
  ports?: SeqViewerSpec.Ports;
  /** The default UI's parts and their places; none for a headless viewer. */
  ui?: SeqViewerSpec.UI;
}
```

- **Registry entries are pure data**, as the plugin's: a record of providers grouped by kind, no behaviors, no config,
  no functions run at registration. Kinds in this version: node kinds, layer types, themes and palettes, residue
  schemes, numbering schemes and codon tables, format readers and writers, map runtimes, backends, query forms and
  properties (`query-form`, `query-property`: decided by the author, 2026-11-20, "make the language extensible so people
  could add generators, filters, etc on their own props"; [locations-and-queries.md](locations-and-queries.md),
  [properties.md](properties.md)), commands. Layer types and section layouts are registry entries (`layer-type`, `layout`; MV, 2026-11-19): the canvas and `renderSvg`/`exportSvg` draw the composition's own, so a layer type the canvas draws the export draws too, and one neither has is said (`seq.render.layout`); a spec whose check finds an error exports what resolves, as the canvas draws it, the first ten errors in the picture's `<desc>` (the author, 2026-11-22, drift svg-on-error), except a spec past the bounds (`seq.lang.bound.*`), exported blank, since resolving it is what the bounds guard against; tracks is the one fixed fallback layout. Keys are names; extension names are namespaced ([extensibility.md](extensibility.md)). Chrome parts
  are internal interfaces until an outside user needs them; track header actions are the first opened (`headerActions`,
  [track-header-parts](track-header-parts.md)).
- **Ports are factories with ownership**, as `PluginSpec.assetManager`:

  ```ts
  interface Ports {
    assets?: () => ({ manager: AssetManager } | { reader: AssetReader }) & { disposeWithViewer?: boolean };
    textMetrics?: () => { measure: (font: string, text: string) => number; disposeWithViewer?: boolean };
  }
  ```

  Each is called once per viewer. Unset, the viewer makes the default adapter and disposes it with itself; a host that
  passes its own (a Mol\* plugin's asset manager) sets `disposeWithViewer: false` and keeps it ([ports.md](ports.md)).
  The assets port takes a manager, or a bare reader (bytes by location) where a test or a CLI has no manager; the
  policy wraps either. Formats, map runtimes, backends and services are registry entries, since they are providers by name; the ports carry
  what is one per viewer. Time and frames are `@molstar/ui/dom`'s `DomHost`, not a port.
- **Keys are matched against the commands' shortcuts**, scoped to the viewer's container, until a command manager and
  the `@molstar/ui` keymap are wired. Why: the viewer works headless of both, and the match moves into the manager
  without changing the commands.
- **`theme.source`** (`root`, default): a mounted viewer reads its root's CSS; `config` reads the core's `theme` item;
  a headless viewer always uses the config. Why: a mounted viewer follows its page's theme without a second setting.
- **`rested`** moves each time the picture comes to rest: `canvas.settleMs` (150) after the view's last move, or after
  the last layout that changed what is drawn while the view was still (data arriving, the content's height). What
  follows the picture (the SVG beside it) redraws on it. Why: a picture that changed while still would otherwise
  never be redrawn.
- **Behaviors** are as the plugin's: a provider and its params, initialized after the registry, able to register more
  entries at run time through `viewer.register(entry)` and to undo them.
- **Config items** are typed settings described with `ParamDefinition`, as the plugin's `PluginConfigItem`s: one typed
  configuration set at composition. User settings layered over it, the effective-value lookup and a generated settings
  form come later ([extensibility.md](extensibility.md)).

### Creating a viewer

- `createSeqViewer(spec, target?)` makes the viewer: a `ReactiveModel` that creates its managers, binds the ports,
  registers `spec.registry`, then initializes the behaviors, then mounts the UI into `target` if one is given. Without a
  target it is headless: models, commands, rendering to SVG.
- **Registration as the plugin's**: atomic (every conflict found first, one error listing them, the viewer unchanged),
  counted, undoable; registering before the managers exist throws; behaviors can rely on registered providers.
- **Defaults are a spec, not magic**: `DefaultSeqViewerSpec()` (`@molstar/seqx-viewer/default-spec`) lists the built-in
  entries, behaviors, bindings and settings. An application spreads it and changes what it needs, or writes a smaller
  spec; what is not imported and listed is not in its bundle.
- **Injected, opt in** (the author, 2026-11-19: "Injectable, opt-in", over lazy imports; measured by SW: a minimal
  viewer 959 kB against the default's 1023 kB, the unused parts kept in by the renderer's built-in tables and the
  viewer's static imports): the renderer's layer types and layouts come from the registry the spec lists, never from
  a built-in table a `Renderer` defaults to; the viewer's history, story, navigator and rows pieces are behaviors or
  parts the default spec lists, never imported by `createSeqViewer` itself; barrels (`render`, `export`) re-export
  nothing a minimal composition does not ask for; the packages say `sideEffects: false`. `DefaultSeqViewerSpec()` and
  `createSeqContext`'s default list everything, so the default composition is unchanged; a smaller spec ships less.
  Why: one rule for every part (listed or absent), and the size a host pays follows what it lists. **Lazy chunks**: a
  listed part used rarely or late (the icons' markup, a header's menu and DOM parts, a runtime such as mmCIF, the
  molecule letters) may load by a dynamic `import()` on first use, which every bundler splits; the size script lists
  each chunk (the author, 2026-11-20: "Allow lazy chunks"). The first frame does not wait for one: what it would draw
  is drawn plainly (a slot's letter for its icon) and drawn again when it lands. Checked by the size script's `viewer, minimal` entry (an esbuild metafile: none of the
  unlisted modules in it).
- **Headless without a viewer**: `createSeqContext(spec)` in the core takes the registry, the config and the ports a
  render needs (assets, text metrics, backends) and gives the resolve passes and `renderSvg`, for the CLI, servers and
  tests.
- **A signal named for the picture follows every cause of the picture**: one named for a state of what is drawn
  (rested, drawn, ready) fires on every cause that changes the picture (a move, data arriving, a resize), not only on
  the input that usually causes it; its test drives each cause. Why: a `rested` that followed only the view's moves
  never fired again for a picture changed by data arriving while the view was still, so what follows it (a story's
  SVG drawn again) missed the change.
- **Disposal**: `viewer.retire()` disposes what the viewer made, its owned ports included, and undoes its
  registrations; a host's ports are left to the host.

### Inside a Mol\* plugin

- A Mol\* plugin composes a Seq\* viewer with its own adapters: the plugin's asset manager (downloads shared), its
  host and frames, its command manager (Seq\* commands in the plugin's palette), its settings (config items in the
  plugin's), and, beside the composition, `linkStructure` between the viewer and the plugin's structures (a library, not
  a port: [structure-syncing](structure-syncing.md)). The Mol\* side is a
  plugin behavior that builds the Seq\* spec from the plugin's context; Seq\* imports nothing of the plugin.

## Consequences

- The Storybook's stories each compose a viewer from a spec in view; the Extensibility page shows a small spec, the
  default spec, and a spec bound to a Mol\* plugin's ports.
- A conformance test composes a viewer from the default spec and from an empty one, and checks registration is atomic,
  counted and undone on retire.

## In the specs

| Decision | Where |
|---|---|
| The spec's fields (registry, behaviors, bindings, config, ports, ui) | viewer › Composition |
| Registry entries are pure data, keyed by name, the kinds | core › Registry entries |
| Layer types and layouts are registry entries; tracks the one fixed fallback | core › The tracks layout is the one fixed fallback |
| An export of a spec with errors draws what resolves, the errors in `<desc>`; past the bounds blank | core › Export |
| Ports are factories with ownership (`disposeWithViewer`); a manager or a bare reader | core › Ports |
| The viewer's default port adapters | viewer › Ports |
| Keys matched against the commands' shortcuts until a command manager | viewer › Keys |
| `theme.source` (`root`, `config`; headless reads the config) | viewer › Tokens |
| `rested` and a signal named for the picture follows every cause | viewer › A signal named for the picture |
| Behaviors: a provider and its params, after the registry; run-time registration and its undo | viewer › Extensions |
| Config items as `[item, value]` pairs | viewer › Configuration |
| Config items described with `ParamDefinition`, user settings, a settings form | design-only: not built (the defaults are a plain object; the form comes later) |
| `createSeqViewer(spec, target?)` and its order; headless without a target | viewer › Composition |
| Registration atomic, counted, undoable | core › Registry entries |
| Defaults are a spec (`DefaultSeqViewerSpec`) | viewer › Composition |
| Injected, opt in; lazy chunks | viewer › Injected, opt in |
| Headless without a viewer (`createSeqContext`) | core › Context |
| Disposal (`retire`) | viewer › Retire |
| Inside a Mol\* plugin: a plugin behavior building the spec from the plugin's adapters | design-only: not built (no Mol\* plugin behavior yet) |
| `linkStructure` beside the composition, a library not a port | viewer › Structure link |
