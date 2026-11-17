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
  schemes, numbering schemes and codon tables, format readers and writers, map runtimes, backends, query functions,
  commands. Keys are names; extension names are namespaced ([extensibility.md](extensibility.md)). Layouts and chrome
  parts are internal interfaces until an outside user needs them.
- **Ports are factories with ownership**, as `PluginSpec.assetManager`:

  ```ts
  interface Ports {
    assets?: () => { manager: AssetManager; disposeWithViewer?: boolean };
    textMetrics?: () => { measure: (font: string, text: string) => number; disposeWithViewer?: boolean };
    structure?: () => { runtime: StructureRuntime; disposeWithViewer?: boolean };
  }
  ```

  Each is called once per viewer. Unset, the viewer makes the default adapter and disposes it with itself; a host that
  passes its own (a Mol\* plugin's asset manager) sets `disposeWithViewer: false` and keeps it ([ports.md](ports.md)).
  Formats, map runtimes, backends and services are registry entries, since they are providers by name; the ports carry
  what is one per viewer. Time and frames are `@molstar/ui/dom`'s `DomHost`, not a port.
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
- **Headless without a viewer**: `createSeqContext(spec)` in the core takes the registry, the config and the ports a
  render needs (assets, text metrics, backends) and gives the resolve passes and `renderSvg`, for the CLI, servers and
  tests.
- **Disposal**: `viewer.retire()` disposes what the viewer made, its owned ports included, and undoes its
  registrations; a host's ports are left to the host.

### Inside a Mol\* plugin

- A Mol\* plugin composes a Seq\* viewer with its own adapters: the plugin's asset manager (downloads shared), its
  host and frames, its command manager (Seq\* commands in the plugin's palette), its settings (config items in the
  plugin's), and the Mol\* structure runtime (structure systems bound to the plugin's structures). The Mol\* side is a
  plugin behavior that builds the Seq\* spec from the plugin's context; Seq\* imports nothing of the plugin.

## Consequences

- The Storybook's stories each compose a viewer from a spec in view; the Extensibility page shows a small spec, the
  default spec, and a spec bound to a Mol\* plugin's ports.
- A conformance test composes a viewer from the default spec and from an empty one, and checks registration is atomic,
  counted and undone on retire.
