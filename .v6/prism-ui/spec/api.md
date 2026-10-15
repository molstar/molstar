# The viewer's API: the managers and builders done again

`viewer.api` is what scripts, embedders, view models, commands, tools and the agent's code all call: one object per
area, each with its state and its operations as methods. It is the plugin's `managers` and `builders` done again
(David, 2026-10-15), with what v6 learned: every change goes through one run path, takes targets, returns a result,
and is watched through rxjs states. [operations.md](operations.md) defines an operation and lists the areas; this file
defines the objects they make up, the types around them (review F24), grouping, extensions, the code skill and what is
promised stable. Status: **accepted** (David, 2026-10-15).

## 1. At a glance

```ts
const viewer = await createPrismViewer({ target, spec });

await viewer.api.load.entry({ source: 'pdb', id: '1tqn' });          // loaders feature
await viewer.api.structure.preset({ preset: 'polymer-cartoon' });
await viewer.api.representation.color({ target: 'ligand', color: 'element-symbol' });
await viewer.api.select.set({ target: 'A:120-130', mode: 'set' });
await viewer.api.camera.frame({ target: 'selection' });

viewer.api.select.current.subscribe((s) => console.log(s.summary));    // state: a ReadonlyState
const r = await viewer.api.component.add({ target: 'resn HEM', label: 'Heme' });
if (!r.ok) console.warn(r.message, r.hint);                             // results, never throws

await viewer.api.group('Prepare the figure', async (api) => {          // one history step, one record group
  await api.canvas.background({ color: 'white' });
  await api.styles.apply({ style: 'illustrative' });
});
```

## 2. One object per area

| `viewer.api.…`   | Today's plugin underneath                                                          | Brought by        |
| ---------------- | ---------------------------------------------------------------------------------- | ----------------- |
| `load`           | `builders.data`, the loaders (`@molstar/plugin-extension/loaders`), state actions  | core; `loaders`   |
| `session`        | `state` (sessions), the remote-states client (today in the classic UI's view)      | core; `remote-states` |
| `structure`      | `managers.structure.hierarchy`, `builders.structure` (hierarchy presets)           | core              |
| `contents`       | the structure's model (entities, chains), the components it hides                  | core              |
| `component`      | `managers.structure.component`                                                     | core              |
| `representation` | `managers.structure.component` (representations, themes), the selection layers     | core              |
| `styles`, `annotations` | representation presets, quick styles; what extensions register to apply      | core              |
| `volume`, `particles` | `managers.volume`, `managers.particles`                                       | core              |
| `select`, `focus`, `highlight` | `managers.structure.selection`, `.focus`, `managers.interactivity`   | core              |
| `measure`        | `managers.structure.measurement`                                                   | `measurements`    |
| `superpose`      | the superposition algorithms (today in the classic UI's view)                      | `superposition`   |
| `snapshot`       | `managers.snapshot` (its list only; state.md §5)                                   | `snapshots`       |
| `animation`, `trajectory` | `managers.animation`, the model-index transforms                          | core; `animation` |
| `camera`, `canvas` | `managers.camera`, `canvas3d` props, illumination, XR                            | core; `canvas-settings` |
| `export`         | `helpers.viewportScreenshot`, session save; extensions' exports                    | `export-*`        |
| `sequence`       | the Seq\* viewer (sequence.md)                                                     | core              |
| `state`          | `state.data` (the tree, applicable actions)                                        | `state-tree`      |
| `settings`       | `config`, behavior params, the viewer's options                                    | core              |
| `history`, `view` | prism-ui's `history`, `navigation`                                                | core              |

The object's name is its area's id, so `viewer.api.component.add` is the operation `component.add` (operations.md §1):
one name in the API, the command menu, the tools, the run log and the docs.

## 3. The shape every area shares

```ts
abstract class ApiArea {                                    // packages/prism/viewer/src/api/area.ts
  abstract readonly id: string;                             // 'component'
  // state: ReadonlyStates read from the sources (state.md §4); happenings: an `events` object (README §3.6)
  protected operation<P extends PD.Params, R>(def: OperationDefinition<P, R>): OperationFn<P, R>;
}

interface OperationFn<P extends PD.Params, R> {
  (values?: Partial<PD.ValuesFor<P>>, options?: RunOptions): Promise<ApiResult<R>>;
  readonly definition: OperationDefinition<P, R>;           // for commands, tools, coverage, the generated typings
}

interface RunOptions { signal?: AbortSignal; copy?: boolean }   // copy: operations.md §1, identity
```

An area declares its operations as fields, so the method and its definition are one thing:

```ts
class ComponentApi extends ApiArea {
  readonly id = 'component';
  readonly list: ReadonlyState<readonly ComponentInfo[]>;   // from the hierarchy source
  readonly add = this.operation({
    id: 'component.add', title: 'Add component',
    params: { target: TargetParam(), label: PD.Text('') },
    identity: (v) => [v.structure, canonicalTarget(v.target)],
    async run(v, ctx) { /* ctx.targets.resolve(v.target) → builders.structure… */ },
  });
  readonly remove = this.operation({ … });
}
```

- **Values are partial**: a missing value is the param's default, so `add({ target: 'resn HEM' })` is a whole call.
- **Results, never throws** (operations.md §1). The code skill's sandbox throws an `ApiError { code, message, hint }`
  instead, as Stories' story code does (§7).
- **Reads** are the states (for views and scripts) and `read: true` operations (`structure.get`, `select.get`), which
  return the same thing as text and data for commands and tools.

## 4. The types around it

```ts
interface ApiContext {                                      // what an operation's run gets
  readonly plugin: PluginContext;                           // only api/ and sources/ use it (contexts.md §7)
  readonly sources: ViewerSources;                          // read-only states over the plugin (state.md §4)
  readonly targets: TargetResolver;                         // resolve(target, structure?) → loci and "what it selected" (targets.md)
  readonly api: ViewerApi;                                  // other areas, bound to the same origin
  readonly origin: RunSource;                               // person, agent (with its turn), script, url, link
  readonly modules: ApiModules;                             // what the host offers (below)
  readonly signal: AbortSignal;
}

type ApiModuleKey = 'canvas' | 'files' | 'network' | 'storage' | 'page';
interface ApiModules { has(key: ApiModuleKey): boolean }    // an operation's `modules` are checked against it

interface Projection<T = unknown> {                         // an optimistic edit (managers.transactions)
  readonly model: string;                                   // which view model's base it projects onto
  apply(base: T): T;
  isAbsorbed(base: T): boolean;                             // the source reported it: drop the projection
}

type ViewerApi = CoreApi & ExtensionApis;                   // CoreApi: the areas of §2; ExtensionApis: §6
```

## 5. Origins, grouping and builders

- **Bound by origin**: `viewer.api` runs as `script`. Each surface gets the same objects bound to its origin:
  `viewer.api.as({ kind: 'person' })` for view models and commands, `{ kind: 'agent', turn }` for tools, `url` for the
  loaders feature, `link` for the sequence. The history and the record read the origin from there.
- **`viewer.api.group(label, work)`**: the calls inside are one history step and one group in the record
  (sequence.md §7.3); a failure inside restores the group's checkpoint, so it is all or nothing. Commands that run
  several operations, scripts and the code skill use it.
- **Builders, the declarative half**: a whole scene is built declaratively with MolViewSpec (`@molstar/mvs-builder`)
  and loaded with `load.mvs({ data })`, one step; the plugin's own `builders` stay inside operations. No second
  builder language.
- **Raw plugin work**: `viewer.run(label, (plugin) => …)` (README §3.2) for what no operation covers yet, as one step;
  what it does is recorded as one opaque run.

## 6. Extensions

- **An extension's areas** come from `ViewerExtension.api` (viewer-spec.md §3) and sit under its id:
  `viewer.api['geo-export'].export({ format: 'glb' })`, operation `geo-export.export`.
- **Typed when imported**: an extension augments the API's types, so a script that imports it gets them:

  ```ts
  declare module '@molstar/prism-viewer/api' {
    interface ExtensionApis { 'geo-export': GeoExportApi }
  }
  ```

- **Lazy extensions**: their areas exist from the start as stubs from `stub.operations`; calling one loads the chunk,
  then runs the real operation (the call is replayed; a failed load is the result).
- **Removed or overridden**: an area whose feature is left out is absent (`viewer.api.measure` is `undefined` without
  `measurements`); commands and tools over it say why they cannot run (viewer-spec.md §4).

## 7. Code: the agent and scripts writing against the API

- **The `viewer-code` skill** (agent.md §3): `run_viewer_code { code }` runs TypeScript against `ViewerApi` in a
  sandbox, on Stories' toolscript (`any-toolscript`: a null-origin iframe in the page, QuickJS in Node; the executor
  convention `async function executor(api) { … }`; limits on time, calls, output).
- **The typings the model reads are generated** from the API's definitions (the areas, their params as TypeScript,
  their states as read-only data), as compact as the MVS builder's, with the prelude that builds the object model in
  the sandbox; only plain data crosses.
- **All or nothing**: a snippet runs in a group (§5), bound to the agent's origin and turn; a throw, a limit or a type
  error restores the group's checkpoint and says what and where. The result the model reads: what the snippet returned,
  what it printed, its calls in one line each (their `did`), what changed.
- The same typings are the API's reference for people (`apps/docs`).

## 8. What is promised

- **The core areas are `@molstar/prism-viewer`'s public contract**: their ids, params and results follow semver once the
  package is released (the operation ids are also what the run log, kept records and the agent's skills name).
  Extensions' areas carry no such promise (README decision 13).
- **The plugin stays reachable** (`viewer.plugin`) and is not part of the promise.
- **v5's `Viewer` methods** map onto areas (the `load*` methods onto `load`, `handleResize` onto nothing: the shell
  sizes the canvas, `structureInteractivity` onto `highlight` and `select` with `events.hover`/`click`), listed in
  `.v6/plans/breaking-v6-changes.md` when the viewer ships.

## 9. Where

```
packages/prism/viewer/src/api/
  area.ts            ApiArea, OperationFn, RunOptions, the run path's entry
  context.ts         ApiContext, ApiModules, TargetResolver binding
  index.ts           createViewerApi(ctx): the core areas; extensions' areas added by the composition
  <area>.ts          one file per area of §2, its operations as fields, its states from sources/
  typings.ts         the generated typings and prelude for the code skill (§7)
```

## 10. Open

1. Whether a script's `viewer.api` should default to `person` when called from the page's own UI code (today: `script`,
   and the UI binds `person` explicitly).
2. Redo (sequence.md §7.3): left out until checkpoints can be kept after a step too.
