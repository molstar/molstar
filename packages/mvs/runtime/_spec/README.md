# MolViewSpec runtime: the specification

`@molstar/mvs` (`packages/mvs/runtime/`) loads MolViewSpec documents (MVSJ, MVSX) into a Mol\* plugin as plugin state
snapshots. The format, its trees and the builder are `@molstar/mvs-builder` (`packages/mvs/builder/`). Code comments
cite sections of these files (`load.md §4`). The code is the truth for how.

| File                         | What it covers                                                                   |
| ---------------------------- | -------------------------------------------------------------------------------- |
| [load.md](load.md)           | the pipeline per snapshot and frame, refs, versions, the virtual builder, replay |
| [animation.md](animation.md) | expanding an animation into frames: semantics, the update plan, color maps       |
| [testing.md](testing.md)     | fixtures, goldens, bench                                                         |

## Pipeline

`loadMVS(plugin, data, options)` (`load.ts`), per snapshot:

```
validate     MVS and animation tree schemas, MolQL
convert      MVS tree → Molstar tree (tree/molstar/conversion.ts)
load         Molstar tree → serialized state tree, never committed (load-generic.ts); camera, canvas
animation    generateStateTransition (helpers/animation.ts): one MVS tree per frame,
             then convert and load each → snapshot.transition
```

The entries replace (or are appended to) the snapshot manager's, and one is applied.

## Entries

| Module                               | What it is                                                                       |
| ------------------------------------ | -------------------------------------------------------------------------------- |
| `load.ts`: `loadMVS`                 | the pipeline; the loading actions per Molstar node kind                          |
| `load-generic.ts`: `loadTreeVirtual` | a tree → a state snapshot (`UpdateTarget`, `VirtualStateBuilder`, `SubtreeMemo`) |
| `tree/molstar/conversion.ts`         | `convertMvsToMolstar`, `convertMvsToMolstarMemo`                                 |
| `helpers/animation.ts`               | `generateStateTransition`                                                        |
| `behavior.ts`: `MolViewSpec`         | the plugin behavior: formats (MVSJ, MVSX), drag and drop, the providers          |
| `registry.ts`: `MVSRuntimeRegistry`  | the representations and themes MVS needs, for a plugin with an empty registry    |
