# MolViewSpec runtime: the specification

`@molstar/mvs` (`packages/mvs/runtime/`): loads MolViewSpec documents (MVSJ, MVSX) into a Mol\* plugin, as plugin state
snapshots. The format itself, its trees and the builder are `@molstar/mvs-builder` (`packages/mvs/builder/`); the
format's documentation for authors is external (molstar.org/mol-view-spec-docs). This spec is for working on the runtime
without reading all of it. Code comments cite its sections (`animation.md §2`). The code is the truth for how. Being
written along the load optimizations (`.v6/mvs-optimizations/`): files are added as their area is reworked.

| File                             | What it covers                                                                                  |
| -------------------------------- | ----------------------------------------------------------------------------------------------- |
| [load.md](load.md)               | the pipeline per snapshot and frame, refs, the virtual builder, reuse and replay between frames |
| [animation.md](animation.md)     | expanding a snapshot's animation into frames: the update plan, frame dedup, identity            |
| [performance.md](performance.md) | the fixtures, the goldens, the bench, the numbers to keep                                       |

## The load pipeline

`loadMVS(plugin, data, options)` (`load.ts`), per snapshot:

```
validate (MVS and animation tree schemas, MolQL)
convertMvsToMolstar (tree/molstar/conversion.ts)  MVS tree → Molstar tree: defaults, URIs, split nodes, condense
molstarTreeToEntry → loadTreeVirtual (load-generic.ts)  Molstar tree → serialized state tree, never committed; camera, canvas
assignStateTransition                             the animation: generateStateTransition (helpers/animation.ts)
                                                  then the two steps above for every frame → snapshot.transition
```

then the entries replace (or are appended to) the snapshot manager's and one is applied.

- Transform refs are deterministic: `!mvs:<hash of parent ref and transformer>:<n>` (`TargetManager`), so the same node
  gets the same ref in every snapshot and frame, and applying a snapshot updates rather than rebuilds.
- Versions are hashes of the params: equal params, equal version, no update when a frame is applied. Frames reuse the
  transforms they share with the previous frame or their snapshot ([load.md](load.md) §3).

## Entries

| Module                               | What it is                                                                    |
| ------------------------------------ | ----------------------------------------------------------------------------- |
| `load.ts`: `loadMVS`                 | the pipeline above; the loading actions per Molstar node kind                 |
| `load-generic.ts`: `loadTreeVirtual` | a tree → a state snapshot through `UpdateTarget` (`VirtualStateBuilder`)      |
| `behavior.ts`: `MolViewSpec`         | the plugin behavior: formats (MVSJ, MVSX), drag and drop, the providers       |
| `registry.ts`: `MVSRuntimeRegistry`  | the representations and themes MVS needs, for a plugin with an empty registry |
| `helpers/animation.ts`               | `generateStateTransition`                                                     |
