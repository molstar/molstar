# Baseline: where MVS load time goes

Measured 2026-10-18 on `v6-mvs-load` (`0050f1464`), Node 22, Apple Silicon, no canvas.

## Fixtures

`packages/mvs/runtime/src/_test/fixtures/*.mvsj`: the `index.mvsj` of three real MVSX files (assets left out; nothing
is downloaded because the bench never applies a snapshot).

| Fixture  | Snapshots | MVS nodes | `interpolate` nodes | Time steps (1 + duration / frame time) | Distinct frames | Mol\* transforms / frame |
| -------- | --------: | --------: | ------------------: | -------------------------------------: | --------------: | -----------------------: |
| features |        48 |       505 |                  14 |                                  6,695 |           6,447 |                     14.3 |
| motm1    |         5 |       156 |                  25 |                                  7,746 |           2,718 |                     40.6 |
| exosome  |         5 |       573 |                  94 |                                  7,983 |           2,886 |                     75.1 |

What the animations interpolate: motm1 — opacity, label opacity, emissive (`custom.molstar_representation_params`),
translation, a `transform_matrix`, `model_index`; exosome — 57× emissive, 35× clip `point.0`, opacity; features — one
long rotation (86 s, 5,161 steps of a 14-node tree), colors, camera vectors.

## Method

`packages/mvs/runtime/src/_test/load-animation.bench.test.ts` (skipped unless `MVS_BENCH=1`):

```sh
MVS_BENCH=1 pnpm test:unit packages/mvs/runtime/src/_test/load-animation.bench.test.ts --silent=false
MVS_BENCH=1 MVS_BENCH_PROFILE=/tmp/mvs-prof pnpm test:unit ...   # also a .cpuprofile per fixture
```

A `PluginContext` with the MVS behavior and no canvas; `loadMVS(..., { keepCamera: true })` with
`State.Snapshots.Apply` swallowed. Best of 3, allocation = `v8` `total_allocated_bytes` of the first run, retained =
heap after GC with the entries in the snapshot manager.

## Numbers

| Fixture  | `generateStateTransition` | `convertMvsToMolstar` (frames only) | `loadMVS` total | Allocated | Retained |
| -------- | ------------------------: | ----------------------------------: | --------------: | --------: | -------: |
| features |           90 ms (247 MB) |                      75 ms (228 MB) |          587 ms |   1.29 GB |  23.5 MB |
| motm1    |          140 ms (312 MB) |                      67 ms (183 MB) |          736 ms |   1.40 GB |  25.4 MB |
| exosome  |          396 ms (795 MB) |                     196 ms (560 MB) |        1,803 ms |   3.27 GB |  77.0 MB |

Frames differ little from their snapshot: of the 75 transforms of an exosome frame, 8 have params that are not the
base snapshot's (motm1: 3.6 of 41; features: 2.6 of 14). About 90% of the per-frame work rebuilds what the base
snapshot already has.

## Profile (exosome, one `loadMVS`, ~1.9 s sampled)

Inclusive, by stage:

| Stage                                                   |  Share | Notes                                                                                     |
| ------------------------------------------------------- | -----: | ----------------------------------------------------------------------------------------- |
| `molstarTreeToEntry` → `loadTreeVirtual` (per frame)    |   57 % | rebuilding the whole Mol\* state tree of every frame                                      |
| ↳ `StateBuilder.apply` + `StateTransform.create`        |   15 % | immutable.js maps, a `UUID.create22()` version per transform that the hash then overwrites |
| ↳ `TargetManager.getChildRef` → `stringHash`            |   12 % | the same deterministic ref recomputed every frame; `split('').reverse().join('')`         |
| ↳ `getTree({ useHashVersion })` → `hashMurmur128o`      |   13 % | hashes the params of every transform of every frame                                      |
| ↳ loading actions, `dfs`, tags, `toJSON`                |  ~15 % | `representationProps`, `clippingForNode`, `makeNearestReprMap`, …                         |
| `generateStateTransition` (`produce` = mutative)        |   19 % | proxies; `select` through drafts alone 5 % self                                           |
| `convertMvsToMolstar` (per frame)                       |    9 % | `addDefaults`, `convertTree`, `condenseTree`: three full copies of the tree               |
| `UUID.create22` (self, all callers)                     |    8 % |                                                                                           |
| `reuseTransformParams`                                  |    1 % | rebuilds the base snapshot's ref map for every frame                                      |

motm1 has the same shape (61 % `loadTreeVirtual`, 15 % mutative, 7 % conversion).

## Not measured here

- **The canvas**: with a canvas, every frame also builds `modifyCanvasProps(plugin.canvas3d.props, …)` (a copy of the
  canvas props) and a camera snapshot, even when the animation has neither `include_canvas` nor `include_camera`; both
  are kept in the entry the frame is built from. To measure in the browser (Prism viewer or the classic app, Performance
  panel) before and after.
- **Playback**: each frame is applied with `State.setSnapshot` → `StateTree.fromJSON` of the whole frame and a diff by
  `version`.
- **Session size**: transitions are serialized frame by frame (`getStateSnapshot` keeps `transition`); a session of
  exosome holds ~2,900 full state trees.
