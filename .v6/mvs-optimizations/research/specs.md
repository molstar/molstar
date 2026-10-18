# `_spec` folders for MVS

Neither MVS package has a `_spec` yet (`packages/mvs/runtime/src/README.md` points to the external docs). Written
alongside the optimization, each file when its area is touched, in the style of the other specs (terse, rule + reason,
code comments cite sections).

## `packages/mvs/builder/_spec/` (`@molstar/mvs-builder`)

| File             | What it covers                                                                                           |
| ---------------- | -------------------------------------------------------------------------------------------------------- |
| `README.md`      | the package: MVS data (`MVSData`, states, snapshots, MVSJ/MVSX), the builder, the tree schemas; entries  |
| `tree.md`        | generic trees: schema, params schema and defaults, `addDefaults`, `convertTree`, `condenseTree`, `dfs`, validation; identity rules (which functions copy, which mutate, what the conversion memo may rely on) |
| `animation.md`   | the animation tree: `interpolate` kinds and params, `start_ms` / grouping / previous-end semantics, frequency and alternation, noise, `frame_time_ms`, `include_camera`/`include_canvas` |

## `packages/mvs/runtime/_spec/` (`@molstar/mvs`)

| File             | What it covers                                                                                           |
| ---------------- | -------------------------------------------------------------------------------------------------------- |
| `README.md`      | the package, the load pipeline in one picture, entries (`loadMVS`, behavior, registry), where the code is |
| `load.md`        | validate → expand → convert → load virtual → entries; deterministic refs (`!mvs:hash:n`), hash versions, `replaceExisting`; the loading actions, their context and their locality (the unit of reuse, molstar-conversion.md step 4) |
| `animation.md`   | expansion: the compiled update plan, slot equality and frame dedup, identity across frames, what a frame holds and shares with the base and the previous frame |
| `testing.md` | budgets per fixture (time, allocation, retained), the bench and how to read it, the goldens               |

## Questions to settle while writing them

- Are precomputed frames a requirement (sessions must replay without the MVS source) or an implementation choice?
  Decides step 5 of molstar-conversion.md.
- Should a frame carry camera/canvas only with `include_camera`/`include_canvas`? Today the frame's entry computes
  both and the transition drops them; nothing else reads them (to confirm).
- `deepEqual` dedup treats two time steps with equal trees as one frame of summed duration. With noise every step is
  distinct; is a frame per step intended there (5,000+ frames for a long noisy animation)?
- `interpolate` on `model_index` (motm1) changes a `model` node's params, which the conversion splits from `structure`:
  the conversion memo must not key `model` on the MVS node alone. Any other interpolatable param that changes the shape
  of the converted tree (a `discrete` scalar picking a format or a branch)? None found; the frame emitter can assert the
  shape (the list of refs) is the base's and fall back to a full rebuild.
