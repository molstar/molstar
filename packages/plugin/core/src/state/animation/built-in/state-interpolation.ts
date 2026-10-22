/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '@molstar/plugin/commands';
import { StateTransform } from '@molstar/core/state';
import { deepEqual } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginStateAnimation } from '../model.js';

export const AnimateStateInterpolation = PluginStateAnimation.create({
  name: 'built-in.animate-state-interpolation',
  display: { name: 'Animate State (experimental)' },
  params: () => ({
    transtionDurationInMs: PD.Numeric(2000, { min: 100, max: 30000, step: 10 }),
  }),
  canApply(plugin) {
    return { canApply: plugin.managers.snapshot.state.entries.size > 1 };
  },
  initialState: () => ({ last: new Map<string, unknown>() }),
  async apply(animState, t, ctx) {
    const last = animState.last ?? new Map<string, unknown>();
    const snapshots = ctx.plugin.managers.snapshot.state.entries;
    if (snapshots.size < 2) return { kind: 'finished' };

    const currentT = (t.current % ctx.params.transtionDurationInMs) / ctx.params.transtionDurationInMs;

    const srcIndex = Math.floor(t.current / ctx.params.transtionDurationInMs) % snapshots.size;
    let tarIndex = Math.ceil(t.current / ctx.params.transtionDurationInMs);
    if (tarIndex === srcIndex) tarIndex++;
    tarIndex = tarIndex % snapshots.size;

    const _src = snapshots.get(srcIndex)!.snapshot,
      _tar = snapshots.get(tarIndex)!.snapshot;

    if (!_src.data || !_tar.data) return { kind: 'skip' };

    const src = _src.data.tree.transforms,
      tar = _tar.data.tree.transforms;

    const state = ctx.plugin.state.data;
    const update = state.build();

    // TODO: support for adding/removing nodes
    const targets = new Map<string, StateTransform.Serialized>();
    for (const t of tar) targets.set(t.ref, t);

    for (const s of src) {
      const t = targets.get(s.ref);
      if (!t || t.version === s.version) continue;

      const cell = state.cells.get(s.ref);
      if (!cell) continue;

      const e = StateTransform.fromJSON(s),
        f = StateTransform.fromJSON(t);

      const newParams = e.transformer.definition.interpolate
        ? e.transformer.definition.interpolate(e.params, f.params, currentT, ctx.plugin)
        : currentT <= 0.5
          ? e.params
          : f.params;

      // the cell's params (normalized), or what this animation set last (as in the snapshot), not the cell itself: a
      // tick that leaves them as they are updates nothing (S8)
      if (deepEqual(last.get(s.ref), newParams) || deepEqual(cell.transform.params, newParams)) continue;
      update.to(s.ref).update(newParams);
      last.set(s.ref, newParams);
    }

    if (update.editInfo.count > 0) {
      await PluginCommands.State.Update(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });
    }

    return { kind: 'next', state: { last } };
  },
});
