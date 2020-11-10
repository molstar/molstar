/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '../../../mol-plugin/commands';
import { StateTransform } from '../../../mol-state';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateAnimation } from '../model';

export const AnimateStateInterpolation = PluginStateAnimation.create({
    name: 'built-in.animate-state-interpolation',
    display: { name: 'Animate State Interpolation' },
    params: () => ({
        transtionDurationInMs: PD.Numeric(2000, { min: 100, max: 30000, step: 10 })
    }),
    canApply(plugin) {
        return { canApply: plugin.managers.snapshot.state.entries.size > 1 };
    },
    initialState: () => ({ }),
    async apply(animState, t, ctx) {

        const snapshots = ctx.plugin.managers.snapshot.state.entries;
        if (snapshots.size < 2) return { kind: 'finished' };

        const currentT = (t.current % ctx.params.transtionDurationInMs) / ctx.params.transtionDurationInMs;

        let srcIndex = Math.floor(t.current / ctx.params.transtionDurationInMs) % snapshots.size;
        let tarIndex = Math.ceil(t.current / ctx.params.transtionDurationInMs);
        if (tarIndex === srcIndex) tarIndex++;
        tarIndex = tarIndex % snapshots.size;

        const _src = snapshots.get(srcIndex)!.snapshot, _tar = snapshots.get(tarIndex)!.snapshot;

        if (!_src.data || !_tar.data) return { kind: 'skip' };

        const src = _src.data.tree.transforms, tar = _tar.data.tree.transforms;

        const state = ctx.plugin.state.data;
        const update = state.build();

        for (const s of src) {
            for (const t of tar) {
                if (t.ref !== s.ref) continue;
                if (t.version === s.version) continue;

                const e = StateTransform.fromJSON(s), f = StateTransform.fromJSON(t);

                if (!e.transformer.definition.interpolate) {
                    update.to(s.ref).update(currentT <= 0.5 ? e.params : f.params);
                } else {
                    update.to(s.ref).update(e.transformer.definition.interpolate(e.params, f.params, currentT, ctx.plugin));
                }
            }
        }

        await PluginCommands.State.Update(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });

        return { kind: 'next', state: { } };
    }
});