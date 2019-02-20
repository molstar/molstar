/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateAnimation } from './model';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { StateSelection } from 'mol-state/state/selection';
import { PluginCommands } from 'mol-plugin/command';
import { ParamDefinition } from 'mol-util/param-definition';

export const AnimateModelIndex = PluginStateAnimation.create({
    name: 'built-in.animate-model-index',
    display: { name: 'Animate Model Index' },
    params: () => ({
        direction: ParamDefinition.Select('forward', [['forward', 'Forward'], ['backward', 'Backward']]),
        maxFPS: ParamDefinition.Numeric(3, { min: 0.5, max: 30, step: 0.5 })
    }),
    initialState: () => ({ }),
    async apply(animState, t, ctx) {
        // limit fps
        if (t.current > 0 && t.current - t.lastApplied < 1000 / ctx.params.maxFPS) {
            return { kind: 'skip' };
        }

        const state = ctx.plugin.state.dataState;
        const models = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Model)
            .filter(c => c.transform.transformer === StateTransforms.Model.ModelFromTrajectory));

        const update = state.build();

        const dir = ctx.params.direction === 'backward' ? -1 : 1;

        for (const m of models) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, [PluginStateObject.Molecule.Trajectory]);
            if (!parent || !parent.obj) continue;
            const traj = parent.obj as PluginStateObject.Molecule.Trajectory;
            update.to(m.transform.ref).update(StateTransforms.Model.ModelFromTrajectory,
                old => {
                    let modelIndex = (old.modelIndex + dir) % traj.data.length;
                    if (modelIndex < 0) modelIndex += traj.data.length;
                    return { modelIndex };
                });
        }

        await PluginCommands.State.Update.dispatch(ctx.plugin, { state, tree: update });
        return { kind: 'next', state: {} };
    }
})