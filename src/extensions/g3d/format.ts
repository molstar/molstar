/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TrajectoryFormatProvider, TrajectoryFormatCategory } from '../../mol-plugin-state/formats/trajectory';
import { PluginStateTransform, PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { G3dHeader, getG3dHeader, getG3dDataBlock } from './data';
import { Task } from '../../mol-task';
import { PluginContext } from '../../mol-plugin/context';
import { ParamDefinition } from '../../mol-util/param-definition';
import { trajectoryFromG3D } from './model';
import { StateObjectRef, StateAction } from '../../mol-state';
import { PluginBehavior } from '../../mol-plugin/behavior';

export const G3dProvider: TrajectoryFormatProvider = {
    label: 'G3D',
    description: 'G3D',
    category: TrajectoryFormatCategory,
    binaryExtensions: ['g3d'],
    parse: async (plugin, data) => {
        const trajectory = await plugin.state.data.build()
            .to(data)
            .apply(G3DHeaderFromFile, {}, { state: { isGhost: true } })
            .apply(G3DTrajectory)
            .commit();

        return { trajectory };
    },
    visuals: defaultVisuals
};

async function defaultVisuals(plugin: PluginContext, data: { trajectory: StateObjectRef<SO.Molecule.Trajectory> }) {
    const builder = plugin.builders.structure;
    const model = await builder.createModel(data.trajectory);
    const modelProperties = await builder.insertModelProperties(model);
    const structure = await builder.createStructure(modelProperties);
    const all = await builder.tryCreateComponentStatic(structure, 'all');

    if (all) {
        await builder.representation.addRepresentation(all, {
            type: 'cartoon',
            color: 'polymer-index',
            size: 'uniform',
            sizeParams: { value: 2 }
        });
    }
}

export class G3dHeaderObject extends SO.Create<{
    header: G3dHeader,
    urlOrData: Uint8Array | string
}>({ name: 'G3D Header', typeClass: 'Data' }) { }

export type G3DHeaderFromFile = typeof G3DHeaderFromFile
export const G3DHeaderFromFile = PluginStateTransform.BuiltIn({
    name: 'g3d-header-from-file',
    display: { name: 'G3D Header', description: 'Parse G3D Header' },
    from: SO.Data.Binary,
    to: G3dHeaderObject
})({
    apply({ a }, plugin: PluginContext) {
        return Task.create('Parse G3D', async () => {
            const header = await getG3dHeader(plugin, a.data);
            return new G3dHeaderObject({ header, urlOrData: a.data }, { label: header.name, description: header.genome });
        });
    }
});

export type G3DHeaderFromUrl = typeof G3DHeaderFromUrl
export const G3DHeaderFromUrl = PluginStateTransform.BuiltIn({
    name: 'g3d-header-from-url',
    display: { name: 'G3D Header', description: 'Parse G3D Header' },
    params: { url: ParamDefinition.Text('') },
    from: SO.Root,
    to: G3dHeaderObject
})({
    apply({ params }, plugin: PluginContext) {
        return Task.create('Parse G3D', async () => {
            const header = await getG3dHeader(plugin, params.url);
            return new G3dHeaderObject({ header, urlOrData: params.url }, { label: header.name, description: header.genome });
        });
    }
});

export type G3DTrajectory = typeof G3DHeaderFromUrl
export const G3DTrajectory = PluginStateTransform.BuiltIn({
    name: 'g3d-trajecotry',
    display: { name: 'G3D Trajectory', description: 'Create G3D Trajectory' },
    params: a => {
        if (!a) return { resolution: ParamDefinition.Numeric(200000) };
        const rs = a.data.header.resolutions;
        return {
            resolution: ParamDefinition.Select(rs[rs.length - 1], rs.map(r => [r, '' + r] as const))
        };
    },
    from: G3dHeaderObject,
    to: SO.Molecule.Trajectory
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('G3D Trajectory', async ctx => {
            const data = await getG3dDataBlock(plugin, a.data.header, a.data.urlOrData, params.resolution);
            const traj = await trajectoryFromG3D(data).runInContext(ctx);
            return new SO.Molecule.Trajectory(traj, { label: a.label, description: a.description });
        });
    }
});

export const LoadG3D = StateAction.build({
    display: { name: 'Load Genome 3D (G3D)', description: 'Load G3D file from the specified URL.' },
    from: SO.Root,
    params: { url: ParamDefinition.Text('') }
})(({ params, state }, ctx: PluginContext) => Task.create('Genome3D', taskCtx => {
    return state.transaction(async () => {
        if (params.url.trim().length === 0) {
            throw new Error('Specify URL');
        }

        ctx.behaviors.layout.leftPanelTabName.next('data');

        const trajectory = await state.build().toRoot()
            .apply(G3DHeaderFromUrl, { url: params.url })
            .apply(G3DTrajectory)
            .commit();

        await defaultVisuals(ctx, { trajectory });
    }).runInContext(taskCtx);
}));

export const G3DFormat = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'g3d',
    category: 'misc',
    display: {
        name: 'G3D',
        description: 'G3D Format Support'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register() {
            this.ctx.state.data.actions.add(LoadG3D);
        }
        unregister() {
            this.ctx.state.data.actions.remove(LoadG3D);
        }
    }
});