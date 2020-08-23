/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Trajectory } from '../../mol-model/structure';
import { TrajectoryFormatCategory, TrajectoryFormatProvider } from '../../mol-plugin-state/formats/trajectory';
import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { DefaultQueryRuntimeTable } from '../../mol-script/runtime/query/base';
import { StateAction, StateObjectRef } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition } from '../../mol-util/param-definition';
import { G3dHeader, getG3dDataBlock, getG3dHeader } from './data';
import { g3dHaplotypeQuery, G3dLabelProvider, trajectoryFromG3D, G3dSymbols, G3dInfoDataProperty } from './model';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { createStructureRepresentationParams } from '../../mol-plugin-state/helpers/structure-representation-params';
import { stringToWords } from '../../mol-util/string';
import { objectForEach } from '../../mol-util/object';

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
    visuals: defaultStructure
};

async function defaultStructure(plugin: PluginContext, data: { trajectory: StateObjectRef<SO.Molecule.Trajectory> }) {
    const builder = plugin.builders.structure;
    const model = await builder.createModel(data.trajectory);

    if (!model) return;
    const structure = await builder.createStructure(model);

    const info = G3dInfoDataProperty.get(model.data!);
    if (!info) return;

    const components = plugin.build().to(structure);

    const repr = createStructureRepresentationParams(plugin, void 0, {
        type: 'cartoon',
        color: 'polymer-index',
        size: 'uniform',
        sizeParams: { value: 0.25 }
    });

    for (const h of info.haplotypes) {
        components
            .apply(StateTransforms.Model.StructureSelectionFromExpression, { expression: g3dHaplotypeQuery(h), label: stringToWords(h) })
            .apply(StateTransforms.Representation.StructureRepresentation3D, repr);
    }

    await components.commit();
}

export class G3dHeaderObject extends SO.Create<{
    header: G3dHeader,
    urlOrData: Uint8Array | string,
    cache: { [resolution: number]: Trajectory | undefined }
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
            return new G3dHeaderObject({ header, urlOrData: a.data, cache: { } }, { label: header.name, description: header.genome });
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
            return new G3dHeaderObject({ header, urlOrData: params.url, cache: { } }, { label: header.name, description: header.genome });
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
            if (a.data.cache[params.resolution]) {
                return new SO.Molecule.Trajectory(a.data.cache[params.resolution]!, { label: a.label, description: a.description });
            }
            const data = await getG3dDataBlock(plugin, a.data.header, a.data.urlOrData, params.resolution);
            const traj = await trajectoryFromG3D(data).runInContext(ctx);
            a.data.cache[params.resolution] = traj;
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

        await defaultStructure(ctx, { trajectory });
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
            objectForEach(G3dSymbols, s => DefaultQueryRuntimeTable.addSymbol(s));
            this.ctx.managers.lociLabels.addProvider(G3dLabelProvider);
        }
        unregister() {
            this.ctx.state.data.actions.remove(LoadG3D);
            objectForEach(G3dSymbols, s => DefaultQueryRuntimeTable.removeSymbol(s));
            this.ctx.managers.lociLabels.removeProvider(G3dLabelProvider);
        }
    }
});