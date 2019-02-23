/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF from 'mol-io/reader/cif';
import { Vec3 } from 'mol-math/linear-algebra';
import { volumeFromCcp4 } from 'mol-model-formats/volume/ccp4';
import { volumeFromDensityServerData } from 'mol-model-formats/volume/density-server';
import { volumeFromDsn6 } from 'mol-model-formats/volume/dsn6';
import { Task } from 'mol-task';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { VolumeStreaming } from 'mol-plugin/behavior/dynamic/volume';
import { PluginContext } from 'mol-plugin/context';
import { StateTransformer } from 'mol-state';
import { VolumeData, VolumeIsoValue } from 'mol-model/volume';
import { BuiltInVolumeRepresentations } from 'mol-repr/volume/registry';
import { createTheme } from 'mol-theme/theme';
import { VolumeRepresentation3DHelpers } from './representation';
import { Color } from 'mol-util/color';

export { VolumeFromCcp4 };
export { VolumeFromDsn6 };
export { VolumeFromDensityServerCif };
type VolumeFromCcp4 = typeof VolumeFromCcp4
const VolumeFromCcp4 = PluginStateTransform.BuiltIn({
    name: 'volume-from-ccp4',
    display: { name: 'Volume from CCP4/MRC/MAP', description: 'Create Volume from CCP4/MRC/MAP data' },
    from: SO.Format.Ccp4,
    to: SO.Volume.Data,
    params(a) {
        return {
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from CCP4/MRC/MAP', async ctx => {
            const volume = await volumeFromCcp4(a.data, params).runInContext(ctx)
            const props = { label: 'Volume' };
            return new SO.Volume.Data(volume, props);
        });
    }
});

type VolumeFromDsn6 = typeof VolumeFromDsn6
const VolumeFromDsn6 = PluginStateTransform.BuiltIn({
    name: 'volume-from-dsn6',
    display: { name: 'Volume from DSN6/BRIX', description: 'Create Volume from DSN6/BRIX data' },
    from: SO.Format.Dsn6,
    to: SO.Volume.Data,
    params(a) {
        return {
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from DSN6/BRIX', async ctx => {
            const volume = await volumeFromDsn6(a.data, params).runInContext(ctx)
            const props = { label: 'Volume' };
            return new SO.Volume.Data(volume, props);
        });
    }
});

type VolumeFromDensityServerCif = typeof VolumeFromDensityServerCif
const VolumeFromDensityServerCif = PluginStateTransform.BuiltIn({
    name: 'volume-from-density-server-cif',
    display: { name: 'Volume from density-server CIF', description: 'Identify and create all separate models in the specified CIF data block' },
    from: SO.Format.Cif,
    to: SO.Volume.Data,
    params(a) {
        if (!a) {
            return {
                blockHeader: PD.makeOptional(PD.Text(void 0, { description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.' }))
            };
        }
        const blocks = a.data.blocks.slice(1); // zero block contains query meta-data
        return {
            blockHeader: PD.makeOptional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' }))
        };
    }
})({
    isApplicable: a => a.data.blocks.length > 0,
    apply({ a, params }) {
        return Task.create('Parse density-server CIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[1].header; // zero block contains query meta-data
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const densityServerCif = CIF.schema.densityServer(block)
            const volume = await volumeFromDensityServerData(densityServerCif).runInContext(ctx)
            const props = { label: densityServerCif.volume_data_3d_info.name.value(0), description: `${densityServerCif.volume_data_3d_info.name.value(0)}` };
            return new SO.Volume.Data(volume, props);
        });
    }
});

export { VolumeStreamingBehavior }
type VolumeStreamingBehavior = typeof VolumeStreamingBehavior
const VolumeStreamingBehavior = PluginStateTransform.BuiltIn({
    name: 'volume-streaming-behavior',
    display: { name: 'Volume Streaming Behavior', description: 'Create Volume Streaming behavior.' },
    from: SO.Molecule.Structure,
    to: VolumeStreaming.Obj,
    params: VolumeStreaming.Params
})({
    canAutoUpdate: ({ oldParams, newParams }) => oldParams.serverUrl === newParams.serverUrl && oldParams.id === newParams.id,
    apply: ({ a, params }, plugin: PluginContext) => Task.create('Volume Streaming', async ctx => {
        const behavior = new VolumeStreaming.Behavior(plugin, params, a.data);
        // get the initial data now so that the child projections dont get empty volumes.
        await behavior.update(behavior.params);
        return new VolumeStreaming.Obj(behavior, { label: 'Volume Streaming' });
    }),
    update({ b, newParams }) {
        return Task.create('Update Volume Streaming', async _ => {
            await b.data.update(newParams);
            return StateTransformer.UpdateResult.Updated;
        });
    }
});

// export { VolumeStreamingData }
// type VolumeStreamingData = typeof VolumeStreamingData
// const VolumeStreamingData = PluginStateTransform.BuiltIn({
//     name: 'volume-streaming-data',
//     display: { name: 'Volume Streaming Data' },
//     from: VolumeStreaming.Obj,
//     to: SO.Volume.Data,
//     params: {
//         channel: PD.Select<keyof VolumeStreaming.ChannelData>('EM', [['EM', 'EM'], ['FO-FC', 'Fo-Fc'], ['2FO-FC', '2Fo-Fc']], { isHidden: true }),
//         level: PD.Text<VolumeStreaming.LevelType>('em')
//     }
// })({
//     apply({ a, params }, plugin: PluginContext) {
//         const data = a.data.currentData[params.channel] || VolumeData.Empty;
//         console.log({ data });
//         return new SO.Volume.Data(a.data.currentData[params.channel] || VolumeData.Empty, { label: params.level });
//     }
// });

export { VolumeStreamingVisual }
type VolumeStreamingVisual = typeof VolumeStreamingVisual
const VolumeStreamingVisual = PluginStateTransform.BuiltIn({
    name: 'volume-streaming-visual',
    display: { name: 'Volume Streaming Visual' },
    from: VolumeStreaming.Obj,
    to: SO.Volume.Representation3D,
    params: {
        channel: PD.Select<keyof VolumeStreaming.ChannelData>('EM', [['EM', 'EM'], ['FO-FC', 'Fo-Fc'], ['2FO-FC', '2Fo-Fc']], { isHidden: true }),
        level: PD.Text<VolumeStreaming.LevelType>('em')
    }
})({
    apply: ({ a, params: srcParams }, plugin: PluginContext) => Task.create('Volume Representation', async ctx => {
        const { data, params } = createVolumeProps(a.data, srcParams.channel, srcParams.level);

        const provider = BuiltInVolumeRepresentations.isosurface;
        const props = params.type.params || {}
        const repr = provider.factory({ webgl: plugin.canvas3d.webgl, ...plugin.volumeRepresentation.themeCtx }, provider.getParams)
        repr.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: data }, params))
        await repr.createOrUpdate(props, data).runInContext(ctx);
        return new SO.Volume.Representation3D(repr, { label: srcParams.level, description: VolumeRepresentation3DHelpers.getDescription(props) });
    }),
    update: ({ a, b, oldParams, newParams }, plugin: PluginContext) => Task.create('Volume Representation', async ctx => {
        // TODO : check if params/underlying data/etc have changed; maybe will need to export "data" or some other "tag" in the Representation for this to work
        const { data, params } = createVolumeProps(a.data, newParams.channel, newParams.level);
        const props = { ...b.data.props, ...params.type.params };
        b.data.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: data }, params))
        await b.data.createOrUpdate(props, data).runInContext(ctx);
        return StateTransformer.UpdateResult.Updated;
    })
});

function createVolumeProps(streaming: VolumeStreaming.Behavior, channel: keyof VolumeStreaming.ChannelData, level: VolumeStreaming.LevelType) {
    const data = streaming.currentData[channel] || VolumeData.One;
    // TODO: createTheme fails when VolumeData.Empty is used for some reason.

    let isoValue: VolumeIsoValue, color: Color;

    if (level === 'em' && streaming.params.levels.name === 'em') {
        isoValue = streaming.params.levels.params.isoValue;
        color = streaming.params.levels.params.color;
    } else if (level !== 'em' && streaming.params.levels.name === 'x-ray') {
        isoValue = streaming.params.levels.params[level].isoValue;
        color = streaming.params.levels.params[level].color;
    } else {
        throw new Error(`Unsupported iso level ${level}.`);
    }

    const params = VolumeRepresentation3DHelpers.getDefaultParamsStatic(streaming.ctx, 'isosurface', { isoValue, alpha: 0.3 }, 'uniform', { value: color });
    return { data, params };
}