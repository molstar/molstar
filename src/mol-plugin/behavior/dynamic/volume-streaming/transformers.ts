/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject as SO, PluginStateTransform } from '../../../state/objects';
import { VolumeServerInfo, VolumeServerHeader } from './model';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Task } from 'mol-task';
import { PluginContext } from 'mol-plugin/context';
import { urlCombine } from 'mol-util/url';
import { createIsoValueParam } from 'mol-repr/volume/isosurface';
import { VolumeIsoValue } from 'mol-model/volume';
import { StateAction, StateObject, StateTransformer } from 'mol-state';
import { getStreamingMethod, getEmdbIdAndContourLevel } from './util';
import { VolumeStreaming } from './behavior';
import { VolumeRepresentation3DHelpers } from 'mol-plugin/state/transforms/representation';
import { BuiltInVolumeRepresentations } from 'mol-repr/volume/registry';
import { createTheme } from 'mol-theme/theme';
// import { PluginContext } from 'mol-plugin/context';

export const InitVolumeStreaming = StateAction.build({
    display: { name: 'Volume Streaming' },
    from: SO.Molecule.Structure,
    params(a) {
        return {
            method: PD.Select<VolumeServerInfo.Kind>(getStreamingMethod(a && a.data), [['em', 'EM'], ['x-ray', 'X-Ray']]),
            id: PD.Text((a && a.data.models[0].label) || ''),
            serverUrl: PD.Text('https://webchem.ncbr.muni.cz/DensityServer')
        };
    }
})(({ ref, state, params }, plugin: PluginContext) => Task.create('Volume Streaming', async taskCtx => {
    // TODO: custom react view for this and the VolumeStreamingBehavior transformer

    let dataId = params.id, emDefaultContourLevel: number | undefined;
    if (params.method === 'em') {
        await taskCtx.update('Getting EMDB info...');
        const emInfo = await getEmdbIdAndContourLevel(plugin, taskCtx, params.id);
        dataId = emInfo.emdbId;
        emDefaultContourLevel = emInfo.contour;
    }

    const infoTree = state.build().to(ref)
        .apply(CreateVolumeStreamingInfo, {
            serverUrl: params.serverUrl,
            source: params.method === 'em'
                ? { name: 'em', params: { isoValue: VolumeIsoValue.absolute(emDefaultContourLevel || 0) } }
                : { name: 'x-ray', params: { } },
            dataId
        });

    const infoObj = await state.updateTree(infoTree).runInContext(taskCtx);

    const behTree = state.build().to(infoTree.ref).apply(CreateVolumeStreamingBehavior, PD.getDefaultValues(VolumeStreaming.createParams(infoObj.data)))
    if (params.method === 'em') {
        behTree.apply(VolumeStreamingVisual, { channel: 'em' });
    } else {
        behTree.apply(VolumeStreamingVisual, { channel: '2fo-fc' });
        behTree.apply(VolumeStreamingVisual, { channel: 'fo-fc(+ve)' });
        behTree.apply(VolumeStreamingVisual, { channel: 'fo-fc(-ve)' });
    }
    await state.updateTree(behTree).runInContext(taskCtx);
}));

export { CreateVolumeStreamingInfo }
type CreateVolumeStreamingInfo = typeof CreateVolumeStreamingInfo
const CreateVolumeStreamingInfo = PluginStateTransform.BuiltIn({
    name: 'create-volume-streaming-info',
    display: { name: 'Volume Streaming Info' },
    from: SO.Molecule.Structure,
    to: VolumeServerInfo,
    params(a) {
        return {
            serverUrl: PD.Text('https://webchem.ncbr.muni.cz/DensityServer'),
            source: PD.MappedStatic('x-ray', {
                'em': PD.Group({
                    isoValue: createIsoValueParam(VolumeIsoValue.relative(1))
                }),
                'x-ray': PD.Group({ })
            }),
            dataId: PD.Text('')
        };
    }
})({
    apply: ({ a, params }, plugin: PluginContext) => Task.create('', async taskCtx => {
        const dataId = params.dataId;
        const emDefaultContourLevel = params.source.name === 'em' ? params.source.params.isoValue : VolumeIsoValue.relative(1);
        await taskCtx.update('Getting server header...');
        const header = await plugin.fetch<VolumeServerHeader>({ url: urlCombine(params.serverUrl, `${params.source.name}/${dataId}`), type: 'json' }).runInContext(taskCtx);
        const data: VolumeServerInfo.Data = {
            serverUrl: params.serverUrl,
            dataId,
            kind: params.source.name,
            header,
            emDefaultContourLevel,
            structure: a.data
        };
        return new VolumeServerInfo(data, { label: `Volume Streaming: ${dataId}` });
    })
});

export { CreateVolumeStreamingBehavior }
type CreateVolumeStreamingBehavior = typeof CreateVolumeStreamingBehavior
const CreateVolumeStreamingBehavior = PluginStateTransform.BuiltIn({
    name: 'create-volume-streaming-behavior',
    display: { name: 'Volume Streaming Behavior' },
    from: VolumeServerInfo,
    to: VolumeStreaming,
    params(a) {
        return VolumeStreaming.createParams(a && a.data);
    }
})({
    canAutoUpdate: ({ oldParams, newParams }) => {
        return oldParams.view === newParams.view;
    },
    apply: ({ a, params }, plugin: PluginContext) => Task.create('Volume streaming', async _ => {
        const behavior = new VolumeStreaming.Behavior(plugin, a.data);
        await behavior.update(params);
        return new VolumeStreaming(behavior, { label: 'Streaming Controls' });
    }),
    update({ b, newParams }) {
        return Task.create('Update Volume Streaming', async _ => {
            return await b.data.update(newParams) ? StateTransformer.UpdateResult.Updated : StateTransformer.UpdateResult.Unchanged;
        });
    }
});


export { VolumeStreamingVisual }
type VolumeStreamingVisual = typeof VolumeStreamingVisual
const VolumeStreamingVisual = PluginStateTransform.BuiltIn({
    name: 'create-volume-streaming-visual',
    display: { name: 'Volume Streaming Visual' },
    from: VolumeStreaming,
    to: SO.Volume.Representation3D,
    params: {
        channel: PD.Select<VolumeStreaming.ChannelType>('em', VolumeStreaming.ChannelTypeOptions, { isHidden: true })
    }
})({
    apply: ({ a, params: srcParams }, plugin: PluginContext) => Task.create('Volume Representation', async ctx => {
        const channel = a.data.channels[srcParams.channel];
        if (!channel) return StateObject.Null;

        const params = createVolumeProps(a.data, srcParams.channel);

        const provider = BuiltInVolumeRepresentations.isosurface;
        const props = params.type.params || {}
        const repr = provider.factory({ webgl: plugin.canvas3d.webgl, ...plugin.volumeRepresentation.themeCtx }, provider.getParams)
        repr.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: channel.data }, params))
        await repr.createOrUpdate(props, channel.data).runInContext(ctx);
        return new SO.Volume.Representation3D(repr, { label: `${Math.round(channel.isoValue.relativeValue * 100) / 100} σ [${srcParams.channel}]` });
    }),
    update: ({ a, b, oldParams, newParams }, plugin: PluginContext) => Task.create('Volume Representation', async ctx => {
        // TODO : check if params/underlying data/etc have changed; maybe will need to export "data" or some other "tag" in the Representation for this to work

        const channel = a.data.channels[newParams.channel];
        // TODO: is this correct behavior?
        if (!channel) return StateTransformer.UpdateResult.Unchanged;

        const params = createVolumeProps(a.data, newParams.channel);
        const props = { ...b.data.props, ...params.type.params };
        b.data.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: channel.data }, params))
        await b.data.createOrUpdate(props, channel.data).runInContext(ctx);
        return StateTransformer.UpdateResult.Updated;
    })
});

function createVolumeProps(streaming: VolumeStreaming.Behavior, channelName: VolumeStreaming.ChannelType) {
    const channel = streaming.channels[channelName]!;
    return VolumeRepresentation3DHelpers.getDefaultParamsStatic(streaming.plugin, 'isosurface',
        { isoValue: channel.isoValue, alpha: channel.opacity }, 'uniform', { value: channel.color });
}