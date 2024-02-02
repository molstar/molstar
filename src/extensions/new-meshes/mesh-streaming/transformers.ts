/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { ShapeRepresentation } from '../../../mol-repr/shape/representation';
import { StateAction, StateTransformer } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { shallowEqualObjects } from '../../../mol-util';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

import { BACKGROUND_OPACITY, FOREROUND_OPACITY, MeshlistData, VolsegTransform } from '../mesh-extension';
import { MeshStreaming, NO_SEGMENT } from './behavior';
import { MeshServerInfo } from './server-info';


// // // // // // // // // // // // // // // // // // // // // // // //

export const MeshServerTransformer = VolsegTransform({
    name: 'mesh-server-info',
    from: PluginStateObject.Root,
    to: MeshServerInfo,
    params: MeshServerInfo.Params,
})({
    apply({ a, params }, plugin: PluginContext) { // `a` is the parent node, `params` are 2nd argument to To.apply()
        params.serverUrl = params.serverUrl.replace(/\/*$/, ''); // trim trailing slash
        const description: string = params.entryId;
        return new MeshServerInfo({ ...params }, { label: 'Mesh Server', description: description });
    }
});

// // // // // // // // // // // // // // // // // // // // // // // //

export const MeshStreamingTransformer = VolsegTransform({
    name: 'mesh-streaming-from-server-info',
    display: { name: 'Mesh Streaming' },
    from: MeshServerInfo,
    to: MeshStreaming,
    params: a => MeshStreaming.Params.create(a!.data),
})({
    canAutoUpdate() { return true; },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Mesh Streaming', async ctx => {
            const behavior = new MeshStreaming.Behavior(plugin, a.data, params);
            await behavior.update(params);
            return new MeshStreaming(behavior, { label: 'Mesh Streaming', description: behavior.getDescription() });
        });
    },
    update({ a, b, oldParams, newParams }) {
        return Task.create('Update Mesh Streaming', async ctx => {
            if (a.data.source !== b.data.parentData.source || a.data.entryId !== b.data.parentData.entryId) {
                return StateTransformer.UpdateResult.Recreate;
            }
            b.data.parentData = a.data;
            await b.data.update(newParams);
            b.description = b.data.getDescription();
            return StateTransformer.UpdateResult.Updated;
        });
    }
});

// // // // // // // // // // // // // // // // // // // // // // // //

interface MeshVisualGroupData {
    opacity: number,
}

// export type MeshVisualGroupTransformer = typeof MeshVisualGroupTransformer;
export const MeshVisualGroupTransformer = VolsegTransform({
    name: 'mesh-visual-group-from-streaming',
    display: { name: 'Mesh Visuals for a Segment' },
    from: MeshStreaming,
    to: PluginStateObject.Group,
    params: {
        /** Shown on the node in GUI */
        label: PD.Text('', { isHidden: true }),
        /** Shown on the node in GUI (gray letters) */
        description: PD.Text(''),
        segmentId: PD.Numeric(NO_SEGMENT, {}, { isHidden: true }),
        opacity: PD.Numeric(-1, { min: 0, max: 1, step: 0.01 }),
    }
})({
    apply({ a, params }, plugin) {
        trySetAutoOpacity(params, a);
        return new PluginStateObject.Group({ opacity: params.opacity }, params);
    },
    update({ a, b, oldParams, newParams }, plugin) {
        if (shallowEqualObjects(oldParams, newParams)) {
            return StateTransformer.UpdateResult.Unchanged;
        }
        newParams.label ||= oldParams.label; // Protect against resetting params to invalid defaults
        if (newParams.segmentId === NO_SEGMENT) newParams.segmentId = oldParams.segmentId; // Protect against resetting params to invalid defaults
        trySetAutoOpacity(newParams, a);
        b.label = newParams.label;
        b.description = newParams.description;
        (b.data as MeshVisualGroupData).opacity = newParams.opacity;
        return StateTransformer.UpdateResult.Updated;
    },
    canAutoUpdate({ oldParams, newParams }, plugin) {
        return newParams.description === oldParams.description;
    },
});

function trySetAutoOpacity(params: StateTransformer.Params<typeof MeshVisualGroupTransformer>, parent: MeshStreaming) {
    if (params.opacity === -1) {
        const isBgSegment = parent.data.backgroundSegments[params.segmentId];
        if (isBgSegment !== undefined) {
            params.opacity = isBgSegment ? BACKGROUND_OPACITY : FOREROUND_OPACITY;
        }
    }
}


// // // // // // // // // // // // // // // // // // // // // // // //

export const MeshVisualTransformer = VolsegTransform({
    name: 'mesh-visual-from-streaming',
    display: { name: 'Mesh Visual from Streaming' },
    from: MeshStreaming,
    to: PluginStateObject.Shape.Representation3D,
    params: {
        /** Must be set to PluginStateObject reference to self */
        ref: PD.Text('', { isHidden: true, isEssential: true }), // QUESTION what is isEssential
        /** Identification of the mesh visual, e.g. 'low-2' */
        tag: PD.Text('', { isHidden: true, isEssential: true }),
        /** Opacity of the visual (not to be set directly, but controlled by the opacity of the parent Group, and by VisualInfo.visible) */
        opacity: PD.Numeric(-1, { min: 0, max: 1, step: 0.01 }, { isHidden: true }),
    }
})({
    apply({ a, params, spine }, plugin: PluginContext) {
        return Task.create('Mesh Visual', async ctx => {
            const visualInfo: MeshStreaming.VisualInfo = a.data.visuals![params.tag];
            if (!visualInfo) throw new Error(`VisualInfo with tag '${params.tag}' is missing.`);
            const groupData = spine.getAncestorOfType(PluginStateObject.Group)?.data as MeshVisualGroupData | undefined;
            params.opacity = visualInfo.visible ? (groupData?.opacity ?? FOREROUND_OPACITY) : 0.0;
            const props = PD.getDefaultValues(Mesh.Params);
            props.flatShaded = true; // `flatShaded: true` is to see the real mesh vertices and triangles (default: false)
            props.alpha = params.opacity;
            const repr = ShapeRepresentation((ctx, meshlist: MeshlistData) => MeshlistData.getShape(meshlist, visualInfo.color), Mesh.Utils);
            await repr.createOrUpdate(props, visualInfo.data ?? MeshlistData.empty()).runInContext(ctx);
            return new PluginStateObject.Shape.Representation3D({ repr, sourceData: visualInfo.data }, { label: 'Mesh Visual', description: params.tag });
        });
    },
    update({ a, b, oldParams, newParams, spine }, plugin: PluginContext) {
        return Task.create('Update Mesh Visual', async ctx => {
            newParams.ref ||= oldParams.ref; // Protect against resetting params to invalid defaults
            newParams.tag ||= oldParams.tag; // Protect against resetting params to invalid defaults
            const visualInfo: MeshStreaming.VisualInfo = a.data.visuals![newParams.tag];
            if (!visualInfo) throw new Error(`VisualInfo with tag '${newParams.tag}' is missing.`);
            const oldData = b.data.sourceData as MeshlistData | undefined;
            if (visualInfo.data?.detail !== oldData?.detail) {
                return StateTransformer.UpdateResult.Recreate;
            }
            const groupData = spine.getAncestorOfType(PluginStateObject.Group)?.data as MeshVisualGroupData | undefined;
            const newOpacity = visualInfo.visible ? (groupData?.opacity ?? FOREROUND_OPACITY) : 0.0; // do not store to newParams directly, because oldParams and newParams might point to the same object!
            if (newOpacity !== oldParams.opacity) {
                newParams.opacity = newOpacity;
                await b.data.repr.createOrUpdate({ alpha: newParams.opacity }).runInContext(ctx);
                return StateTransformer.UpdateResult.Updated;
            } else {
                return StateTransformer.UpdateResult.Unchanged;
            }
        });
    },
    canAutoUpdate(params, globalCtx) {
        return true;
    },
    dispose({ b, params }, plugin) {
        b?.data.repr.destroy(); // QUESTION is this correct?
    },
});

// // // // // // // // // // // // // // // // // // // // // // // //

export const InitMeshStreaming = StateAction.build({
    display: { name: 'Mesh Streaming' },
    from: PluginStateObject.Root,
    params: MeshServerInfo.Params,
    isApplicable: (a, _, plugin: PluginContext) => true
})(function (p, plugin: PluginContext) {
    return Task.create('Mesh Streaming', async ctx => {
        const { params } = p;
        // p.ref
        const serverNode = await plugin.build().to(p.ref).apply(MeshServerTransformer, params).commit();
        // const serverNode = await plugin.build().toRoot().apply(MeshServerTransformer, params).commit();
        const streamingNode = await plugin.build().to(serverNode).apply(MeshStreamingTransformer, {}).commit();
        const visuals = streamingNode.data?.visuals ?? {};
        const bgSegments = streamingNode.data?.backgroundSegments ?? {};

        const segmentGroups: { [segid: number]: string } = {};
        for (const tag in visuals) {
            const segid = visuals[tag].segmentId;
            if (!segmentGroups[segid]) {
                let description = visuals[tag].segmentName;
                if (bgSegments[segid]) description += ' (background)';
                const group = await plugin.build().to(streamingNode).apply(MeshVisualGroupTransformer, { label: `Segment ${segid}`, description: description, segmentId: segid }, { state: { isCollapsed: true } }).commit();
                segmentGroups[segid] = group.ref;
            }
        }
        const visualsUpdate = plugin.build();
        for (const tag in visuals) {
            const ref = `${streamingNode.ref}-${tag}`;
            const segid = visuals[tag].segmentId;
            visualsUpdate.to(segmentGroups[segid]).apply(MeshVisualTransformer, { ref: ref, tag: tag }, { ref: ref }); // ref - hack to allow the node make itself invisible
        }
        await plugin.state.data.updateTree(visualsUpdate).runInContext(ctx); // QUESTION what is really the difference between this and `visualsUpdate.commit()`?
    });
});

// TODO make available in GUI, in left panel or in right panel like Volume Streaming in src/mol-plugin-ui/structure/volume.tsx?
