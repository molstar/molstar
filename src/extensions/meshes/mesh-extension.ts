/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

/** Defines new types of State tree transformers for dealing with mesh data. */


import { BaseGeometry, VisualQuality, VisualQualityOptions } from '../../mol-geo/geometry/base';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { CifFile } from '../../mol-io/reader/cif';
import { Box3D } from '../../mol-math/geometry';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Shape } from '../../mol-model/shape';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { Color } from '../../mol-util/color';
import { Material } from '../../mol-util/material';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import * as MeshUtils from './mesh-utils';


export const VolsegTransform: StateTransformer.Builder.Root = StateTransformer.builderFactory('volseg');


// // // // // // // // // // // // // // // // // // // // // // // //
// Parsed data

/** Data type for `MeshlistStateObject` - list of meshes */
export interface MeshlistData {
    segmentId: number,
    segmentName: string,
    detail: number,
    meshIds: number[],
    mesh: Mesh,
    /** Reference to the object which created this meshlist (e.g. `MeshStreaming.Behavior`) */
    ownerId?: string,
}

export namespace MeshlistData {
    export function empty(): MeshlistData {
        return {
            segmentId: 0,
            segmentName: 'Empty',
            detail: 0,
            meshIds: [],
            mesh: Mesh.createEmpty(),
        };
    };
    export async function fromCIF(data: CifFile, segmentId: number, segmentName: string, detail: number): Promise<MeshlistData> {
        const { mesh, meshIds } = await MeshUtils.meshFromCif(data);
        return {
            segmentId,
            segmentName,
            detail,
            meshIds,
            mesh,
        };
    }
    export function stats(meshListData: MeshlistData): string {
        return `Meshlist "${meshListData.segmentName}" (detail ${meshListData.detail}): ${meshListData.meshIds.length} meshes, ${meshListData.mesh.vertexCount} vertices, ${meshListData.mesh.triangleCount} triangles`;
    }
    export function getShape(data: MeshlistData, color: Color): Shape<Mesh> {
        const mesh = data.mesh;
        const meshShape: Shape<Mesh> = Shape.create(data.segmentName, data, mesh,
            () => color,
            () => 1,
            // group => `${data.segmentName} | Segment ${data.segmentId} | Detail ${data.detail} | Mesh ${group}`,
            group => data.segmentName,
        );
        return meshShape;
    }

    export function combineBBoxes(boxes: (Box3D | null)[]): Box3D | null {
        let result = null;
        for (const box of boxes) {
            if (!box) continue;
            if (result) {
                Vec3.min(result.min, result.min, box.min);
                Vec3.max(result.max, result.max, box.max);
            } else {
                result = Box3D.zero();
                Box3D.copy(result, box);
            }
        }
        return result;
    }
    export function bbox(data: MeshlistData): Box3D | null {
        return MeshUtils.bbox(data.mesh);
    }

    export function allVerticesUsed(data: MeshlistData): boolean {
        const unusedVertices = new Set();
        for (let i = 0; i < data.mesh.vertexCount; i++) {
            unusedVertices.add(i);
        }
        for (let i = 0; i < 3 * data.mesh.triangleCount; i++) {
            const v = data.mesh.vertexBuffer.ref.value[i];
            unusedVertices.delete(v);
        }
        return unusedVertices.size === 0;
    }
}



// // // // // // // // // // // // // // // // // // // // // // // //
// Raw Data -> Parsed data

export class MeshlistStateObject extends PluginStateObject.Create<MeshlistData>({ name: 'Parsed Meshlist', typeClass: 'Object' }) { }
// QUESTION: is typeClass just for color, or does do something?

export const ParseMeshlistTransformer = VolsegTransform({
    name: 'meshlist-from-string',
    from: PluginStateObject.Format.Cif,
    to: MeshlistStateObject,
    params: {
        label: PD.Text(MeshlistStateObject.type.name, { isHidden: true }), // QUESTION: Is this the right way to pass a value to apply() without exposing it in GUI?
        segmentId: PD.Numeric(1, {}, { isHidden: true }),
        segmentName: PD.Text('Segment'),
        detail: PD.Numeric(1, {}, { isHidden: true }),
        /** Reference to the object which manages this meshlist (e.g. `MeshStreaming.Behavior`) */
        ownerId: PD.Text('', { isHidden: true }),
    }
})({
    apply({ a, params }, globalCtx) { // `a` is the parent node, params are 2nd argument to To.apply(), `globalCtx` is the plugin
        return Task.create('Create Parsed Meshlist', async ctx => {
            const meshlistData = await MeshlistData.fromCIF(a.data, params.segmentId, params.segmentName, params.detail);
            meshlistData.ownerId = params.ownerId;
            const es = meshlistData.meshIds.length === 1 ? '' : 'es';
            return new MeshlistStateObject(meshlistData, { label: params.label, description: `${meshlistData.segmentName} (${meshlistData.meshIds.length} mesh${es})` });
        });
    }
});


// // // // // // // // // // // // // // // // // // // // // // // //
// Parsed data -> Shape

/** Data type for PluginStateObject.Shape.Provider */
type MeshShapeProvider = ShapeProvider<MeshlistData, Mesh, Mesh.Params>;
namespace MeshShapeProvider {
    export function fromMeshlistData(meshlist: MeshlistData, color?: Color): MeshShapeProvider {
        const theColor = color ?? MeshUtils.ColorGenerator.next().value;
        return {
            label: 'Mesh',
            data: meshlist,
            params: meshParamDef, // TODO how to pass the real params correctly?
            geometryUtils: Mesh.Utils,
            getShape: (ctx, data: MeshlistData) => MeshlistData.getShape(data, theColor),
        };
    }
}

/** Params for MeshShapeTransformer */
const meshShapeParamDef = {
    color: PD.Value<Color | undefined>(undefined), // undefined means random color
};

const meshParamDef: Mesh.Params = {
    ...Mesh.Params,
    // These are basically original MS.Mesh.Params:
    // BaseGeometry.Params
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity', isEssential: true, description: 'How opaque/transparent the representation is rendered.' }),
    quality: PD.Select<VisualQuality>('custom', VisualQualityOptions, { isEssential: true, description: 'Visual/rendering quality of the representation.' }), // use 'custom' when wanting to apply doubleSided
    material: Material.getParam(),
    clip: Mesh.Params.clip, // PD.Group(MS.Clip.Params),
    instanceGranularity: PD.Boolean(false, { description: 'Use instance granularity for marker, transparency, clipping, overpaint, substance data to save memory.' }),
    // Mesh.Params
    doubleSided: PD.Boolean(true, BaseGeometry.CustomQualityParamInfo), // default: false (set true, to show at least something in weird cases)
    flipSided: PD.Boolean(false, BaseGeometry.ShadingCategory),
    flatShaded: PD.Boolean(false, BaseGeometry.ShadingCategory), // default: false (set true to see the real mesh vertices and triangles)
    ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
    xrayShaded: PD.Boolean(false, BaseGeometry.ShadingCategory), // this is like better opacity (angle-dependent), nice
    transparentBackfaces: PD.Select('off', PD.arrayToOptions(['off', 'on', 'opaque']), BaseGeometry.ShadingCategory),
    bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, BaseGeometry.ShadingCategory),
    // TODO when I change values here, it has effect, but not if I change them in GUI
};

export const MeshShapeTransformer = VolsegTransform({
    name: 'shape-from-meshlist',
    display: { name: 'Shape from Meshlist', description: 'Create Shape from Meshlist data' },
    from: MeshlistStateObject,
    to: PluginStateObject.Shape.Provider,
    params: meshShapeParamDef
})({
    apply({ a, params }) {
        // you can look for example at ShapeFromPly in mol-plugin-state/tansforms/model.ts as an example
        const shapeProvider = MeshShapeProvider.fromMeshlistData(a.data, params.color);
        return new PluginStateObject.Shape.Provider(shapeProvider, { label: PluginStateObject.Shape.Provider.type.name, description: a.description });
    }
});


// // // // // // // // // // // // // // // // // // // // // // // //
// Shape -> Repr

// type MeshRepr = MS.PluginStateObject.Representation3DData<MS.ShapeRepresentation<MS.ShapeProvider<any,any,any>, MS.Mesh, MS.Mesh.Params>, any>;

// export const CustomMeshReprTransformer = VolsegTransform({
//     name: 'custom-repr',
//     from: MS.PluginStateObject.Shape.Provider, // later we can change this
//     to: MS.PluginStateObject.Shape.Representation3D,
// })({
//     apply({ a }, globalCtx) {
//         const repr: MeshRepr = createRepr(a.data); // TODO implement createRepr
//         // have a look at MS.StateTransforms.Representation.ShapeRepresentation3D if you want to try implementing yourself
//         return new MS.PluginStateObject.Shape.Representation3D(repr)
//     },
// })

// export async function createMeshRepr(plugin: MS.PluginContext, data: any) {
//     await plugin.build()
//         .toRoot()
//         .apply(CreateMyShapeTransformer, { data })
//         .apply(MS.StateTransforms.Representation.ShapeRepresentation3D) // this should work
//         // or .apply(CustomMeshRepr)
//         .commit();
// }

// export function createRepr(reprData: MS.ShapeProvider<any,any,any>): MeshRepr {
//     throw new Error('NotImplemented');
//     return {} as MeshRepr;
// }
