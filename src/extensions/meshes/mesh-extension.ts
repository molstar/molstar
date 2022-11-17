/** Defines new types of State tree transformers for dealing with mesh data. */

import * as MS from './molstar-lib-imports';
import PD = MS.ParamDefinition;

import * as MeshUtils from './mesh-utils';


export const CellStarTransform: MS.StateTransformer.Builder.Root = MS.StateTransformer.builderFactory('cellstar');


// // // // // // // // // // // // // // // // // // // // // // // //
// Parsed data

/** Data type for `MeshlistStateObject` - list of meshes */
export interface MeshlistData {
    segmentId: number,
    segmentName: string,
    detail: number,
    meshIds: number[],
    mesh: MS.Mesh,
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
            mesh: MS.Mesh.createEmpty(),
        };
    };
    export async function fromCIF(data: MS.CifFile, segmentId: number, segmentName: string, detail: number): Promise<MeshlistData> {
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
    export function getShape(data: MeshlistData, color: MS.Color): MS.Shape<MS.Mesh> {
        const mesh = data.mesh;
        const meshShape: MS.Shape<MS.Mesh> = MS.Shape.create('MyShape', data, mesh,
            () => color,
            () => 1, (group) => `${data.segmentName} | Segment ${data.segmentId} | Detail ${data.detail} | Mesh ${group}`);
        return meshShape;
    }

    export function combineBBoxes(boxes: (MS.Box3D | null)[]): MS.Box3D | null {
        let result = null;
        for (const box of boxes) {
            if (!box) continue;
            if (result) {
                MS.Vec3.min(result.min, result.min, box.min);
                MS.Vec3.max(result.max, result.max, box.max);
            } else {
                result = MS.Box3D.zero();
                MS.Box3D.copy(result, box);
            }
        }
        return result;
    }
    export function bbox(data: MeshlistData): MS.Box3D | null {
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

export class MeshlistStateObject extends MS.PluginStateObject.Create<MeshlistData>({ name: 'Parsed Meshlist', typeClass: 'Object' }) { }
// QUESTION: is typeClass just for color, or does do something?

export const ParseMeshlistTransformer = CellStarTransform({
    name: 'meshlist-from-string',
    from: MS.PluginStateObject.Format.Cif,
    to: MeshlistStateObject,
    params: {
        label: PD.Text(MeshlistStateObject.type.name, { isHidden: true }), // QUESTION: Is this the right way to pass a value to apply() without exposing it in GUI?
        segmentId: PD.Numeric(1, {}, { isHidden: true }),
        segmentName: PD.Text('Segment'),
        detail: PD.Numeric(1, {}, { isHidden: true }),
    }
})({
    apply({ a, params }, globalCtx) { // `a` is the parent node, params are 2nd argument to To.apply(), `globalCtx` is the plugin
        return MS.Task.create('Create Parsed Meshlist', async ctx => {
            const meshlistData = await MeshlistData.fromCIF(a.data, params.segmentId, params.segmentName, params.detail);
            const es = meshlistData.meshIds.length === 1 ? '' : 'es';
            return new MeshlistStateObject(meshlistData, { label: params.label, description: `${meshlistData.segmentName} (${meshlistData.meshIds.length} mesh${es})` });
        });
    }
});


// // // // // // // // // // // // // // // // // // // // // // // //
// Parsed data -> Shape

/** Data type for PluginStateObject.Shape.Provider */
type MeshShapeProvider = MS.ShapeProvider<MeshlistData, MS.Mesh, MS.Mesh.Params>;
namespace MeshShapeProvider {
    export function fromMeshlistData(meshlist: MeshlistData, color?: MS.Color): MeshShapeProvider {
        const theColor = color ?? MeshUtils.ColorGenerator.next().value;
        return {
            label: 'Mesh',
            data: meshlist,
            params: meshParamDef, // TODO how to pass the real params correctly?
            geometryUtils: MS.Mesh.Utils,
            getShape: (ctx, data: MeshlistData) => MeshlistData.getShape(data, theColor),
        };
    }
}

/** Params for MeshShapeTransformer */
const meshShapeParamDef = {
    color: PD.Value<MS.Color | undefined>(undefined), // undefined means random color
};

const meshParamDef: MS.Mesh.Params = {
    // These are basically original MS.Mesh.Params:
    // BaseGeometry.Params
    alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity', isEssential: true, description: 'How opaque/transparent the representation is rendered.' }),
    quality: PD.Select<MS.VisualQuality>('auto', MS.VisualQualityOptions, { isEssential: true, description: 'Visual/rendering quality of the representation.' }),
    material: MS.Material.getParam(),
    clip: MS.Mesh.Params.clip, // PD.Group(MS.Clip.Params),
    instanceGranularity: PD.Boolean(false, { description: 'Use instance granularity for marker, transparency, clipping, overpaint, substance data to save memory.' }),
    // Mesh.Params
    doubleSided: PD.Boolean(false, MS.BaseGeometry.CustomQualityParamInfo),
    flipSided: PD.Boolean(false, MS.BaseGeometry.ShadingCategory),
    flatShaded: PD.Boolean(true, MS.BaseGeometry.ShadingCategory), // CHANGED, default: false (set true to see the real mesh vertices and triangles)
    ignoreLight: PD.Boolean(false, MS.BaseGeometry.ShadingCategory),
    xrayShaded: PD.Boolean(false, MS.BaseGeometry.ShadingCategory), // this is like better opacity (angle-dependent), nice
    transparentBackfaces: PD.Select('off', PD.arrayToOptions(['off', 'on', 'opaque']), MS.BaseGeometry.ShadingCategory),
    bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, MS.BaseGeometry.ShadingCategory),
    bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, MS.BaseGeometry.ShadingCategory),
    // TODO when I change values here, it has effect, but not if I change them in GUI
};

export const MeshShapeTransformer = CellStarTransform({
    name: 'shape-from-meshlist',
    display: { name: 'Shape from Meshlist', description: 'Create Shape from Meshlist data' },
    from: MeshlistStateObject,
    to: MS.PluginStateObject.Shape.Provider,
    params: meshShapeParamDef
})({
    apply({ a, params }) {
        // you can look for example at ShapeFromPly in mol-plugin-state/tansforms/model.ts as an example
        const shapeProvider = MeshShapeProvider.fromMeshlistData(a.data, params.color);
        return new MS.PluginStateObject.Shape.Provider(shapeProvider, { label: MS.PluginStateObject.Shape.Provider.type.name, description: a.description });
    }
});


// // // // // // // // // // // // // // // // // // // // // // // //
// Shape -> Repr

// type MeshRepr = MS.PluginStateObject.Representation3DData<MS.ShapeRepresentation<MS.ShapeProvider<any,any,any>, MS.Mesh, MS.Mesh.Params>, any>;

// export const CustomMeshReprTransformer = CellStarTransform({
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
