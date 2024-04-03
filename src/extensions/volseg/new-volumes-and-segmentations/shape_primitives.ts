// import { StateTransformer } from '../../mol-state';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Shape } from '../../../mol-model/shape';
import { Color } from '../../../mol-util/color';
import { StateTransformer } from '../../../mol-state';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Box } from '../../../mol-geo/primitive/box';
import { Pyramid, TriangularPyramid } from '../../../mol-geo/primitive/pyramid';
import { Primitive } from '../../../mol-geo/primitive/primitive';
import { polygon } from '../../../mol-geo/primitive/polygon';
import { addEllipsoid } from '../../../mol-geo/geometry/mesh/builder/ellipsoid';
import { DescriptionData, SegmentAnnotationData, Cylinder, ShapePrimitiveData, Ellipsoid, PyramidPrimitive, Sphere, BoxPrimitive, Vector4 } from './volseg-api/data';



export class VolsegShapePrimitivesData {
    constructor(public shapePrimitiveData: ShapePrimitiveData) {
    }
}

// 'Data', try 'Object' later
export class VolsegGeometricSegmentation extends PluginStateObject.Create<VolsegShapePrimitivesData>({ name: 'Vol & Seg Geometric Segmentation', typeClass: 'Data' }) { }


let cone: Primitive;
export function Cone() {
    if (!cone) cone = Pyramid(polygon(48, true));
    return cone;
}

function rgbaToHex(rgbaNormalized: Vector4) {
    const rgba = rgbaNormalized.map(i => Math.round(i * 255));
    const [red, green, blue, opacity] = rgba;

    let r = Math.round(red).toString(16);
    let g = Math.round(green).toString(16);
    let b = Math.round(blue).toString(16);
    let a = Math.round(opacity).toString(16);

    if (r.length === 1)
        r = '0' + r;
    if (g.length === 1)
        g = '0' + g;
    if (b.length === 1)
        b = '0' + b;
    if (a.length === 1)
        a = '0' + a;

    const hexString = '#' + r + g + b + a;
    // const hexNumber = parseInt(hexString.replace(/^#/, ''), 16);
    // console.log(`Hex number ${hexNumber} for hex string ${hexString} and RGBA ${red}, ${green}, ${blue}, ${opacity}`);
    return hexString;
}



// export type ShapePrimitive =
//     | { kind: 'sphere', center: number[], radius: number, label: string, color: number }
//     | { kind: 'cylinder', start: number[], end: number[], radius: number, label: string, color: number }
//     | { kind: 'box', translation: number[], scaling: number[], label: string, color: number }
//     | { kind: 'pyramid', translation: number[], scaling: number[], label: string, color: number }
//     | { kind: 'cone', translation: number[], scaling: number[], label: string, color: number }
//     | { kind: 'ellipsoid', dir_major: number[], dir_minor: number[], center: number[], radius_scale: number[], label: string, color: number }

// export type ShapePrimitivesData = ShapePrimitive[]

function addBox(state: MeshBuilder.State,
    translation: Vec3 = Vec3.create(0.5, 0.5, 0.5),
    scaling: Vec3 = Vec3.create(1, 1, 1)) {
    const mat4 = Mat4.identity();
    Mat4.scale(mat4, mat4, scaling);
    Mat4.translate(mat4, mat4, translation);

    MeshBuilder.addPrimitive(state, mat4, Box());
}

function addTriangularPyramid(state: MeshBuilder.State,
    translation: Vec3 = Vec3.create(0.5, 0.5, 0.5),
    scaling: Vec3 = Vec3.create(1, 1, 1)) {
    const mat4 = Mat4.identity();
    Mat4.scale(mat4, mat4, scaling);
    Mat4.translate(mat4, mat4, translation);
    MeshBuilder.addPrimitive(state, mat4, TriangularPyramid());
}

// function addCone(state: MeshBuilder.State,
//     translation: Vec3 = [0.5, 0.5, 0.5] as Vec3,
//     scaling: Vec3 = [1, 1, 1] as Vec3) {
//     const mat4 = Mat4.identity();
//     Mat4.scale(mat4, mat4, scaling);
//     Mat4.translate(mat4, mat4, translation);
//     MeshBuilder.addPrimitive(state, mat4, Cone());
// }

export const isShapePrimitiveParamsValues = (value: CreateShapePrimitiveProviderParamsValues): value is CreateShapePrimitiveProviderParamsValues => !!value?.segmentAnnotations;

export type CreateShapePrimitiveProviderParamsValues = PD.Values<typeof CreateShapePrimitiveProviderParams>;
export const CreateShapePrimitiveProviderParams = {
    // data: PD.Value<ShapePrimitiveData>([] as any, { isHidden: true }),
    // data: PD.Value<BoxPrimitive | Sphere | Cylinder | Ellipsoid | PyramidPrimitive>([] as any, { isHidden: true }),
    segmentAnnotations: PD.Value<SegmentAnnotationData[]>([] as any, { isHidden: true }),
    descriptions: PD.Value<DescriptionData[]>([] as any, { isHidden: true }),
    segmentationId: PD.Text(''),
    segmentId: PD.Numeric(0),
};

const Transform = StateTransformer.builderFactory('msvolseg');
export const CreateShapePrimitiveProvider = Transform({
    name: 'create-shape-primitive-provider',
    display: { name: 'Shape Primitives' },
    from: VolsegGeometricSegmentation,
    to: PluginStateObject.Shape.Provider,
    params: CreateShapePrimitiveProviderParams
})({
    apply({ a, params }) {
        return new PluginStateObject.Shape.Provider({
            label: 'Shape Primitives',
            data: params,
            // data: params.data,
            params: Mesh.Params,
            geometryUtils: Mesh.Utils,
            getShape: (_, data) => createShapePrimitive(a.data.shapePrimitiveData, params)
        }, { label: 'Shape Primitives' });
    }
});

export const CreateShapePrimitiveProviderCVSX = Transform({
    name: 'create-shape-primitive-provider-cvsx',
    display: { name: 'Shape Primitives' },
    from: PluginStateObject.Data.String,
    to: PluginStateObject.Shape.Provider,
    params: CreateShapePrimitiveProviderParams
})({
    apply({ a, params }) {
        const shapePrimitiveData: ShapePrimitiveData = JSON.parse(a.data);
        return new PluginStateObject.Shape.Provider({
            label: 'Shape Primitives',
            data: params,
            // data: params.data,
            params: Mesh.Params,
            geometryUtils: Mesh.Utils,
            getShape: (_, data) => createShapePrimitive(shapePrimitiveData, params)
        }, { label: 'Shape Primitives' });
    }
});


function _get_target_segment_name(allDescriptions: DescriptionData[], segment_id: number) {
    // NOTE: for now single description
    const description = allDescriptions.filter(d => d.target_id && d.target_id.segment_id === segment_id);
    console.log(`Target segment name is ${description[0].name!}`);
    return description[0].name!;
}

function _get_target_segment_color_as_hex(allSegmentAnnotations: SegmentAnnotationData[], segment_id: number) {
    // NOTE: for now single annotation, should be single one
    const annotation = allSegmentAnnotations.filter(a => a.segment_id === segment_id);
    const colorAsArray = annotation[0].color!;
    const colorAsHex = rgbaToHex(colorAsArray);
    return colorAsHex;
}

function createShapePrimitive(data: ShapePrimitiveData, params: CreateShapePrimitiveProviderParamsValues) {
    const builder = MeshBuilder.createState(512, 512);
    const descriptions = params.descriptions;
    const segmentAnnotations = params.segmentAnnotations;
    // TODO: instead of data, should be specific BoxPrimitive | Sphere | Cylinder | Ellipsoid | PyramidPrimitive
    // selected based on params.segmentId
    const p = data.shape_primitive_list.find(s => s.id === params.segmentId);
    builder.currentGroup = 0;
    switch (p!.kind) {
        case 'sphere':
            addSphere(builder, (p as Sphere).center, (p as Sphere).radius, 2);
            break;
        case 'cylinder':
            addCylinder(builder, (p as Cylinder).start as Vec3, (p as Cylinder).end as Vec3, 1, {
                radiusTop: (p as Cylinder).radius_top,
                radiusBottom: (p as Cylinder).radius_bottom,
                bottomCap: true,
                topCap: true,
            });

            break;
        case 'box':
            addBox(
                builder,
                (p as BoxPrimitive).translation as Vec3,
                (p as BoxPrimitive).scaling as Vec3
            );
            break;
        case 'pyramid':
            addTriangularPyramid(
                builder,
                (p as PyramidPrimitive).translation as Vec3,
                (p as PyramidPrimitive).scaling as Vec3
            );
            break;
        case 'ellipsoid':
            addEllipsoid(
                builder,
                (p as Ellipsoid).center as Vec3,
                (p as Ellipsoid).dir_major as Vec3,
                (p as Ellipsoid).dir_minor as Vec3,
                (p as Ellipsoid).radius_scale as Vec3,
                2
            );
            break;
    }
    // }


    return Shape.create(
        'Shape Primitives',
        params,
        MeshBuilder.getMesh(builder),
        // g => Color(data[g].color),
        g => Color.fromHexStyle(_get_target_segment_color_as_hex(segmentAnnotations, p!.id)),
        () => 1,
        // g => data[g].label,
        g => _get_target_segment_name(descriptions, p!.id)
    );
}