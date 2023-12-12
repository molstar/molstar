// import { StateTransformer } from '../../mol-state';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Shape } from '../../mol-model/shape';
import { Color } from '../../mol-util/color';
import { StateTransformer } from '../../mol-state';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { Box } from '../../mol-geo/primitive/box';
import { Pyramid, TriangularPyramid } from '../../mol-geo/primitive/pyramid';
import { Primitive } from '../../mol-geo/primitive/primitive';
import { polygon } from '../../mol-geo/primitive/polygon';
import { addEllipsoid } from '../../mol-geo/geometry/mesh/builder/ellipsoid';


let cone: Primitive;
export function Cone() {
    if (!cone) cone = Pyramid(polygon(48, true));
    return cone;
}


export type ShapePrimitive =
    | { kind: 'sphere', center: number[], radius: number, label: string, color: number }
    | { kind: 'cylinder', start: number[], end: number[], radius: number, label: string, color: number }
    | { kind: 'box', translation: number[], scaling: number[], label: string, color: number }
    | { kind: 'pyramid', translation: number[], scaling: number[], label: string, color: number }
    | { kind: 'cone', translation: number[], scaling: number[], label: string, color: number }
    | { kind: 'ellipsoid', dir_major: number[], dir_minor: number[], center: number[], radius_scale: number[], label: string, color: number }

export type ShapePrimitivesData = ShapePrimitive[]

function addBox(state: MeshBuilder.State,
    translation: Vec3 = [0.5, 0.5, 0.5] as Vec3,
    scaling: Vec3 = [1, 1, 1] as Vec3) {
    const mat4 = Mat4.identity();
    Mat4.scale(mat4, mat4, scaling);
    Mat4.translate(mat4, mat4, translation);

    MeshBuilder.addPrimitive(state, mat4, Box());
}

function addTriangularPyramid(state: MeshBuilder.State,
    translation: Vec3 = [0.5, 0.5, 0.5] as Vec3,
    scaling: Vec3 = [1, 1, 1] as Vec3) {
    const mat4 = Mat4.identity();
    Mat4.scale(mat4, mat4, scaling);
    Mat4.translate(mat4, mat4, translation);
    MeshBuilder.addPrimitive(state, mat4, TriangularPyramid());
}

function addCone(state: MeshBuilder.State,
    translation: Vec3 = [0.5, 0.5, 0.5] as Vec3,
    scaling: Vec3 = [1, 1, 1] as Vec3) {
    const mat4 = Mat4.identity();
    Mat4.scale(mat4, mat4, scaling);
    Mat4.translate(mat4, mat4, translation);
    MeshBuilder.addPrimitive(state, mat4, Cone());
}

const Transform = StateTransformer.builderFactory('msvolseg');
export const CreateShapePrimitivesProvider = Transform({
    name: 'create-sphere-provider',
    display: { name: 'Spheres' },
    from: PluginStateObject.Root, // EntryData
    to: PluginStateObject.Shape.Provider,
    params: {
        data: PD.Value<ShapePrimitivesData>([] as any, { isHidden: true })
    }
})({
    apply({ params }) {
        return new PluginStateObject.Shape.Provider({
            label: 'Shape Primitives',
            data: params.data,
            params: Mesh.Params,
            geometryUtils: Mesh.Utils,
            getShape: (_, data) => createShapePrimitives(data)
        }, { label: 'Shape Primitives' });
    }
});

function createShapePrimitives(data: ShapePrimitivesData) {

    const builder = MeshBuilder.createState(512, 512);

    for (let i = 0; i < data.length; i++) {
        const p = data[i];
        builder.currentGroup = i;

        switch (p.kind) {
            case 'sphere':
                addSphere(builder, p.center as Vec3, p.radius, 2);
                break;
            case 'cylinder':
                addCylinder(builder, p.start as Vec3, p.end as Vec3, 1, {
                    radiusTop: p.radius,
                    radiusBottom: p.radius,
                    bottomCap: true,
                    topCap: true,
                });

                break;
            case 'box':
                addBox(
                    builder,
                    p.translation as Vec3,
                    p.scaling as Vec3
                );
                break;
            case 'pyramid':
                addTriangularPyramid(
                    builder,
                    p.translation as Vec3,
                    p.scaling as Vec3
                );
                break;
            case 'cone':
                addCone(
                    builder,
                    p.translation as Vec3,
                    p.scaling as Vec3
                );
                break;
            case 'ellipsoid':
                addEllipsoid(
                    builder,
                    p.center as Vec3,
                    p.dir_major as Vec3,
                    p.dir_minor as Vec3,
                    p.radius_scale as Vec3,
                    2
                );
                break;
        }
    }

    return Shape.create(
        'Shape Primitives',
        {},
        MeshBuilder.getMesh(builder),
        g => Color(data[g].color),
        () => 1,
        g => data[g].label,
    );
}