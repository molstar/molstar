/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { ShapeRepresentation } from '../../mol-repr/shape/representation';
import { Shape } from '../../mol-model/shape';
import { ColorNames } from '../../mol-util/color/names';
import { RuntimeContext } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
// import { Polyhedron, DefaultPolyhedronProps } from '../../mol-geo/primitive/polyhedron';
// import { Icosahedron } from '../../mol-geo/primitive/icosahedron';
import { Sphere } from '../../mol-geo/primitive/sphere';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { RepresentationParamsGetter, Representation, RepresentationContext } from '../../mol-repr/representation';


interface MembraneSphereData {
    radius: number
    center: Vec3
}


const MembraneSphereParams = {
    ...Mesh.Params,
    cellColor: PD.Color(ColorNames.orange),
    cellScale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 }),
    radius: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 }),
    center: PD.Vec3(Vec3.create(0, 0, 0)),
    quality: { ...Mesh.Params.quality, isEssential: false },
};

type MeshParams = typeof MembraneSphereParams

const MembraneSphereVisuals = {
    'mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<MembraneSphereData, MeshParams>) => ShapeRepresentation(getMBShape, Mesh.Utils),
};

export const MBParams = {
    ...MembraneSphereParams
};
export type MBParams = typeof MBParams
export type UnitcellProps = PD.Values<MBParams>

function getMBMesh(data: MembraneSphereData, props: UnitcellProps, mesh?: Mesh) {
    const state = MeshBuilder.createState(256, 128, mesh);
    const radius = props.radius;
    const asphere = Sphere(3);
    const trans: Mat4 = Mat4.identity();
    Mat4.fromScaling(trans, Vec3.create(radius, radius, radius));
    state.currentGroup = 1;
    MeshBuilder.addPrimitive(state, trans, asphere);
    const m = MeshBuilder.getMesh(state);
    return m;
}

function getMBShape(ctx: RuntimeContext, data: MembraneSphereData, props: UnitcellProps, shape?: Shape<Mesh>) {
    const geo = getMBMesh(data, props, shape && shape.geometry);
    const label = 'mb';
    return Shape.create(label, data, geo, () => props.cellColor, () => 1, () => label);
}

export type MBRepresentation = Representation<MembraneSphereData, MBParams>
export function MBRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<MembraneSphereData, MBParams>): MBRepresentation {
    return Representation.createMulti('MB', ctx, getParams, Representation.StateBuilder, MembraneSphereVisuals as unknown as Representation.Def<MembraneSphereData, MBParams>);
}