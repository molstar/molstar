/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../mol-repr/representation';
import { Structure } from '../../mol-model/structure';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { StructureRepresentationProvider, StructureRepresentation, StructureRepresentationStateBuilder } from '../../mol-repr/structure/representation';
import { MembraneOrientation } from './membrane-orientation';
import { ThemeRegistryContext } from '../../mol-theme/theme';
import { ShapeRepresentation } from '../../mol-repr/shape/representation';
import { Shape } from '../../mol-model/shape';
import { ColorNames } from '../../mol-util/color/names';
import { RuntimeContext } from '../../mol-task';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
import { Circle } from '../../mol-geo/primitive/circle';
import { transformPrimitive } from '../../mol-geo/primitive/primitive';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { MembraneOrientationProvider } from './membrane-orientation';

const SharedParams = {
    color: PD.Color(ColorNames.lightgrey),
    radiusFactor: PD.Numeric(0.8333, { min: 0.1, max: 3.0, step: 0.01 }, { description: 'Scale the radius of the membrane layer' })
};

const BilayerSpheresParams = {
    ...SharedParams,
    ...Spheres.Params,
    sphereSize: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }, { description: 'Size of spheres that represent membrane planes' }),
    density: PD.Numeric(1, { min: 0.25, max: 10, step: 0.25 }, { description: 'Distance between spheres'})
};
export type BilayerSpheresParams = typeof BilayerSpheresParams
export type BilayerSpheresProps = PD.Values<BilayerSpheresParams>

const BilayerPlanesParams = {
    ...SharedParams,
    ...Mesh.Params,
    sectorOpacity: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 })
};
export type BilayerPlanesParams = typeof BilayerPlanesParams
export type BilayerPlanesProps = PD.Values<BilayerPlanesParams>

const BilayerRimsParams = {
    ...SharedParams,
    ...Lines.Params,
    lineSizeAttenuation: PD.Boolean(true),
    linesSize: PD.Numeric(0.5, { min: 0.01, max: 50, step: 0.01 }),
    dashedLines: PD.Boolean(true)
};
export type BilayerRimsParams = typeof BilayerRimsParams
export type BilayerRimsProps = PD.Values<BilayerRimsParams>

const MembraneOrientationVisuals = {
    'bilayer-spheres': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<MembraneOrientation, BilayerSpheresParams>) => ShapeRepresentation(getBilayerSpheres, Spheres.Utils, { modifyState: s => ({ ...s, pickable: false }) }),
    'bilayer-planes': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<MembraneOrientation, BilayerPlanesParams>) => ShapeRepresentation(getBilayerPlanes, Mesh.Utils, { modifyProps: p => ({ ...p, alpha: p.sectorOpacity }), modifyState: s => ({ ...s, pickable: false }) }),
    'bilayer-rims': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<MembraneOrientation, BilayerRimsParams>) => ShapeRepresentation(getBilayerRims, Lines.Utils, { modifyState: s => ({ ...s, pickable: false }) })
};

export const MembraneOrientationParams = {
    ...BilayerSpheresParams,
    ...BilayerPlanesParams,
    ...BilayerRimsParams,
    visuals: PD.MultiSelect(['bilayer-planes', 'bilayer-rims'], PD.objectToOptions(MembraneOrientationVisuals)),
};
export type MembraneOrientationParams = typeof MembraneOrientationParams
export type MembraneOrientationProps = PD.Values<MembraneOrientationParams>

export function getMembraneOrientationParams(ctx: ThemeRegistryContext, structure: Structure) {
    console.log('rep @ getParams');
    return PD.clone(MembraneOrientationParams);
}

export type MembraneOrientationRepresentation = StructureRepresentation<MembraneOrientationParams>
export function MembraneOrientationRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MembraneOrientationParams>): MembraneOrientationRepresentation {
    console.log('rep @ create');
    return Representation.createMulti('Membrane Orientation', ctx, getParams, StructureRepresentationStateBuilder, MembraneOrientationVisuals as unknown as Representation.Def<Structure, MembraneOrientationParams>);
}

export const MembraneOrientationRepresentationProvider = StructureRepresentationProvider({
    name: 'membrane-orientation',
    label: 'Membrane Orientation',
    description: 'Displays a grid of points representing membrane layers.',
    factory: MembraneOrientationRepresentation,
    getParams: getMembraneOrientationParams,
    defaultValues: PD.getDefaultValues(MembraneOrientationParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});

function getBilayerRims(ctx: RuntimeContext, data: Structure, props: BilayerRimsProps): Shape<Lines> {
    console.log('rims');
    const { p1, p2, centroid, normal, radius } = MembraneOrientationProvider.get(data).value!;
    const scaledRadius = props.radiusFactor * radius;
    const builder = LinesBuilder.create(128, 64);
    getLayerCircle(builder, p1, centroid, normal, scaledRadius, props);
    getLayerCircle(builder, p2, centroid, normal, scaledRadius, props);
    return Shape.create(name, data, builder.getLines(), () => props.color, () => props.linesSize, () => '');
}

function getLayerCircle(builder: LinesBuilder, p: Vec3, centroid: Vec3, normal: Vec3, radius: number, props: BilayerRimsProps) {
    const circle = getCircle(p, centroid, normal, radius);
    const { indices, vertices } = circle;
    for (let j = 0, jl = indices.length; j < jl; j += 3) {
        if (props.dashedLines && j % 2 === 1) continue; // draw every other segment to get dashes
        const start = indices[j] * 3;
        const end = indices[j + 1] * 3;
        const startX = vertices[start];
        const startY = vertices[start + 1];
        const startZ = vertices[start + 2];
        const endX = vertices[end];
        const endY = vertices[end + 1];
        const endZ = vertices[end + 2];
        builder.add(startX, startY, startZ, endX, endY, endZ, 0);
    }
}

const tmpMat = Mat4();
function getCircle(p: Vec3, centroid: Vec3, normal: Vec3, radius: number) {
    Mat4.targetTo(tmpMat, p, centroid, normal);
    Mat4.setTranslation(tmpMat, p);
    Mat4.mul(tmpMat, tmpMat, Mat4.rotX90);

    const circle = Circle({ radius });
    return transformPrimitive(circle, tmpMat);
}

function getBilayerPlanes(ctx: RuntimeContext, data: Structure, props: BilayerPlanesProps, shape?: Shape<Mesh>): Shape<Mesh> {
    console.log('planes');
    const { p1, p2, centroid, normal, radius } = MembraneOrientationProvider.get(data).value!;
    const state = MeshBuilder.createState(128, 64, shape && shape.geometry);
    const scaledRadius = props.radiusFactor * radius;
    getLayerPlane(state, p1, centroid, normal, scaledRadius);
    getLayerPlane(state, p2, centroid, normal, scaledRadius);
    return Shape.create(name, data, MeshBuilder.getMesh(state), () => props.color, () => 1, () => '');
}

function getLayerPlane(state: MeshBuilder.State, p: Vec3, centroid: Vec3, normal: Vec3, radius: number) {
    const circle = getCircle(p, centroid, normal, radius);
    state.currentGroup = 0;
    MeshBuilder.addPrimitive(state, Mat4.id, circle);
    MeshBuilder.addPrimitiveFlipped(state, Mat4.id, circle);
}

function getBilayerSpheres(ctx: RuntimeContext, data: Structure, props: BilayerSpheresProps): Shape<Spheres> {
    const { density } = props;
    const { radius, p1, p2, normal } = MembraneOrientationProvider.get(data).value!;
    const scaledRadius = (props.radiusFactor * radius) * (props.radiusFactor * radius);

    const spheresBuilder = SpheresBuilder.create();
    getLayerSpheres(spheresBuilder, p1, normal, density, scaledRadius);
    getLayerSpheres(spheresBuilder, p2, normal, density, scaledRadius);
    return Shape.create(name, data, spheresBuilder.getSpheres(), () => props.color, () => props.sphereSize, () => '');
}

function getLayerSpheres(spheresBuilder: SpheresBuilder, point: Vec3, normalVector: Vec3, density: number, sqRadius: number) {
    Vec3.normalize(normalVector, normalVector);
    const d = -Vec3.dot(normalVector, point);
    const rep = Vec3();
    for (let i = -1000, il = 1000; i < il; i += density) {
        for (let j = -1000, jl = 1000; j < jl; j += density) {
            Vec3.set(rep, i, j, normalVector[2] === 0 ? 0 : -(d + i * normalVector[0] + j * normalVector[1]) / normalVector[2]);
            if (Vec3.squaredDistance(rep, point) < sqRadius) {
                spheresBuilder.add(rep[0], rep[1], rep[2], 0);
            }
        }
    }
}