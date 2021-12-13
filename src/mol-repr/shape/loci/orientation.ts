/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { Shape } from '../../../mol-model/shape';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { structureElementLociLabelMany } from '../../../mol-theme/label';
import { addAxes } from '../../../mol-geo/geometry/mesh/builder/axes';
import { addOrientedBox } from '../../../mol-geo/geometry/mesh/builder/box';
import { addEllipsoid } from '../../../mol-geo/geometry/mesh/builder/ellipsoid';
import { Axes3D } from '../../../mol-math/geometry';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { MarkerActions } from '../../../mol-util/marker-action';
import { StructureElement } from '../../../mol-model/structure';

export interface OrientationData {
    locis: StructureElement.Loci[]
}

const SharedParams = {
    color: PD.Color(ColorNames.orange),
    scaleFactor: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }),
    radiusScale: PD.Numeric(2, { min: 0.1, max: 10, step: 0.1 })
};

const AxesParams = {
    ...Mesh.Params,
    ...SharedParams
};
type AxesParams = typeof AxesParams

const BoxParams = {
    ...Mesh.Params,
    ...SharedParams
};
type BoxParams = typeof BoxParams

const EllipsoidParams = {
    ...Mesh.Params,
    ...SharedParams
};
type EllipsoidParams = typeof EllipsoidParams

const OrientationVisuals = {
    'axes': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, AxesParams>) => ShapeRepresentation(getAxesShape, Mesh.Utils),
    'box': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, BoxParams>) => ShapeRepresentation(getBoxShape, Mesh.Utils),
    'ellipsoid': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, EllipsoidParams>) => ShapeRepresentation(getEllipsoidShape, Mesh.Utils),
};

export const OrientationParams = {
    ...AxesParams,
    ...BoxParams,
    ...EllipsoidParams,
    visuals: PD.MultiSelect(['box'], PD.objectToOptions(OrientationVisuals)),
};
export type OrientationParams = typeof OrientationParams
export type OrientationProps = PD.Values<OrientationParams>

//

function getAxesName(locis: StructureElement.Loci[]) {
    const label = structureElementLociLabelMany(locis, { countsOnly: true });
    return `Principal Axes of ${label}`;
}

function buildAxesMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh);
    const principalAxes = StructureElement.Loci.getPrincipalAxesMany(data.locis);
    Axes3D.scale(principalAxes.momentsAxes, principalAxes.momentsAxes, props.scaleFactor);

    state.currentGroup = 0;
    addAxes(state, principalAxes.momentsAxes, props.radiusScale, 2, 20);
    return MeshBuilder.getMesh(state);
}

function getAxesShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const mesh = buildAxesMesh(data, props, shape && shape.geometry);
    const name = getAxesName(data.locis);
    return Shape.create(name, data, mesh, () => props.color, () => 1, () => name);
}

//

function getBoxName(locis: StructureElement.Loci[]) {
    const label = structureElementLociLabelMany(locis, { countsOnly: true });
    return `Oriented Box of ${label}`;
}

function buildBoxMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh);
    const principalAxes = StructureElement.Loci.getPrincipalAxesMany(data.locis);
    Axes3D.scale(principalAxes.boxAxes, principalAxes.boxAxes, props.scaleFactor);

    state.currentGroup = 0;
    addOrientedBox(state, principalAxes.boxAxes, props.radiusScale, 2, 20);
    return MeshBuilder.getMesh(state);
}

function getBoxShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const mesh = buildBoxMesh(data, props, shape && shape.geometry);
    const name = getBoxName(data.locis);
    return Shape.create(name, data, mesh, () => props.color, () => 1, () => name);
}

//

function getEllipsoidName(locis: StructureElement.Loci[]) {
    const label = structureElementLociLabelMany(locis, { countsOnly: true });
    return `Oriented Ellipsoid of ${label}`;
}

function buildEllipsoidMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh);
    const principalAxes = StructureElement.Loci.getPrincipalAxesMany(data.locis);

    const axes = principalAxes.boxAxes;
    const { origin, dirA, dirB } = axes;
    const size = Axes3D.size(Vec3(), axes);
    Vec3.scale(size, size, 0.5 * props.scaleFactor);
    const radiusScale = Vec3.create(size[2], size[1], size[0]);

    state.currentGroup = 0;
    addEllipsoid(state, origin, dirA, dirB, radiusScale, 2);
    return MeshBuilder.getMesh(state);
}

function getEllipsoidShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const mesh = buildEllipsoidMesh(data, props, shape && shape.geometry);
    const name = getEllipsoidName(data.locis);
    return Shape.create(name, data, mesh, () => props.color, () => 1, () => name);
}

//

export type OrientationRepresentation = Representation<OrientationData, OrientationParams>
export function OrientationRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, OrientationParams>): OrientationRepresentation {
    const repr = Representation.createMulti('Orientation', ctx, getParams, Representation.StateBuilder, OrientationVisuals as unknown as Representation.Def<OrientationData, OrientationParams>);
    repr.setState({ markerActions: MarkerActions.Highlighting });
    return repr;
}