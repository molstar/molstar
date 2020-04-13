/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { RuntimeContext } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { Shape } from '../../../mol-model/shape';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { lociLabel } from '../../../mol-theme/label';
import { addAxes } from '../../../mol-geo/geometry/mesh/builder/axes';
import { addOrientedBox } from '../../../mol-geo/geometry/mesh/builder/box';
import { addEllipsoid } from '../../../mol-geo/geometry/mesh/builder/ellipsoid';
import { Axes3D } from '../../../mol-math/geometry';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { MarkerActions } from '../../../mol-util/marker-action';

export interface OrientationData {
    locis: Loci[]
}

const SharedParams = {
    color: PD.Color(ColorNames.orange),
    scale: PD.Numeric(2, { min: 0.1, max: 10, step: 0.1 })
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
    visuals: PD.MultiSelect(['box'], PD.objectToOptions(OrientationVisuals)),
    color: PD.Color(ColorNames.orange),
    scale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 })
};
export type OrientationParams = typeof OrientationParams
export type OrientationProps = PD.Values<OrientationParams>

//

function orientationLabel(loci: Loci) {
    const label = lociLabel(loci, { countsOnly: true });
    return `Principal Axes of ${label}`;
}

function getOrientationName(data: OrientationData) {
    return data.locis.length === 1 ? orientationLabel(data.locis[0]) : `${data.locis.length} Orientations`;
}

//

function buildAxesMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh);
    for (let i = 0, il = data.locis.length; i < il; ++i) {
        const principalAxes = Loci.getPrincipalAxes(data.locis[i]);
        if (principalAxes) {
            state.currentGroup = i;
            addAxes(state, principalAxes.momentsAxes, props.scale, 2, 20);
        }
    }
    return MeshBuilder.getMesh(state);
}

function getAxesShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const mesh = buildAxesMesh(data, props, shape && shape.geometry);
    const name = getOrientationName(data);
    const getLabel = function (groupId: number ) {
        return orientationLabel(data.locis[groupId]);
    };
    return Shape.create(name, data, mesh, () => props.color, () => 1, getLabel);
}

//

function buildBoxMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh);
    for (let i = 0, il = data.locis.length; i < il; ++i) {
        const principalAxes = Loci.getPrincipalAxes(data.locis[i]);
        if (principalAxes) {
            state.currentGroup = i;
            addOrientedBox(state, principalAxes.boxAxes, props.scale, 2, 20);
        }
    }
    return MeshBuilder.getMesh(state);
}

function getBoxShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const mesh = buildBoxMesh(data, props, shape && shape.geometry);
    const name = getOrientationName(data);
    const getLabel = function (groupId: number ) {
        return orientationLabel(data.locis[groupId]);
    };
    return Shape.create(name, data, mesh, () => props.color, () => 1, getLabel);
}

//

function buildEllipsoidMesh(data: OrientationData, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh);
    for (let i = 0, il = data.locis.length; i < il; ++i) {
        const principalAxes = Loci.getPrincipalAxes(data.locis[i]);
        if (principalAxes) {
            const axes = principalAxes.boxAxes;
            const { origin, dirA, dirB } = axes;
            const size = Axes3D.size(Vec3(), axes);
            Vec3.scale(size, size, 0.5);
            const radiusScale = Vec3.create(size[2], size[1], size[0]);

            state.currentGroup = i;
            addEllipsoid(state, origin, dirA, dirB, radiusScale, 2);
        }
    }
    return MeshBuilder.getMesh(state);
}

function getEllipsoidShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const mesh = buildEllipsoidMesh(data, props, shape && shape.geometry);
    const name = getOrientationName(data);
    const getLabel = function (groupId: number ) {
        return orientationLabel(data.locis[groupId]);
    };
    return Shape.create(name, data, mesh, () => props.color, () => 1, getLabel);
}

//

export type OrientationRepresentation = Representation<OrientationData, OrientationParams>
export function OrientationRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, OrientationParams>): OrientationRepresentation {
    const repr = Representation.createMulti('Orientation', ctx, getParams, Representation.StateBuilder, OrientationVisuals as unknown as Representation.Def<OrientationData, OrientationParams>);
    repr.setState({ markerActions: MarkerActions.Highlighting });
    return repr;
}