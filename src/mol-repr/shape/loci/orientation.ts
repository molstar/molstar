/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { RuntimeContext } from '../../../mol-task';
import { stringToWords } from '../../../mol-util/string';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { Shape } from '../../../mol-model/shape';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { lociLabel } from '../../../mol-theme/label';
import { addAxes } from '../../../mol-geo/geometry/mesh/builder/axes';
import { addOrientedBox } from '../../../mol-geo/geometry/mesh/builder/box';
import { addEllipsoid } from '../../../mol-geo/geometry/mesh/builder/ellipsoid';
import { Axes3D } from '../../../mol-math/geometry';
import { Vec3 } from '../../../mol-math/linear-algebra';

export interface OrientationData {
    loci: Loci
}

const SharedParams = {
    color: PD.Color(ColorNames.orange),
    scale: PD.Numeric(2, { min: 0.1, max: 10, step: 0.1 })
}

const AxesParams = {
    ...Mesh.Params,
    ...SharedParams
}
type AxesParams = typeof AxesParams

const BoxParams = {
    ...Mesh.Params,
    ...SharedParams
}
type BoxParams = typeof BoxParams

const EllipsoidParams = {
    ...Mesh.Params,
    ...SharedParams
}
type EllipsoidParams = typeof EllipsoidParams

const OrientationVisuals = {
    'axes': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, AxesParams>) => ShapeRepresentation(getAxesShape, Mesh.Utils),
    'box': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, BoxParams>) => ShapeRepresentation(getBoxShape, Mesh.Utils),
    'ellipsoid': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, EllipsoidParams>) => ShapeRepresentation(getEllipsoidShape, Mesh.Utils),
}
type OrientationVisualName = keyof typeof OrientationVisuals
const OrientationVisualOptions = Object.keys(OrientationVisuals).map(name => [name, stringToWords(name)] as [OrientationVisualName, string])

export const OrientationParams = {
    ...AxesParams,
    ...BoxParams,
    visuals: PD.MultiSelect<OrientationVisualName>(['box'], OrientationVisualOptions),
    color: PD.Color(ColorNames.orange),
    scale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 })
}
export type OrientationParams = typeof OrientationParams
export type OrientationProps = PD.Values<OrientationParams>

//

function buildAxesMesh(principalAxes: PrincipalAxes, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh)
    state.currentGroup = 1
    addAxes(state, principalAxes.momentsAxes, props.scale, 2, 20)
    return MeshBuilder.getMesh(state)
}

function getAxesShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const label = lociLabel(data.loci, { countsOnly: true })
    const principalAxes = Loci.getPrincipalAxes(data.loci)
    const mesh = principalAxes ? buildAxesMesh(principalAxes, props, shape && shape.geometry) : Mesh.createEmpty(shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return `Principal Axes of ${label}`
    }
    return Shape.create('Principal Axes', data, mesh, () => props.color, () => 1, getLabel)
}

//

function buildBoxMesh(principalAxes: PrincipalAxes, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh)
    state.currentGroup = 1
    addOrientedBox(state, principalAxes.boxAxes, props.scale, 2, 20)
    return MeshBuilder.getMesh(state)
}

function getBoxShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const label = lociLabel(data.loci, { countsOnly: true })
    const principalAxes = Loci.getPrincipalAxes(data.loci)
    const mesh = principalAxes ? buildBoxMesh(principalAxes, props, shape && shape.geometry) : Mesh.createEmpty(shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return `Oriented Box of ${label}`
    }
    return Shape.create('Oriented Box', data, mesh, () => props.color, () => 1, getLabel)
}

//

function buildEllipsoidMesh(principalAxes: PrincipalAxes, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh)

    const axes = principalAxes.boxAxes
    const { origin, dirA, dirB } = axes
    const size = Axes3D.size(Vec3(), axes)
    Vec3.scale(size, size, 0.5)
    const radiusScale = Vec3.create(size[2], size[1], size[0])

    state.currentGroup = 1
    addEllipsoid(state, origin, dirA, dirB, radiusScale, 2)
    return MeshBuilder.getMesh(state)
}

function getEllipsoidShape(ctx: RuntimeContext, data: OrientationData, props: OrientationProps, shape?: Shape<Mesh>) {
    const label = lociLabel(data.loci, { countsOnly: true })
    const principalAxes = Loci.getPrincipalAxes(data.loci)
    const mesh = principalAxes ? buildEllipsoidMesh(principalAxes, props, shape && shape.geometry) : Mesh.createEmpty(shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return `Ellipsoid of ${label}`
    }
    return Shape.create('Ellipsoid', data, mesh, () => props.color, () => 1, getLabel)
}

//

export type OrientationRepresentation = Representation<OrientationData, OrientationParams>
export function OrientationRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, OrientationParams>): OrientationRepresentation {
    return Representation.createMulti('Orientation', ctx, getParams, Representation.StateBuilder, OrientationVisuals as unknown as Representation.Def<OrientationData, OrientationParams>)
}