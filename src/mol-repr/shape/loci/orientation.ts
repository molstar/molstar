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
import { Vec3, Mat4 } from '../../../mol-math/linear-algebra';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { createCage } from '../../../mol-geo/primitive/cage';
import { Axes3D } from '../../../mol-math/geometry';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { lociLabel } from '../../../mol-theme/label';

export interface OrientationData {
    loci: Loci
}

const SharedParams = {
    color: PD.Color(ColorNames.orange),
    scale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 })
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

const OrientationVisuals = {
    'axes': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, AxesParams>) => ShapeRepresentation(getAxesShape, Mesh.Utils),
    'box': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, BoxParams>) => ShapeRepresentation(getBoxShape, Mesh.Utils),
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

const tmpAxesVec = Vec3()

function buildAxesMesh(principalAxes: PrincipalAxes, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh)

    const { origin, dirA, dirB, dirC } = principalAxes.momentsAxes

    const vertices = new Float32Array(6 * 3)
    const edges = new Uint8Array([0, 1, 2, 3, 4, 5])
    Vec3.add(tmpAxesVec, origin, dirA)
    Vec3.toArray(Vec3.add(tmpAxesVec, origin, dirA), vertices, 0)
    Vec3.toArray(Vec3.sub(tmpAxesVec, origin, dirA), vertices, 3)
    Vec3.toArray(Vec3.add(tmpAxesVec, origin, dirB), vertices, 6)
    Vec3.toArray(Vec3.sub(tmpAxesVec, origin, dirB), vertices, 9)
    Vec3.toArray(Vec3.add(tmpAxesVec, origin, dirC), vertices, 12)
    Vec3.toArray(Vec3.sub(tmpAxesVec, origin, dirC), vertices, 15)

    const matrix = Mat4.identity()

    const cage = createCage(vertices, edges)
    const volume = Axes3D.volume(principalAxes.boxAxes)
    const radius = (Math.cbrt(volume) / 300) * props.scale
    state.currentGroup = 1
    MeshBuilder.addCage(state, matrix, cage, radius, 2, 20)

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

const tmpBoxVecCorner = Vec3()
const tmpBoxVecA = Vec3()
const tmpBoxVecB = Vec3()
const tmpBoxVecC = Vec3()

function buildBoxMesh(principalAxes: PrincipalAxes, props: OrientationProps, mesh?: Mesh): Mesh {
    const state = MeshBuilder.createState(256, 128, mesh)

    const { origin, dirA, dirB, dirC } = principalAxes.boxAxes
    const negDirA = Vec3.negate(tmpBoxVecA, dirA)
    const negDirB = Vec3.negate(tmpBoxVecB, dirB)
    const negDirC = Vec3.negate(tmpBoxVecC, dirC)

    const vertices = new Float32Array(8 * 3)
    const edges = new Uint8Array([
        0, 1, 0, 3, 0, 6, 1, 2, 1, 7, 2, 3,
        2, 4, 3, 5, 4, 5, 4, 7, 5, 6, 6, 7
    ])

    let offset = 0
    const addCornerHelper = function (v1: Vec3, v2: Vec3, v3: Vec3) {
        Vec3.copy(tmpBoxVecCorner, origin)
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v1)
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v2)
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v3)
        Vec3.toArray(tmpBoxVecCorner, vertices, offset)
        offset += 3
    }
    addCornerHelper(dirA, dirB, dirC)
    addCornerHelper(dirA, dirB, negDirC)
    addCornerHelper(dirA, negDirB, negDirC)
    addCornerHelper(dirA, negDirB, dirC)
    addCornerHelper(negDirA, negDirB, negDirC)
    addCornerHelper(negDirA, negDirB, dirC)
    addCornerHelper(negDirA, dirB, dirC)
    addCornerHelper(negDirA, dirB, negDirC)

    const matrix = Mat4.identity()

    const cage = createCage(vertices, edges)
    const volume = Axes3D.volume(principalAxes.boxAxes)
    const radius = (Math.cbrt(volume) / 300) * props.scale
    state.currentGroup = 1
    MeshBuilder.addCage(state, matrix, cage, radius, 2, 20)

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

export type OrientationRepresentation = Representation<OrientationData, OrientationParams>
export function OrientationRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<OrientationData, OrientationParams>): OrientationRepresentation {
    return Representation.createMulti('Orientation', ctx, getParams, Representation.StateBuilder, OrientationVisuals as unknown as Representation.Def<OrientationData, OrientationParams>)
}