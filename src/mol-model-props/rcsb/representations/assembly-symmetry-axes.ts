/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure } from '../../../mol-model/structure';
import { AssemblySymmetryProvider, AssemblySymmetryValue } from '../assembly-symmetry';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { getSymmetrySelectParam } from '../util';
import { RuntimeContext } from '../../../mol-task';
import { Shape } from '../../../mol-model/shape';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../../../mol-repr/shape/representation';

export const AssemblySymmetryAxesParams = {
    ...Mesh.Params,
    axesColor: PD.Color(ColorNames.orange),
    axesScale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 }),
    symmetryIndex: getSymmetrySelectParam(),
}
export type AssemblySymmetryAxesParams = typeof AssemblySymmetryAxesParams
export type AssemblySymmetryAxesProps = PD.Values<AssemblySymmetryAxesParams>

function getAxesMesh(data: AssemblySymmetryValue, props: PD.Values<AssemblySymmetryAxesParams>, mesh?: Mesh) {

    const { symmetryIndex, axesScale } = props

    const rotation_axes = data?.[symmetryIndex]?.rotation_axes
    if (!rotation_axes) return Mesh.createEmpty(mesh)

    const axis = rotation_axes[0]!
    const start = axis.start as Vec3
    const end = axis.end as Vec3
    const radius = (Vec3.distance(start, end) / 300) * axesScale

    const cylinderProps = { radiusTop: radius, radiusBottom: radius }
    const builderState = MeshBuilder.createState(256, 128, mesh)

    for (let i = 0, il = rotation_axes.length; i < il; ++i) {
        const axis = rotation_axes[i]!
        const start = axis.start as Vec3
        const end = axis.end as Vec3
        builderState.currentGroup = i
        addSphere(builderState, start, radius, 2)
        addSphere(builderState, end, radius, 2)
        addCylinder(builderState, start, end, 1, cylinderProps)
    }
    return MeshBuilder.getMesh(builderState)
}

export async function getAssemblySymmetryAxesRepresentation(ctx: RuntimeContext, structure: Structure, params: AssemblySymmetryAxesProps, prev?: ShapeRepresentation<AssemblySymmetryValue, Mesh, Mesh.Params>) {
    const repr = prev || ShapeRepresentation(getAxesShape, Mesh.Utils);
    const data = AssemblySymmetryProvider.getValue(structure).value
    await repr.createOrUpdate(params, data).runInContext(ctx);
    return repr;
}

function getAxesShape(ctx: RuntimeContext, data: AssemblySymmetryValue, props: AssemblySymmetryAxesProps, shape?: Shape<Mesh>) {
    const geo = getAxesMesh(data, props, shape && shape.geometry);
    const getLabel = (groupId: number) => {
        const { symbol, kind, rotation_axes } = data[props.symmetryIndex]
        return `Axes ${groupId + 1} of ${symbol} ${kind} with Order ${rotation_axes![groupId]?.order}`
    }
    return Shape.create('Unitcell', data, geo, () => props.axesColor, () => 1, getLabel)
}