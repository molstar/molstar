/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure } from '../../../mol-model/structure';
import { AssemblySymmetryProvider, AssemblySymmetryValue, getSymmetrySelectParam } from '../assembly-symmetry';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3, Mat4 } from '../../../mol-math/linear-algebra';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { RuntimeContext } from '../../../mol-task';
import { Shape } from '../../../mol-model/shape';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../../../mol-repr/shape/representation';
import { MarkerActions } from '../../../mol-util/marker-action';
import { Prism } from '../../../mol-geo/primitive/prism';
import { Wedge } from '../../../mol-geo/primitive/wedge';
import { Primitive, transformPrimitive } from '../../../mol-geo/primitive/primitive';
import { memoize1 } from '../../../mol-util/memoize';
import { polygon } from '../../../mol-geo/primitive/polygon';
import { ColorMap, Color } from '../../../mol-util/color';
import { TableLegend } from '../../../mol-util/legend';

const OrderColors = ColorMap({
    '2': ColorNames.deepskyblue,
    '3': ColorNames.lime,
    'N': ColorNames.red,
})
const OrderColorsLegend = TableLegend(Object.keys(OrderColors).map(name => {
    return [name, (OrderColors as any)[name] as Color] as [string, Color]
}))

function axesColorHelp(value: { name: string, params: {} }) {
    return value.name === 'byOrder'
        ? { description: 'Color axes by their order', legend: OrderColorsLegend }
        : {}
}

export const AssemblySymmetryAxesParams = {
    ...Mesh.Params,
    axesColor: PD.MappedStatic('byOrder', {
        byOrder: PD.EmptyGroup(),
        uniform: PD.Group({
            colorValue: PD.Color(ColorNames.orange),
        }, { isFlat: true })
    }, { help: axesColorHelp }),
    axesScale: PD.Numeric(2, { min: 0.1, max: 5, step: 0.1 }),
    symmetryIndex: getSymmetrySelectParam(),
}
export type AssemblySymmetryAxesParams = typeof AssemblySymmetryAxesParams
export type AssemblySymmetryAxesProps = PD.Values<AssemblySymmetryAxesParams>

const t = Mat4.identity()
const upY = Vec3.create(0, 1, 0)
const upX = Vec3.create(1, 0, 0)
const tmpV = Vec3()
const center = Vec3()
const scale = Vec3()

const getPrimitive = memoize1((order: number): Primitive | undefined => {
    if (order < 2) {
        return Prism(polygon(48, false))
    } else if (order === 2) {
        const lens = Prism(polygon(48, false))
        Mat4.setIdentity(t)
        Mat4.scale(t, t, Vec3.set(scale, 1, 0.35, 1))
        transformPrimitive(lens, t)
        return lens
    } else if (order === 3) {
        return Wedge()
    } else {
        return Prism(polygon(order, false))
    }
})

function getAxesMesh(data: AssemblySymmetryValue, props: PD.Values<AssemblySymmetryAxesParams>, mesh?: Mesh) {

    const { symmetryIndex, axesScale } = props

    const rotation_axes = data?.[symmetryIndex]?.rotation_axes
    if (!rotation_axes) return Mesh.createEmpty(mesh)

    const axis = rotation_axes[0]!
    const start = axis.start as Vec3
    const end = axis.end as Vec3
    const radius = (Vec3.distance(start, end) / 500) * axesScale

    const cylinderProps = { radiusTop: radius, radiusBottom: radius }
    const builderState = MeshBuilder.createState(256, 128, mesh)

    for (let i = 0, il = rotation_axes.length; i < il; ++i) {
        const axis = rotation_axes[i]!
        const start = axis.start as Vec3
        const end = axis.end as Vec3
        builderState.currentGroup = i
        addCylinder(builderState, start, end, 1, cylinderProps)

        const primitive = getPrimitive(axis.order!)
        if (primitive) {
            Vec3.scale(center, Vec3.add(center, start, end), 0.5)
            if (Vec3.dot(upY, Vec3.sub(tmpV, start, center)) === 0) {
                Mat4.targetTo(t, start, center, upY)
            } else {
                Mat4.targetTo(t, start, center, upX)
            }
            Mat4.scale(t, t, Vec3.set(scale, radius * 7, radius * 7, radius * 0.4))

            Mat4.setTranslation(t, start)
            MeshBuilder.addPrimitive(builderState, t, primitive)
            Mat4.setTranslation(t, end)
            MeshBuilder.addPrimitive(builderState, t, primitive)
        }
    }
    return MeshBuilder.getMesh(builderState)
}

export async function getAssemblySymmetryAxesRepresentation(ctx: RuntimeContext, structure: Structure, params: AssemblySymmetryAxesProps, prev?: ShapeRepresentation<AssemblySymmetryValue, Mesh, Mesh.Params>) {
    const repr = prev || ShapeRepresentation(getAxesShape, Mesh.Utils);
    const data = AssemblySymmetryProvider.get(structure).value
    await repr.createOrUpdate(params, data).runInContext(ctx);
    repr.setState({ markerActions: MarkerActions.Highlighting })
    return repr;
}

function getAxesShape(ctx: RuntimeContext, data: AssemblySymmetryValue, props: AssemblySymmetryAxesProps, shape?: Shape<Mesh>) {
    const geo = getAxesMesh(data, props, shape && shape.geometry);
    const getColor = (groupId: number) => {
        if (props.axesColor.name === 'byOrder') {
            const { rotation_axes } = data[props.symmetryIndex]
            const order = rotation_axes![groupId]?.order
            if (order === 2) return OrderColors[2]
            else if (order === 3) return OrderColors[3]
            else return OrderColors.N
        } else {
            return props.axesColor.params.colorValue
        }
    }
    const getLabel = (groupId: number) => {
        const { symbol, kind, rotation_axes } = data[props.symmetryIndex]
        return `Axes ${groupId + 1} of ${symbol} ${kind} with Order ${rotation_axes![groupId]?.order}`
    }
    return Shape.create('Unitcell', data, geo, getColor, () => 1, getLabel)
}