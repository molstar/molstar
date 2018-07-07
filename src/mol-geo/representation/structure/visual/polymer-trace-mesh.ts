/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Element } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '../index';
import { RuntimeContext } from 'mol-task'
import { createTransforms, createColors } from './util/common';
import { markElement } from './util/element';
import { deepEqual } from 'mol-util';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { OrderedSet } from 'mol-data/int';
import { createMarkers, MarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerElementCount, PolymerTraceIterator } from './util/polymer';
import { Vec3 } from 'mol-math/linear-algebra';

export function reflect(target: Vec3, p1: Vec3, p2: Vec3, amount: number) {
    target[0] = p1[0] - amount * (p2[0] - p1[0])
    target[1] = p1[1] - amount * (p2[1] - p1[1])
    target[2] = p1[2] - amount * (p2[2] - p1[2])
}

async function createPolymerTraceMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    console.log('polymerElementCount', polymerElementCount)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)
    const linearSegments = 5
    const radialSegments = 8
    const tension = 0.9

    const tA = Vec3.zero()
    const tB = Vec3.zero()
    const dA = Vec3.zero()
    const dB = Vec3.zero()
    const torsionVec = Vec3.zero()
    const initialTorsionVec = Vec3.zero()
    const tangentVec = Vec3.zero()
    const normalVec = Vec3.zero()

    const tmp = Vec3.zero()
    const reflectedControlPoint = Vec3.zero()

    const pn = (linearSegments + 1) * 3
    const controlPoints = new Float32Array(pn)
    const torsionVectors = new Float32Array(pn)
    const normalVectors = new Float32Array(pn)

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()
        builder.setId(v.index)

        Vec3.spline(tB, v.t1, v.t2, v.t3, v.t4, 0.5, tension)
        Vec3.spline(dA, v.d12, v.d23, v.d34, v.d45, 0.5, tension)

        Vec3.normalize(initialTorsionVec, Vec3.sub(initialTorsionVec, tB, dB))

        Vec3.toArray(tB, controlPoints, 0)
        Vec3.normalize(torsionVec, Vec3.sub(torsionVec, tB, dB))
        Vec3.toArray(torsionVec, torsionVectors, 0)
        // approximate tangent as direction to previous control point
        Vec3.normalize(tangentVec, Vec3.sub(tangentVec, tB, tA))
        Vec3.normalize(normalVec, Vec3.cross(normalVec, tangentVec, torsionVec))
        Vec3.toArray(normalVec, normalVectors, 0)

        //

        const t12 = Vec3.zero()
        const t23 = Vec3.zero()
        const t34 = Vec3.zero()
        const t45 = Vec3.zero()
        Vec3.spline(t12, v.t0, v.t1, v.t2, v.t3, 0.5, tension)
        Vec3.spline(t23, v.t1, v.t2, v.t3, v.t4, 0.5, tension)
        Vec3.spline(t34, v.t2, v.t3, v.t4, v.t5, 0.5, tension)
        Vec3.spline(t45, v.t3, v.t4, v.t5, v.t6, 0.5, tension)

        // const dp12 = Vec3.zero()
        // const dp23 = Vec3.zero()
        // const dp34 = Vec3.zero()
        // const dp45 = Vec3.zero()
        // Vec3.projectPointOnVector(dp12, v.d12, t12, v.t1)
        // Vec3.projectPointOnVector(dp23, v.d23, t23, v.t2)
        // Vec3.projectPointOnVector(dp34, v.d34, t34, v.t3)
        // Vec3.projectPointOnVector(dp45, v.d45, t45, v.t4)

        const td12 = Vec3.zero()
        const td23 = Vec3.zero()
        const td34 = Vec3.zero()
        const td45 = Vec3.zero()
        Vec3.normalize(td12, Vec3.sub(td12, t12, v.d12))
        Vec3.scaleAndAdd(v.d12, t12, td12, 1)
        Vec3.normalize(td23, Vec3.sub(td23, t23, v.d23))
        if (Vec3.dot(td12, td23) < 0) {
            Vec3.scaleAndAdd(v.d23, t23, td23, -1)
            console.log('foo td0 td1')
        } else {
            Vec3.scaleAndAdd(v.d23, t23, td23, 1)
        }
        Vec3.normalize(td34, Vec3.sub(td34, t34, v.d34))
        if (Vec3.dot(td12, td34) < 0) {
            Vec3.scaleAndAdd(v.d34, t34, td34, -1)
            console.log('foo td1 td2')
        } else {
            Vec3.scaleAndAdd(v.d34, t34, td34, 1)
        }
        Vec3.normalize(td45, Vec3.sub(td45, t45, v.d45))
        if (Vec3.dot(td12, td45) < 0) {
            Vec3.scaleAndAdd(v.d45, t45, td45, -1)
            console.log('foo td2 td3')
        } else {
            Vec3.scaleAndAdd(v.d45, t45, td45, 1)
        }

        // console.log(td0, td1, td2, td3)

        builder.addIcosahedron(t12, 0.3, 1)
        builder.addIcosahedron(t23, 0.3, 1)
        builder.addIcosahedron(t34, 0.3, 1)
        builder.addIcosahedron(t45, 0.3, 1)

        // builder.addIcosahedron(dp12, 0.3, 1)
        // builder.addIcosahedron(dp23, 0.3, 1)
        // builder.addIcosahedron(dp34, 0.3, 1)
        // builder.addIcosahedron(dp45, 0.3, 1)

        builder.addIcosahedron(v.d12, 0.3, 1)
        builder.addIcosahedron(v.d23, 0.3, 1)
        builder.addIcosahedron(v.d34, 0.3, 1)
        builder.addIcosahedron(v.d45, 0.3, 1)

        for (let j = 1; j <= linearSegments; ++j) {
            const t = j * 1.0 / linearSegments;
            Vec3.copy(tA, tB)
            // if ((v.last && t > 0.5) || (v.first && t < 0.5)) break

            if (t < 0.5) {
                Vec3.spline(tB, v.t1, v.t2, v.t3, v.t4, t + 0.5, tension)
            } else {
                Vec3.spline(tB, v.t2, v.t3, v.t4, v.t5, t - 0.5, tension)
            }
            Vec3.spline(dB, v.d12, v.d23, v.d34, v.d45, t, tension)

            // reflect(reflectedControlPoint, tB, tA, 1)
            Vec3.toArray(tB, controlPoints, j * 3)

            Vec3.normalize(torsionVec, Vec3.sub(torsionVec, tB, dB))
            // if (Vec3.dot(initialTorsionVec, torsionVec) < 0) Vec3.scale(torsionVec, torsionVec, -1)
            Vec3.toArray(torsionVec, torsionVectors, j * 3)

            // approximate tangent as direction to previous control point
            Vec3.normalize(tangentVec, Vec3.sub(tangentVec, tB, tA))
            Vec3.normalize(normalVec, Vec3.cross(normalVec, tangentVec, torsionVec))
            Vec3.toArray(normalVec, normalVectors, j * 3)

            // TODO size theme
            // builder.addCylinder(tA, tB, 1.0, { radiusTop: 0.3, radiusBottom: 0.3 })

            builder.addIcosahedron(dB, 0.1, 1)

            builder.addCylinder(tB, Vec3.add(tmp, tB, torsionVec), 1.0, { radiusTop: 0.1, radiusBottom: 0.1 })
            // builder.addCylinder(tB, Vec3.add(tmp, tB, normalVec), 1.0, { radiusTop: 0.1, radiusBottom: 0.1 })

            console.log(tA, tB)
        }

        builder.addTube(controlPoints, torsionVectors, normalVectors, linearSegments, radialSegments)

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Polymer trace mesh', current: i, max: polymerElementCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerTraceProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type PolymerTraceProps = Partial<typeof DefaultPolymerTraceProps>

export function PolymerTraceVisual(): UnitsVisual<PolymerTraceProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultPolymerTraceProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: PolymerTraceProps = {}) {
            currentProps = Object.assign({}, DefaultPolymerTraceProps, props)
            currentGroup = group

            const { colorTheme, unitKinds } = { ...DefaultPolymerTraceProps, ...props }
            const instanceCount = group.units.length
            const elementCount = group.elements.length
            const unit = group.units[0]

            mesh = unitKinds.includes(unit.kind)
                ? await createPolymerTraceMesh(ctx, unit, mesh)
                : Mesh.createEmpty(mesh)

            const transforms = createTransforms(group)
            const color = createColors(group, elementCount, colorTheme)
            const marker = createMarkers(instanceCount * elementCount)

            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...color,
                ...marker,
                aTransform: transforms,
                elements: mesh.indexBuffer,
                ...createMeshValues(currentProps, counts),
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: PolymerTraceProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                mesh = await createPolymerTraceMesh(ctx, unit, mesh)
                ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                const elementCount = currentGroup.elements.length
                if (ctx.shouldUpdate) await ctx.update('Computing trace colors');
                createColors(currentGroup, elementCount, newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (renderObject.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const indices = OrderedSet.ofSingleton(elementId as Element.Index);
                return Element.Loci([{ unit, indices }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(renderObject.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
