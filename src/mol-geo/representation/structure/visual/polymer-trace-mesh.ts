/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '..';
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
import { SecondaryStructureType } from 'mol-model/structure/model/types';
// import { radToDeg } from 'mol-math/misc';

const tmpNormal = Vec3.zero()
const tangentVec = Vec3.zero()
const normalVec = Vec3.zero()
const binormalVec = Vec3.zero()

const prevNormal = Vec3.zero()

const orthogonalizeTmpVec = Vec3.zero()
/** Get a vector that is similar to b but orthogonal to a */
function orthogonalize(out: Vec3, a: Vec3, b: Vec3) {
    Vec3.normalize(orthogonalizeTmpVec, Vec3.cross(orthogonalizeTmpVec, a, b))
    Vec3.normalize(out, Vec3.cross(out, orthogonalizeTmpVec, a))
    return out
}

function interpolateNormals(controlPoints: Helpers.NumberArray, tangentVectors: Helpers.NumberArray, normalVectors: Helpers.NumberArray, binormalVectors: Helpers.NumberArray, firstNormalVector: Vec3, lastNormalVector: Vec3) {
    const n = controlPoints.length / 3
    // const n1 = n - 1

    // const angle = radToDeg(Math.acos(Vec3.dot(firstNormalVector, lastNormalVector)))
    // console.log('angle', angle)
    if (Vec3.dot(firstNormalVector, lastNormalVector) < 0) {
        Vec3.scale(lastNormalVector, lastNormalVector, -1)
        // console.log('flipped last normal vector')
    }

    Vec3.copy(prevNormal, firstNormalVector)

    for (let i = 0; i < n; ++i) {
        const t = i === 0 ? 0 : 1 / (n - i)
        Vec3.normalize(tmpNormal, Vec3.slerp(tmpNormal, prevNormal, lastNormalVector, t))

        Vec3.fromArray(tangentVec, tangentVectors, i * 3)

        orthogonalize(normalVec, tangentVec, tmpNormal)
        Vec3.toArray(normalVec, normalVectors, i * 3)

        // const deltaAngle = radToDeg(Math.acos(Vec3.dot(prevNormal, normalVec)))
        // if (deltaAngle > (angle / n1) * 5 && deltaAngle > 20) {
        //     console.warn(i, 'large delta angle', deltaAngle)
        // }
        // if (Vec3.dot(normalVec, prevNormal) < 0) {
        //     console.warn(i, 'flip compared to prev', radToDeg(Math.acos(Vec3.dot(prevNormal, normalVec))))
        // }
        Vec3.copy(prevNormal, normalVec)

        Vec3.normalize(binormalVec, Vec3.cross(binormalVec, tangentVec, normalVec))
        Vec3.toArray(binormalVec, binormalVectors, i * 3)
    }
}

async function createPolymerTraceMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    console.log('polymerElementCount trace', polymerElementCount)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)
    const linearSegments = 8
    const radialSegments = 12
    const tension = 0.9

    const tanA = Vec3.zero()
    const tanB = Vec3.zero()

    const tB = Vec3.zero()
    const tangentVec = Vec3.zero()

    const tmp = Vec3.zero()

    const pn = (linearSegments + 1) * 3
    const controlPoints = new Float32Array(pn)
    const tangentVectors = new Float32Array(pn)
    const normalVectors = new Float32Array(pn)
    const binormalVectors = new Float32Array(pn)

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()
        // builder.setId(elements[v.center.element])
        builder.setId(v.center.element)

        for (let j = 0; j <= linearSegments; ++j) {
            const t = j * 1.0 / linearSegments;
            // if ((v.last && t > 0.5) || (v.first && t < 0.5)) break

            if (t < 0.5) {
                Vec3.spline(tB, v.t0, v.t1, v.t2, v.t3, t + 0.5, tension)
                Vec3.spline(tanA, v.t0, v.t1, v.t2, v.t3, t + 0.5 + 0.01, tension)
                Vec3.spline(tanB, v.t0, v.t1, v.t2, v.t3, t + 0.5 - 0.01, tension)
            } else {
                Vec3.spline(tB, v.t1, v.t2, v.t3, v.t4, t - 0.5, tension)
                Vec3.spline(tanA, v.t1, v.t2, v.t3, v.t4, t - 0.5 + 0.01, tension)
                Vec3.spline(tanB, v.t1, v.t2, v.t3, v.t4, t - 0.5 - 0.01, tension)
            }
            Vec3.toArray(tB, controlPoints, j * 3)
            Vec3.normalize(tangentVec, Vec3.sub(tangentVec, tanA, tanB))
            Vec3.toArray(tangentVec, tangentVectors, j * 3)
        }

        const firstControlPoint = Vec3.zero()
        const lastControlPoint = Vec3.zero()
        const firstTangentVec = Vec3.zero()
        const lastTangentVec = Vec3.zero()
        const firstNormalVec = Vec3.zero()
        const lastNormalVec = Vec3.zero()
        const firstDirPoint = Vec3.zero()
        const lastDirPoint = Vec3.zero()

        Vec3.fromArray(firstControlPoint, controlPoints, 0)
        Vec3.fromArray(lastControlPoint, controlPoints, linearSegments * 3)
        Vec3.fromArray(firstTangentVec, tangentVectors, 0)
        Vec3.fromArray(lastTangentVec, tangentVectors, linearSegments * 3)
        Vec3.copy(firstDirPoint, v.d12)
        Vec3.copy(lastDirPoint, v.d23)

        Vec3.normalize(tmpNormal, Vec3.sub(tmp, firstControlPoint, firstDirPoint))
        orthogonalize(firstNormalVec, firstTangentVec, tmpNormal)

        Vec3.normalize(tmpNormal, Vec3.sub(tmp, lastControlPoint, lastDirPoint))
        orthogonalize(lastNormalVec, lastTangentVec, tmpNormal)

        // console.log('ELEMENT', i)
        interpolateNormals(controlPoints, tangentVectors, normalVectors, binormalVectors, firstNormalVec, lastNormalVec)

        // const controlPoint = Vec3.zero()
        // for (let j = 0; j <= linearSegments; ++j) {
        //     Vec3.fromArray(controlPoint, controlPoints, j * 3)
        //     Vec3.fromArray(normalVec, normalVectors, j * 3)
        //     Vec3.fromArray(binormalVec, binormalVectors, j * 3)
        //     Vec3.fromArray(tangentVec, tangentVectors, j * 3)
        //     builder.addIcosahedron(controlPoint, 0.25, 1)
        //     builder.addCylinder(
        //         controlPoint,
        //         Vec3.add(tmp, controlPoint, normalVec),
        //         1.5,
        //         { radiusTop: 0.07, radiusBottom: 0.07 }
        //     )
        //     builder.addCylinder(
        //         controlPoint,
        //         Vec3.add(tmp, controlPoint, binormalVec),
        //         0.8,
        //         { radiusTop: 0.07, radiusBottom: 0.07 }
        //     )
        //     builder.addCylinder(
        //         controlPoint,
        //         Vec3.add(tmp, controlPoint, tangentVec),
        //         j === 0 ? 2 : 1.5,
        //         { radiusTop: 0.03, radiusBottom: 0.03 }
        //     )
        // }

        // builder.addIcosahedron(v.t0, 0.25, 1)
        // builder.addIcosahedron(v.t1, 0.25, 1)
        // builder.addIcosahedron(v.t2, 0.25, 1)
        // builder.addIcosahedron(v.t3, 0.25, 1)
        // builder.addIcosahedron(v.t4, 0.25, 1)

        let width = 0.2, height = 0.2

        // TODO size theme
        if (SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Beta)) {
            width = 0.15; height = 1.0
            const arrowHeight = v.secStrucChange ? 1.7 : 0
            builder.addSheet(controlPoints, normalVectors, binormalVectors, linearSegments, width, height, arrowHeight, true, true)
        } else {
            if (SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Helix)) {
                width = 0.2; height = 1.0
            }
            builder.addTube(controlPoints, normalVectors, binormalVectors, linearSegments, radialSegments, width, height, 1, true, true)
        }

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
                const indices = OrderedSet.ofSingleton(elementId as StructureElement.UnitIndex);
                return StructureElement.Loci([{ unit, indices }])
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
