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
import { getPolymerElementCount, PolymerTraceIterator, interpolateNormals } from './util/polymer';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { SecondaryStructureType, MoleculeType } from 'mol-model/structure/model/types';

// TODO handle polymer ends properly
// TODO avoid allocating Vec3, use global temp vars

const t = Mat4.identity()

async function createPolymerTraceMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    console.log('polymerElementCount trace', polymerElementCount)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)
    const linearSegments = 8
    const radialSegments = 12

    const tanA = Vec3.zero()
    const tanB = Vec3.zero()

    const tB = Vec3.zero()
    const tangentVec = Vec3.zero()

    const pn = (linearSegments + 1) * 3
    const controlPoints = new Float32Array(pn)
    const tangentVectors = new Float32Array(pn)
    const normalVectors = new Float32Array(pn)
    const binormalVectors = new Float32Array(pn)

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()
        builder.setId(v.center.element)

        const isNucleic = v.moleculeType === MoleculeType.DNA || v.moleculeType === MoleculeType.RNA
        const isSheet = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Beta)
        const isHelix = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Helix)
        const tension = (isNucleic || isSheet) ? 0.5 : 0.9

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

        // console.log('ELEMENT', i)
        interpolateNormals(controlPoints, tangentVectors, normalVectors, binormalVectors, v.d12, v.d23)

        let width = 0.2, height = 0.2

        // TODO size theme
        if (isSheet) {
            width = 0.15; height = 1.0
            const arrowHeight = v.secStrucChange ? 1.7 : 0
            builder.addSheet(controlPoints, normalVectors, binormalVectors, linearSegments, width, height, arrowHeight, true, true)
        } else {
            if (isHelix) {
                width = 0.2; height = 1.0
            } else if (isNucleic) {
                width = 1.5; height = 0.3
            }
            builder.addTube(controlPoints, normalVectors, binormalVectors, linearSegments, radialSegments, width, height, 1, true, true)
        }

        if ((isSheet && !v.secStrucChange) || !isSheet) {
            const upVec = Vec3.zero()
            let width = 0.5, height = 1.2, depth = 0.6
            if (isNucleic) {
                Vec3.fromArray(upVec, binormalVectors, Math.round(linearSegments / 2) * 3)
                depth = 0.9
            } else {
                Vec3.fromArray(upVec, normalVectors, Math.round(linearSegments / 2) * 3)
            }

            Mat4.targetTo(t, v.t3, v.t1, upVec)
            Mat4.mul(t, t, Mat4.rotXY90)
            Mat4.setTranslation(t, v.t2)
            builder.addWedge(t, { width, height, depth })
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
