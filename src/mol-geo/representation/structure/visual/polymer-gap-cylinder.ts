/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { createTransforms, createColors } from './util/common';
import { deepEqual } from 'mol-util';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { createMarkers, MarkerAction } from '../../../util/marker-data';
import { Loci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerGapCount, PolymerGapIterator } from './util/polymer';
import { getElementLoci, markElement } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';

async function createPolymerGapCylinderMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    const polymerGapCount = getPolymerGapCount(unit)
    if (!polymerGapCount) return Mesh.createEmpty(mesh)
    console.log('polymerGapCount', polymerGapCount)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerGapCount * 30, polymerGapCount * 30 / 2, mesh)

    const { elements } = unit
    const pos = unit.conformation.invariantPosition
    const pA = Vec3.zero()
    const pB = Vec3.zero()

    let i = 0
    const polymerGapIt = PolymerGapIterator(unit)
    while (polymerGapIt.hasNext) {
        // TODO size theme
        const { centerA, centerB } = polymerGapIt.move()
        if (centerA.element === centerB.element) {
            builder.setId(centerA.element)
            pos(elements[centerA.element], pA)
            builder.addSphere(pA, 0.6, 0)
        } else {
            pos(elements[centerA.element], pA)
            pos(elements[centerB.element], pB)
            builder.setId(centerA.element)
            builder.addFixedCountDashedCylinder(pA, pB, 0.5, 10, { radiusTop: 0.2, radiusBottom: 0.2 })
            builder.setId(centerB.element)
            builder.addFixedCountDashedCylinder(pB, pA, 0.5, 10, { radiusTop: 0.2, radiusBottom: 0.2 })
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Gap mesh', current: i, max: polymerGapCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerGapProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type PolymerGapProps = Partial<typeof DefaultPolymerGapProps>

export function PolymerGapVisual(): UnitsVisual<PolymerGapProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultPolymerGapProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: PolymerGapProps = {}) {
            currentProps = Object.assign({}, DefaultPolymerGapProps, props)
            currentGroup = group

            const { colorTheme, unitKinds } = { ...DefaultPolymerGapProps, ...props }
            const instanceCount = group.units.length
            const elementCount = group.elements.length
            const unit = group.units[0]

            mesh = unitKinds.includes(unit.kind)
                ? await createPolymerGapCylinderMesh(ctx, unit, mesh)
                : Mesh.createEmpty(mesh)
            // console.log(mesh)

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
                aColor: ValueCell.create(new Float32Array(mesh.vertexCount * 3))
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: PolymerGapProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                mesh = await createPolymerGapCylinderMesh(ctx, unit, mesh)
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
            return getElementLoci(renderObject.id, currentGroup, pickingId)
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(renderObject.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
