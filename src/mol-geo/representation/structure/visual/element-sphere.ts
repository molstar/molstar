/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Element } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '../index';
import { RuntimeContext } from 'mol-task'
import { createTransforms, createColors } from '../visual/util/common';
import { createElementSphereMesh, markElement, getElementRadius } from '../visual/util/element';
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

export const DefaultElementSphereProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type ElementSphereProps = Partial<typeof DefaultElementSphereProps>

export function ElementSphereVisual(): UnitsVisual<ElementSphereProps> {
    const renderObjects: RenderObject[] = []
    let spheres: MeshRenderObject
    let currentProps: typeof DefaultElementSphereProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        renderObjects,
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: ElementSphereProps = {}) {
            currentProps = Object.assign({}, DefaultElementSphereProps, props)

            renderObjects.length = 0 // clear
            currentGroup = group

            const { detail, colorTheme, sizeTheme, unitKinds } = { ...DefaultElementSphereProps, ...props }
            const instanceCount = group.units.length
            const elementCount = group.elements.length
            const unit = group.units[0]

            const radius = getElementRadius(unit, sizeTheme)
            mesh = unitKinds.includes(unit.kind)
                ? await createElementSphereMesh(ctx, unit, radius, detail, mesh)
                : Mesh.createEmpty(mesh)

            if (ctx.shouldUpdate) await ctx.update('Computing spacefill transforms');
            const transforms = createTransforms(group)

            if (ctx.shouldUpdate) await ctx.update('Computing spacefill colors');
            const color = createColors(group, elementCount, colorTheme)

            if (ctx.shouldUpdate) await ctx.update('Computing spacefill marks');
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

            spheres = createMeshRenderObject(values, state)
            renderObjects.push(spheres)
        },
        async update(ctx: RuntimeContext, props: ElementSphereProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!spheres) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                const radius = getElementRadius(unit, newProps.sizeTheme)
                mesh = await createElementSphereMesh(ctx, unit, radius, newProps.detail, mesh)
                ValueCell.update(spheres.values.drawCount, mesh.triangleCount * 3)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                const elementCount = currentGroup.elements.length
                if (ctx.shouldUpdate) await ctx.update('Computing spacefill colors');
                createColors(currentGroup, elementCount, newProps.colorTheme, spheres.values)
            }

            updateMeshValues(spheres.values, newProps)
            updateRenderableState(spheres.state, newProps)

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (spheres.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const indices = OrderedSet.ofSingleton(elementId as Element.Index);
                return Element.Loci([{ unit, indices }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(spheres.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
