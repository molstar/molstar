/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'

import { MeshRenderObject } from 'mol-gl/render-object'
import { Unit } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { createColors, createUnitsMeshRenderObject } from './util/common';
import { createElementSphereMesh, markElement, getElementRadius, getElementLoci } from './util/element';
import { deepEqual } from 'mol-util';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { MarkerAction } from '../../../util/marker-data';
import { Loci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { updateMeshValues, updateRenderableState, DefaultMeshProps } from '../../util';
import { StructureElementIterator } from './util/location-iterator';

export const DefaultElementSphereProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type ElementSphereProps = Partial<typeof DefaultElementSphereProps>

export function ElementSphereVisual(): UnitsVisual<ElementSphereProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultElementSphereProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: ElementSphereProps = {}) {
            currentProps = Object.assign({}, DefaultElementSphereProps, props)
            currentGroup = group

            const { detail, sizeTheme, unitKinds } = { ...DefaultElementSphereProps, ...props }
            const unit = group.units[0]

            const radius = getElementRadius(unit, sizeTheme)
            mesh = unitKinds.includes(unit.kind)
                ? await createElementSphereMesh(ctx, unit, radius, detail, mesh)
                : Mesh.createEmpty(mesh)

            const locationIt = StructureElementIterator.fromGroup(group)
            renderObject = createUnitsMeshRenderObject(group, mesh, locationIt, currentProps)
        },
        async update(ctx: RuntimeContext, props: ElementSphereProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                const radius = getElementRadius(unit, newProps.sizeTheme)
                mesh = await createElementSphereMesh(ctx, unit, radius, newProps.detail, mesh)
                ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                createColors(StructureElementIterator.fromGroup(currentGroup), newProps.colorTheme, renderObject.values)
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
