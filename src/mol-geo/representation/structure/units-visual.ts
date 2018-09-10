/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { RepresentationProps, Visual } from '..';
import { DefaultStructureMeshProps, MeshUpdateState } from '.';
import { RuntimeContext } from 'mol-task';
import { PickingId } from '../../util/picking';
import { LocationIterator } from '../../util/location-iterator';
import { Mesh } from '../../mesh/mesh';
import { MarkerAction, applyMarkerAction, createMarkers } from '../../util/marker-data';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MeshRenderObject } from 'mol-gl/render-object';
import { createUnitsMeshRenderObject, createColors } from './visual/util/common';
import { deepEqual, ValueCell, UUID } from 'mol-util';
import { updateMeshValues, updateRenderableState } from '../util';
import { Interval } from 'mol-data/int';
import { createTransforms } from '../../util/transform-data';

export type StructureGroup = { structure: Structure, group: Unit.SymmetryGroup }

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> { }

export const DefaultUnitsMeshProps = {
    ...DefaultStructureMeshProps,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type UnitsMeshProps = typeof DefaultUnitsMeshProps

export interface UnitsMeshVisualBuilder<P extends UnitsMeshProps> {
    defaultProps: P
    createMesh(ctx: RuntimeContext, unit: Unit, structure: Structure, props: P, mesh?: Mesh): Promise<Mesh>
    createLocationIterator(group: Unit.SymmetryGroup): LocationIterator
    getLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number): Loci
    mark(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: MeshUpdateState, newProps: P, currentProps: P): void
}

export function UnitsMeshVisual<P extends UnitsMeshProps>(builder: UnitsMeshVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createMesh, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = MeshUpdateState.create()

    let renderObject: MeshRenderObject | undefined
    let currentProps: P
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup
    let currentStructure: Structure
    let locationIt: LocationIterator
    let currentConformationId: UUID

    async function create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: Partial<P> = {}) {
        currentProps = Object.assign({}, defaultProps, props)
        currentProps.colorTheme.structure = currentStructure
        currentGroup = group

        const unit = group.units[0]
        currentConformationId = Unit.conformationId(unit)
        mesh = currentProps.unitKinds.includes(unit.kind)
            ? await createMesh(ctx, unit, currentStructure, currentProps, mesh)
            : Mesh.createEmpty(mesh)

        // TODO create empty location iterator when not in unitKinds
        locationIt = createLocationIterator(group)
        renderObject = await createUnitsMeshRenderObject(ctx, group, mesh, locationIt, currentProps)
    }

    async function update(ctx: RuntimeContext, props: Partial<P> = {}) {
        if (!renderObject) return

        const newProps = Object.assign({}, currentProps, props)
        newProps.colorTheme.structure = currentStructure
        const unit = currentGroup.units[0]

        locationIt.reset()
        MeshUpdateState.reset(updateState)
        setUpdateState(updateState, newProps, currentProps)

        const newConformationId = Unit.conformationId(unit)
        if (newConformationId !== currentConformationId) {
            currentConformationId = newConformationId
            updateState.createMesh = true
        }

        if (currentGroup.units.length !== locationIt.instanceCount) updateState.updateTransform = true

        if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) updateState.createMesh = true
        if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) updateState.updateColor = true
        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) updateState.createMesh = true

        //

        if (updateState.updateTransform) {
            locationIt = createLocationIterator(currentGroup)
            const { instanceCount, groupCount } = locationIt
            createTransforms(currentGroup, renderObject.values)
            createMarkers(instanceCount * groupCount, renderObject.values)
            updateState.updateColor = true
        }

        if (updateState.createMesh) {
            mesh = newProps.unitKinds.includes(unit.kind)
                ? await createMesh(ctx, unit, currentStructure, newProps, mesh)
                : Mesh.createEmpty(mesh)
            ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
            updateState.updateColor = true
        }

        if (updateState.updateColor) {
            await createColors(ctx, locationIt, newProps.colorTheme, renderObject.values)
        }

        updateMeshValues(renderObject.values, newProps)
        updateRenderableState(renderObject.state, newProps)

        currentProps = newProps
    }

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: Partial<P> = {}, structureGroup?: StructureGroup) {
            if (structureGroup) currentStructure = structureGroup.structure
            const group = structureGroup ? structureGroup.group : undefined
            if (!group && !currentGroup) {
                throw new Error('missing group')
            } else if (group && (!currentGroup || !renderObject)) {
                await create(ctx, group, props)
            } else if (group && group.hashCode !== currentGroup.hashCode) {
                await create(ctx, group, props)
            } else {
                if (group && !areGroupsIdentical(group, currentGroup)) {
                    currentGroup = group
                }
                await update(ctx, props)
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentGroup, renderObject.id) : EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            if (!renderObject) return false
            const { tMarker } = renderObject.values
            const { groupCount, instanceCount } = locationIt

            function apply(interval: Interval) {
                const start = Interval.start(interval)
                const end = Interval.end(interval)
                return applyMarkerAction(tMarker.ref.value.array, start, end, action)
            }

            let changed = false
            if (isEveryLoci(loci)) {
                changed = apply(Interval.ofBounds(0, groupCount * instanceCount))
            } else {
                changed = mark(loci, currentGroup, apply)
            }
            if (changed) {
                ValueCell.update(tMarker, tMarker.ref.value)
            }
            return changed
        },
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}

function areGroupsIdentical(groupA: Unit.SymmetryGroup, groupB: Unit.SymmetryGroup) {
    return (
        groupA.units.length === groupB.units.length &&
        Unit.conformationId(groupA.units[0]) === Unit.conformationId(groupB.units[0])
    )
}