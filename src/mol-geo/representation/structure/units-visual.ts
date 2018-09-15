/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { RepresentationProps, Visual } from '../';
import { DefaultStructureMeshProps, VisualUpdateState, DefaultStructurePointProps } from '.';
import { RuntimeContext } from 'mol-task';
import { PickingId } from '../../geometry/picking';
import { LocationIterator } from '../../util/location-iterator';
import { Mesh } from '../../geometry/mesh/mesh';
import { MarkerAction, applyMarkerAction, createMarkers } from '../../geometry/marker-data';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MeshRenderObject, PointRenderObject } from 'mol-gl/render-object';
import { createUnitsMeshRenderObject, createUnitsPointRenderObject, createUnitsTransform } from './visual/util/common';
import { deepEqual, ValueCell, UUID } from 'mol-util';
import { Interval } from 'mol-data/int';
import { Point } from '../../geometry/point/point';
import { updateRenderableState } from '../../geometry/geometry';
import { createColors } from '../../geometry/color-data';
import { createSizes } from '../../geometry/size-data';

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
    setUpdateState(state: VisualUpdateState, newProps: P, currentProps: P): void
}

export function UnitsMeshVisual<P extends UnitsMeshProps>(builder: UnitsMeshVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createMesh, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = VisualUpdateState.create()

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
        VisualUpdateState.reset(updateState)
        setUpdateState(updateState, newProps, currentProps)

        const newConformationId = Unit.conformationId(unit)
        if (newConformationId !== currentConformationId) {
            currentConformationId = newConformationId
            updateState.createGeometry = true
        }

        if (currentGroup.units.length !== locationIt.instanceCount) updateState.updateTransform = true

        if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) updateState.createGeometry = true
        if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) updateState.updateColor = true
        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) updateState.createGeometry = true

        //

        if (updateState.updateTransform) {
            locationIt = createLocationIterator(currentGroup)
            const { instanceCount, groupCount } = locationIt
            createUnitsTransform(currentGroup, renderObject.values)
            createMarkers(instanceCount * groupCount, renderObject.values)
            updateState.updateColor = true
        }

        if (updateState.createGeometry) {
            mesh = newProps.unitKinds.includes(unit.kind)
                ? await createMesh(ctx, unit, currentStructure, newProps, mesh)
                : Mesh.createEmpty(mesh)
            ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
            updateState.updateColor = true
        }

        if (updateState.updateColor) {
            await createColors(ctx, locationIt, newProps.colorTheme, renderObject.values)
        }

        // TODO why do I need to cast here?
        Mesh.updateValues(renderObject.values, newProps as UnitsMeshProps)
        updateRenderableState(renderObject.state, newProps as UnitsMeshProps)

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
                // console.log('unit-visual first create')
                await create(ctx, group, props)
            } else if (group && group.hashCode !== currentGroup.hashCode) {
                // console.log('unit-visual group.hashCode !== currentGroup.hashCode')
                await create(ctx, group, props)
            } else {
                // console.log('unit-visual update')
                if (group && !sameGroupConformation(group, currentGroup)) {
                    // console.log('unit-visual new conformation')
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

function sameGroupConformation(groupA: Unit.SymmetryGroup, groupB: Unit.SymmetryGroup) {
    return (
        groupA.units.length === groupB.units.length &&
        Unit.conformationId(groupA.units[0]) === Unit.conformationId(groupB.units[0])
    )
}

//

export const DefaultUnitsPointProps = {
    ...DefaultStructurePointProps,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type UnitsPointProps = typeof DefaultUnitsPointProps

export interface UnitsPointVisualBuilder<P extends UnitsPointProps> {
    defaultProps: P
    createPoint(ctx: RuntimeContext, unit: Unit, structure: Structure, props: P, point?: Point): Promise<Point>
    createLocationIterator(group: Unit.SymmetryGroup): LocationIterator
    getLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number): Loci
    mark(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: P, currentProps: P): void
}

export function UnitsPointVisual<P extends UnitsPointProps>(builder: UnitsPointVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createPoint, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: PointRenderObject | undefined
    let currentProps: P
    let point: Point
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
        point = currentProps.unitKinds.includes(unit.kind)
            ? await createPoint(ctx, unit, currentStructure, currentProps, point)
            : Point.createEmpty(point)

        // TODO create empty location iterator when not in unitKinds
        locationIt = createLocationIterator(group)
        renderObject = await createUnitsPointRenderObject(ctx, group, point, locationIt, currentProps)
    }

    async function update(ctx: RuntimeContext, props: Partial<P> = {}) {
        if (!renderObject) return

        const newProps = Object.assign({}, currentProps, props)
        newProps.colorTheme.structure = currentStructure
        const unit = currentGroup.units[0]

        locationIt.reset()
        VisualUpdateState.reset(updateState)
        setUpdateState(updateState, newProps, currentProps)

        const newConformationId = Unit.conformationId(unit)
        if (newConformationId !== currentConformationId) {
            currentConformationId = newConformationId
            updateState.createGeometry = true
        }

        if (currentGroup.units.length !== locationIt.instanceCount) updateState.updateTransform = true

        if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) updateState.updateSize = true
        if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) updateState.updateColor = true
        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) updateState.createGeometry = true

        //

        if (updateState.updateTransform) {
            locationIt = createLocationIterator(currentGroup)
            const { instanceCount, groupCount } = locationIt
            createUnitsTransform(currentGroup, renderObject.values)
            createMarkers(instanceCount * groupCount, renderObject.values)
            updateState.updateColor = true
        }

        if (updateState.createGeometry) {
            point = newProps.unitKinds.includes(unit.kind)
                ? await createPoint(ctx, unit, currentStructure, newProps, point)
                : Point.createEmpty(point)
            ValueCell.update(renderObject.values.drawCount, point.vertexCount)
            updateState.updateColor = true
        }

        if (updateState.updateSize) {
            await createSizes(ctx, locationIt, newProps.sizeTheme, renderObject.values)
        }

        if (updateState.updateColor) {
            await createColors(ctx, locationIt, newProps.colorTheme, renderObject.values)
        }

        // TODO why do I need to cast here?
        Point.updateValues(renderObject.values, newProps as UnitsPointProps)
        updateRenderableState(renderObject.state, newProps as UnitsPointProps)

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
                // console.log('unit-visual first create')
                await create(ctx, group, props)
            } else if (group && group.hashCode !== currentGroup.hashCode) {
                // console.log('unit-visual group.hashCode !== currentGroup.hashCode')
                await create(ctx, group, props)
            } else {
                // console.log('unit-visual update')
                if (group && !sameGroupConformation(group, currentGroup)) {
                    // console.log('unit-visual new conformation')
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