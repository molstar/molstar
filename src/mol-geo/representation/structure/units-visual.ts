/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// TODO refactor to make DRY

import { Unit, Structure } from 'mol-model/structure';
import { RepresentationProps, Visual } from '../';
import { VisualUpdateState, StructureMeshParams, StructurePointsParams, StructureLinesParams, StructureDirectVolumeParams, StructureProps } from '.';
import { RuntimeContext } from 'mol-task';
import { PickingId } from '../../geometry/picking';
import { LocationIterator } from '../../util/location-iterator';
import { Mesh } from '../../geometry/mesh/mesh';
import { MarkerAction, applyMarkerAction, createMarkers } from '../../geometry/marker-data';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MeshRenderObject, PointsRenderObject, LinesRenderObject, DirectVolume2dRenderObject, DirectVolume3dRenderObject } from 'mol-gl/render-object';
import { createUnitsMeshRenderObject, createUnitsPointsRenderObject, createUnitsTransform, createUnitsLinesRenderObject, createUnitsDirectVolumeRenderObject } from './visual/util/common';
import { deepEqual, ValueCell, UUID } from 'mol-util';
import { Interval } from 'mol-data/int';
import { Points } from '../../geometry/points/points';
import { updateRenderableState, Geometry } from '../../geometry/geometry';
import { createColors, ColorProps } from '../../geometry/color-data';
import { createSizes, SizeProps } from '../../geometry/size-data';
import { Lines } from '../../geometry/lines/lines';
import { MultiSelectParam, paramDefaultValues } from 'mol-view/parameter';
import { DirectVolume2d, DirectVolume3d } from '../../geometry/direct-volume/direct-volume';

export const UnitKindInfo = {
    'atomic': {},
    'spheres': {},
    'gaussians': {},
}
export type UnitKind = keyof typeof UnitKindInfo
export const UnitKindNames = Object.keys(UnitKindInfo)
export const UnitKindOptions = UnitKindNames.map(n => [n, n] as [UnitKind, string])

export type StructureGroup = { structure: Structure, group: Unit.SymmetryGroup }

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> { }

function sameGroupConformation(groupA: Unit.SymmetryGroup, groupB: Unit.SymmetryGroup) {
    return (
        groupA.units.length === groupB.units.length &&
        Unit.conformationId(groupA.units[0]) === Unit.conformationId(groupB.units[0])
    )
}

function includesUnitKind(unitKinds: UnitKind[], unit: Unit) {
    for (let i = 0, il = unitKinds.length; i < il; ++i) {
        if (Unit.isAtomic(unit) && unitKinds[i] === 'atomic') return true
        if (Unit.isSpheres(unit) && unitKinds[i] === 'spheres') return true
        if (Unit.isGaussians(unit) && unitKinds[i] === 'gaussians') return true
    }
    return false
}

function sizeChanged(oldProps: SizeProps, newProps: SizeProps) {
    return (
        oldProps.sizeTheme !== newProps.sizeTheme ||
        oldProps.sizeValue !== newProps.sizeValue ||
        oldProps.sizeFactor !== newProps.sizeFactor
    )
}

function colorChanged(oldProps: ColorProps, newProps: ColorProps) {
    return (
        oldProps.colorTheme !== newProps.colorTheme ||
        oldProps.colorValue !== newProps.colorValue
    )
}

const UnitsParams = {
    unitKinds: MultiSelectParam<UnitKind>('Unit Kind', '', ['atomic', 'spheres'], UnitKindOptions),
}

interface UnitsVisualBuilder<P extends StructureProps, G extends Geometry> {
    defaultProps: P
    createGeometry(ctx: RuntimeContext, unit: Unit, structure: Structure, props: P, geometry?: G): Promise<G>
    createLocationIterator(group: Unit.SymmetryGroup): LocationIterator
    getLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number): Loci
    mark(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: P, currentProps: P): void
}

// mesh

export const UnitsMeshParams = {
    ...StructureMeshParams,
    ...UnitsParams,
}
export const DefaultUnitsMeshProps = paramDefaultValues(UnitsMeshParams)
export type UnitsMeshProps = typeof DefaultUnitsMeshProps
export interface UnitsMeshVisualBuilder<P extends UnitsMeshProps> extends UnitsVisualBuilder<P, Mesh> { }

export function UnitsMeshVisual<P extends UnitsMeshProps>(builder: UnitsMeshVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: MeshRenderObject | undefined
    let currentProps: P
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup
    let currentStructure: Structure
    let locationIt: LocationIterator
    let currentConformationId: UUID

    async function create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: Partial<P> = {}) {
        currentProps = Object.assign({}, defaultProps, props, { structure: currentStructure })
        currentGroup = group

        const unit = group.units[0]
        currentConformationId = Unit.conformationId(unit)
        mesh = includesUnitKind(currentProps.unitKinds, unit)
            ? await createGeometry(ctx, unit, currentStructure, currentProps, mesh)
            : Mesh.createEmpty(mesh)

        // TODO create empty location iterator when not in unitKinds
        locationIt = createLocationIterator(group)
        renderObject = await createUnitsMeshRenderObject(ctx, group, mesh, locationIt, currentProps)
    }

    async function update(ctx: RuntimeContext, props: Partial<P> = {}) {
        if (!renderObject) return

        const newProps = Object.assign({}, currentProps, props, { structure: currentStructure })
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

        if (sizeChanged(currentProps, newProps)) updateState.createGeometry = true
        if (colorChanged(currentProps, newProps)) updateState.updateColor = true
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
            mesh = includesUnitKind(newProps.unitKinds, unit)
                ? await createGeometry(ctx, unit, currentStructure, newProps, mesh)
                : Mesh.createEmpty(mesh)
            ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
            updateState.updateColor = true
        }

        if (updateState.updateColor) {
            await createColors(ctx, locationIt, newProps, renderObject.values)
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

// points

export const UnitsPointsParams = {
    ...StructurePointsParams,
    ...UnitsParams,
}
export const DefaultUnitsPointsProps = paramDefaultValues(UnitsPointsParams)
export type UnitsPointsProps = typeof DefaultUnitsPointsProps
export interface UnitsPointVisualBuilder<P extends UnitsPointsProps> extends UnitsVisualBuilder<P, Points> { }

export function UnitsPointsVisual<P extends UnitsPointsProps>(builder: UnitsPointVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: PointsRenderObject | undefined
    let currentProps: P
    let points: Points
    let currentGroup: Unit.SymmetryGroup
    let currentStructure: Structure
    let locationIt: LocationIterator
    let currentConformationId: UUID

    async function create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: Partial<P> = {}) {
        currentProps = Object.assign({}, defaultProps, props, { structure: currentStructure })
        currentGroup = group

        const unit = group.units[0]
        currentConformationId = Unit.conformationId(unit)
        points = includesUnitKind(currentProps.unitKinds, unit)
            ? await createGeometry(ctx, unit, currentStructure, currentProps, points)
            : Points.createEmpty(points)

        // TODO create empty location iterator when not in unitKinds
        locationIt = createLocationIterator(group)
        renderObject = await createUnitsPointsRenderObject(ctx, group, points, locationIt, currentProps)
    }

    async function update(ctx: RuntimeContext, props: Partial<P> = {}) {
        if (!renderObject) return

        const newProps = Object.assign({}, currentProps, props, { structure: currentStructure })
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

        if (sizeChanged(currentProps, newProps)) updateState.updateSize = true
        if (colorChanged(currentProps, newProps)) updateState.updateColor = true
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
            points = includesUnitKind(newProps.unitKinds, unit)
                ? await createGeometry(ctx, unit, currentStructure, newProps, points)
                : Points.createEmpty(points)
            ValueCell.update(renderObject.values.drawCount, points.pointCount)
            updateState.updateColor = true
        }

        if (updateState.updateSize) {
            await createSizes(ctx, locationIt, newProps, renderObject.values)
        }

        if (updateState.updateColor) {
            await createColors(ctx, locationIt, newProps, renderObject.values)
        }

        // TODO why do I need to cast here?
        Points.updateValues(renderObject.values, newProps as UnitsPointsProps)
        updateRenderableState(renderObject.state, newProps as UnitsPointsProps)

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

// lines

export const UnitsLinesParams = {
    ...StructureLinesParams,
    ...UnitsParams,
}
export const DefaultUnitsLinesProps = paramDefaultValues(UnitsLinesParams)
export type UnitsLinesProps = typeof DefaultUnitsLinesProps
export interface UnitsLinesVisualBuilder<P extends UnitsLinesProps> extends UnitsVisualBuilder<P, Lines> { }

export function UnitsLinesVisual<P extends UnitsLinesProps>(builder: UnitsLinesVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: LinesRenderObject | undefined
    let currentProps: P
    let lines: Lines
    let currentGroup: Unit.SymmetryGroup
    let currentStructure: Structure
    let locationIt: LocationIterator
    let currentConformationId: UUID

    async function create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: Partial<P> = {}) {
        currentProps = Object.assign({}, defaultProps, props, { structure: currentStructure })
        currentGroup = group

        const unit = group.units[0]
        currentConformationId = Unit.conformationId(unit)
        lines = includesUnitKind(currentProps.unitKinds, unit)
            ? await createGeometry(ctx, unit, currentStructure, currentProps, lines)
            : Lines.createEmpty(lines)

        // TODO create empty location iterator when not in unitKinds
        locationIt = createLocationIterator(group)
        renderObject = await createUnitsLinesRenderObject(ctx, group, lines, locationIt, currentProps)
    }

    async function update(ctx: RuntimeContext, props: Partial<P> = {}) {
        if (!renderObject) return

        const newProps = Object.assign({}, currentProps, props, { structure: currentStructure })
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

        if (sizeChanged(currentProps, newProps)) updateState.updateSize = true
        if (colorChanged(currentProps, newProps)) updateState.updateColor = true
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
            lines = includesUnitKind(newProps.unitKinds, unit)
                ? await createGeometry(ctx, unit, currentStructure, newProps, lines)
                : Lines.createEmpty(lines)
            ValueCell.update(renderObject.values.drawCount, lines.lineCount * 2 * 3)
            updateState.updateColor = true
        }

        if (updateState.updateSize) {
            await createSizes(ctx, locationIt, newProps, renderObject.values)
        }

        if (updateState.updateColor) {
            await createColors(ctx, locationIt, newProps, renderObject.values)
        }

        // TODO why do I need to cast here?
        Lines.updateValues(renderObject.values, newProps as UnitsLinesProps)
        updateRenderableState(renderObject.state, newProps as UnitsLinesProps)

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

// direct-volume

export const UnitsDirectVolumeParams = {
    ...StructureDirectVolumeParams,
    ...UnitsParams,
}
export const DefaultUnitsDirectVolumeProps = paramDefaultValues(UnitsDirectVolumeParams)
export type UnitsDirectVolumeProps = typeof DefaultUnitsDirectVolumeProps
export interface UnitsDirectVolumeVisualBuilder<P extends UnitsDirectVolumeProps> extends UnitsVisualBuilder<P, DirectVolume2d | DirectVolume3d> { }

export function UnitsDirectVolumeVisual<P extends UnitsDirectVolumeProps>(builder: UnitsDirectVolumeVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, setUpdateState } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: DirectVolume2dRenderObject | DirectVolume3dRenderObject | undefined
    let currentProps: P
    let directVolume: DirectVolume2d | DirectVolume3d
    let currentGroup: Unit.SymmetryGroup
    let currentStructure: Structure
    let locationIt: LocationIterator
    let currentConformationId: UUID

    async function create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: Partial<P> = {}) {
        const { webgl } = props
        if (webgl === undefined) throw new Error('UnitsDirectVolumeVisual requires `webgl` in props')

        currentProps = Object.assign({}, defaultProps, props, { structure: currentStructure })
        currentGroup = group

        const unit = group.units[0]
        currentConformationId = Unit.conformationId(unit)
        directVolume = includesUnitKind(currentProps.unitKinds, unit)
            ? await createGeometry(ctx, unit, currentStructure, currentProps, directVolume)
            : (webgl.isWebGL2 ? 
                DirectVolume2d.createEmpty(directVolume as DirectVolume2d) :
                DirectVolume3d.createEmpty(directVolume as DirectVolume3d))

        console.log('directVolume', directVolume)

        // TODO create empty location iterator when not in unitKinds
        locationIt = createLocationIterator(group)
        renderObject = await createUnitsDirectVolumeRenderObject(ctx, group, directVolume, locationIt, currentProps)
    }

    async function update(ctx: RuntimeContext, props: Partial<P> = {}) {
        const { webgl } = props
        if (webgl === undefined) throw new Error('UnitsDirectVolumeVisual requires `webgl` in props')

        if (!renderObject) return

        const newProps = Object.assign({}, currentProps, props, { structure: currentStructure })
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

        if (sizeChanged(currentProps, newProps)) updateState.createGeometry = true
        if (colorChanged(currentProps, newProps)) updateState.updateColor = true
        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) updateState.createGeometry = true

        //

        // if (updateState.updateTransform) {
        //     locationIt = createLocationIterator(currentGroup)
        //     const { instanceCount, groupCount } = locationIt
        //     createUnitsTransform(currentGroup, renderObject.values)
        //     createMarkers(instanceCount * groupCount, renderObject.values)
        //     updateState.updateColor = true
        // }

        if (updateState.createGeometry) {
            directVolume = includesUnitKind(newProps.unitKinds, unit)
                ? await createGeometry(ctx, unit, currentStructure, newProps, directVolume)
                : (webgl.isWebGL2 ? 
                    DirectVolume2d.createEmpty(directVolume as DirectVolume2d) :
                    DirectVolume3d.createEmpty(directVolume as DirectVolume3d))
            updateState.updateColor = true
        }

        // if (updateState.updateColor) {
        //     await createColors(ctx, locationIt, newProps, renderObject.values)
        // }

        if (renderObject.type === 'direct-volume-2d') {
            DirectVolume2d.updateValues(renderObject.values, newProps)
        } else {
            DirectVolume3d.updateValues(renderObject.values, newProps)
        }
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
            // TODO
            return false
        },
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}