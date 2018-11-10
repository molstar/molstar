/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { RepresentationProps, Visual, VisualContext } from '../representation';
import { StructureMeshParams, StructurePointsParams, StructureLinesParams, StructureDirectVolumeParams, StructureParams } from './representation';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MeshRenderObject, PointsRenderObject, LinesRenderObject, DirectVolumeRenderObject } from 'mol-gl/render-object';
import { createUnitsMeshRenderObject, createUnitsPointsRenderObject, createUnitsTransform, createUnitsLinesRenderObject, createUnitsDirectVolumeRenderObject, UnitKind, UnitKindOptions, includesUnitKind } from './visual/util/common';
import { deepEqual, ValueCell, UUID } from 'mol-util';
import { Interval } from 'mol-data/int';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { RenderableValues } from 'mol-gl/renderable/schema';
import { Geometry, updateRenderableState } from 'mol-geo/geometry/geometry';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { createMarkers, MarkerAction, applyMarkerAction } from 'mol-geo/geometry/marker-data';
import { createSizes } from 'mol-geo/geometry/size-data';
import { createColors } from 'mol-geo/geometry/color-data';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Points } from 'mol-geo/geometry/points/points';
import { Lines } from 'mol-geo/geometry/lines/lines';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';
import { VisualUpdateState } from 'mol-repr/util';
import { Theme } from 'mol-theme/theme';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';

export type StructureGroup = { structure: Structure, group: Unit.SymmetryGroup }

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> { }

function sameGroupConformation(groupA: Unit.SymmetryGroup, groupB: Unit.SymmetryGroup) {
    return (
        groupA.units.length === groupB.units.length &&
        Unit.conformationId(groupA.units[0]) === Unit.conformationId(groupB.units[0])
    )
}

const UnitsParams = {
    ...StructureParams,
    unitKinds: PD.MultiSelect<UnitKind>('Unit Kind', '', ['atomic', 'spheres'], UnitKindOptions),
}
const DefaultUnitsProps = PD.getDefaultValues(UnitsParams)
type UnitsProps = typeof DefaultUnitsProps

type UnitsRenderObject = MeshRenderObject | LinesRenderObject | PointsRenderObject | DirectVolumeRenderObject

interface UnitsVisualBuilder<P extends UnitsProps, G extends Geometry> {
    defaultProps: P
    createGeometry(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: P, geometry?: G): Promise<G>
    createLocationIterator(group: Unit.SymmetryGroup): LocationIterator
    getLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number): Loci
    mark(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: P, currentProps: P, newTheme: Theme, currentTheme: Theme): void
}

interface UnitsVisualGeometryBuilder<P extends UnitsProps, G extends Geometry> extends UnitsVisualBuilder<P, G> {
    createEmptyGeometry(geometry?: G): G
    createRenderObject(ctx: VisualContext, group: Unit.SymmetryGroup, geometry: Geometry, locationIt: LocationIterator, theme: Theme, currentProps: P): Promise<UnitsRenderObject>
    updateValues(values: RenderableValues, newProps: P): void
}

export function UnitsVisual<P extends UnitsProps>(builder: UnitsVisualGeometryBuilder<P, Geometry>): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const { createEmptyGeometry, createRenderObject, updateValues } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: UnitsRenderObject | undefined
    let currentProps: P
    let currentTheme: Theme
    let geometry: Geometry
    let currentGroup: Unit.SymmetryGroup
    let currentStructure: Structure
    let locationIt: LocationIterator
    let currentConformationId: UUID

    async function create(ctx: VisualContext, group: Unit.SymmetryGroup, theme: Theme, props: Partial<P> = {}) {
        currentProps = Object.assign({}, defaultProps, props, { structure: currentStructure })
        currentTheme = theme
        currentGroup = group

        const unit = group.units[0]
        currentConformationId = Unit.conformationId(unit)
        geometry = includesUnitKind(currentProps.unitKinds, unit)
            ? await createGeometry(ctx, unit, currentStructure, theme, currentProps, geometry)
            : createEmptyGeometry(geometry)

        // TODO create empty location iterator when not in unitKinds
        locationIt = createLocationIterator(group)
        renderObject = await createRenderObject(ctx, group, geometry, locationIt, theme, currentProps)
    }

    async function update(ctx: VisualContext, theme: Theme, props: Partial<P> = {}) {
        if (!renderObject) return

        const newProps = Object.assign({}, currentProps, props, { structure: currentStructure })
        const unit = currentGroup.units[0]

        locationIt.reset()
        VisualUpdateState.reset(updateState)
        setUpdateState(updateState, newProps, currentProps, theme, currentTheme)

        const newConformationId = Unit.conformationId(unit)
        if (newConformationId !== currentConformationId) {
            currentConformationId = newConformationId
            updateState.createGeometry = true
        }

        if (currentGroup.units.length !== locationIt.instanceCount) updateState.updateTransform = true

        if (ColorTheme.areEqual(theme.color, currentTheme.color)) updateState.updateColor = true
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
            geometry = includesUnitKind(newProps.unitKinds, unit)
                ? await createGeometry(ctx, unit, currentStructure, theme, newProps, geometry)
                : createEmptyGeometry(geometry)
            ValueCell.update(renderObject.values.drawCount, Geometry.getDrawCount(geometry))
            updateState.updateColor = true
        }

        if (updateState.updateSize) {
            // not all geometries have size data, so check here
            if ('uSize' in renderObject.values) {
                await createSizes(ctx.runtime, locationIt, theme.size, renderObject.values)
            }
        }

        if (updateState.updateColor) {
            await createColors(ctx.runtime, locationIt, theme.color, renderObject.values)
        }

        updateValues(renderObject.values, newProps)
        updateRenderableState(renderObject.state, newProps)

        currentProps = newProps
        currentTheme = theme
    }

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<P> = {}, structureGroup?: StructureGroup) {
            if (structureGroup) currentStructure = structureGroup.structure
            const group = structureGroup ? structureGroup.group : undefined
            if (!group && !currentGroup) {
                throw new Error('missing group')
            } else if (group && (!currentGroup || !renderObject)) {
                // console.log('unit-visual first create')
                await create(ctx, group, theme, props)
            } else if (group && group.hashCode !== currentGroup.hashCode) {
                // console.log('unit-visual group.hashCode !== currentGroup.hashCode')
                await create(ctx, group, theme, props)
            } else {
                // console.log('unit-visual update')
                if (group && !sameGroupConformation(group, currentGroup)) {
                    // console.log('unit-visual new conformation')
                    currentGroup = group
                }
                await update(ctx, theme, props)
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, { structure: currentStructure, group: currentGroup }, renderObject.id) : EmptyLoci
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
                changed = mark(loci, { structure: currentStructure, group: currentGroup }, apply)
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

// mesh

export const UnitsMeshParams = {
    ...StructureMeshParams,
    ...UnitsParams,
}
export const DefaultUnitsMeshProps = PD.getDefaultValues(UnitsMeshParams)
export type UnitsMeshProps = typeof DefaultUnitsMeshProps
export interface UnitsMeshVisualBuilder<P extends UnitsMeshProps> extends UnitsVisualBuilder<P, Mesh> { }

export function UnitsMeshVisual<P extends UnitsMeshProps>(builder: UnitsMeshVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: P, currentProps: P, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        createEmptyGeometry: Mesh.createEmpty,
        createRenderObject: createUnitsMeshRenderObject,
        updateValues: Mesh.updateValues
    })
}

// points

export const UnitsPointsParams = {
    ...StructurePointsParams,
    ...UnitsParams,
}
export const DefaultUnitsPointsProps = PD.getDefaultValues(UnitsPointsParams)
export type UnitsPointsProps = typeof DefaultUnitsPointsProps
export interface UnitsPointVisualBuilder<P extends UnitsPointsProps> extends UnitsVisualBuilder<P, Points> { }

export function UnitsPointsVisual<P extends UnitsPointsProps>(builder: UnitsPointVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual({
        ...builder,
        createEmptyGeometry: Points.createEmpty,
        createRenderObject: createUnitsPointsRenderObject,
        setUpdateState: (state: VisualUpdateState, newProps: P, currentProps: P, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true
        },
        updateValues: Points.updateValues
    })
}

// lines

export const UnitsLinesParams = {
    ...StructureLinesParams,
    ...UnitsParams,
}
export const DefaultUnitsLinesProps = PD.getDefaultValues(UnitsLinesParams)
export type UnitsLinesProps = typeof DefaultUnitsLinesProps
export interface UnitsLinesVisualBuilder<P extends UnitsLinesProps> extends UnitsVisualBuilder<P, Lines> { }

export function UnitsLinesVisual<P extends UnitsLinesProps>(builder: UnitsLinesVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual({
        ...builder,
        createEmptyGeometry: Lines.createEmpty,
        createRenderObject: createUnitsLinesRenderObject,
        setUpdateState: (state: VisualUpdateState, newProps: P, currentProps: P, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true
        },
        updateValues: Lines.updateValues
    })
}

// direct-volume

export const UnitsDirectVolumeParams = {
    ...StructureDirectVolumeParams,
    ...UnitsParams,
}
export const DefaultUnitsDirectVolumeProps = PD.getDefaultValues(UnitsDirectVolumeParams)
export type UnitsDirectVolumeProps = typeof DefaultUnitsDirectVolumeProps
export interface UnitsDirectVolumeVisualBuilder<P extends UnitsDirectVolumeProps> extends UnitsVisualBuilder<P, DirectVolume> { }

export function UnitsDirectVolumeVisual<P extends UnitsDirectVolumeProps>(builder: UnitsDirectVolumeVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual({
        ...builder,
        createEmptyGeometry: DirectVolume.createEmpty,
        createRenderObject: createUnitsDirectVolumeRenderObject,
        setUpdateState: (state: VisualUpdateState, newProps: P, currentProps: P, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        updateValues: DirectVolume.updateValues
    })
}