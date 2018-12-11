/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { RepresentationProps, Visual, VisualContext } from '../representation';
import { StructureMeshParams, StructurePointsParams, StructureLinesParams, StructureDirectVolumeParams } from './representation';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MeshRenderObject, PointsRenderObject, LinesRenderObject, DirectVolumeRenderObject } from 'mol-gl/render-object';
import { createUnitsMeshRenderObject, createUnitsPointsRenderObject, createUnitsTransform, createUnitsLinesRenderObject, createUnitsDirectVolumeRenderObject, includesUnitKind } from './visual/util/common';
import { deepEqual, ValueCell } from 'mol-util';
import { Interval } from 'mol-data/int';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { RenderableValues } from 'mol-gl/renderable/schema';
import { Geometry } from 'mol-geo/geometry/geometry';
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
import { Theme, createEmptyTheme } from 'mol-theme/theme';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { UnitsParams } from './units-representation';
import { RenderableState } from 'mol-gl/renderable';
import { Mat4 } from 'mol-math/linear-algebra';
import { setTransform } from 'mol-geo/geometry/transform-data';

export type StructureGroup = { structure: Structure, group: Unit.SymmetryGroup }

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> { }

type UnitsRenderObject = MeshRenderObject | LinesRenderObject | PointsRenderObject | DirectVolumeRenderObject

interface UnitsVisualBuilder<P extends UnitsParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(group: Unit.SymmetryGroup): LocationIterator
    getLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number): Loci
    mark(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme): void
}

interface UnitsVisualGeometryBuilder<P extends UnitsParams, G extends Geometry> extends UnitsVisualBuilder<P, G> {
    createEmptyGeometry(geometry?: G): G
    createRenderObject(group: Unit.SymmetryGroup, geometry: Geometry, locationIt: LocationIterator, theme: Theme, currentProps: PD.Values<P>): UnitsRenderObject
    updateValues(values: RenderableValues, newProps: Partial<PD.Values<P>>): void
    updateBoundingSphere(values: RenderableValues, geometry: Geometry): void
    updateRenderableState(state: RenderableState, props: Partial<PD.Values<P>>): void
}

export function UnitsVisual<P extends UnitsParams>(builder: UnitsVisualGeometryBuilder<P, Geometry>): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const { createEmptyGeometry, createRenderObject, updateValues, updateBoundingSphere, updateRenderableState } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: UnitsRenderObject | undefined

    let newProps: PD.Values<P> = Object.assign({}, defaultProps)
    let newTheme: Theme = createEmptyTheme()
    let newStructureGroup: StructureGroup

    let currentProps: PD.Values<P>
    let currentTheme: Theme
    let currentStructureGroup: StructureGroup

    let geometry: Geometry
    let locationIt: LocationIterator

    function prepareUpdate(theme: Theme, props: Partial<PD.Values<P>> = {}, structureGroup: StructureGroup) {
        if (!structureGroup && !currentStructureGroup) {
            throw new Error('missing structureGroup')
        }

        newProps = Object.assign({}, currentProps, props)
        newTheme = theme
        newStructureGroup = structureGroup

        VisualUpdateState.reset(updateState)

        if (!renderObject) {
            updateState.createNew = true
        } else if (!currentStructureGroup || newStructureGroup.group.hashCode !== currentStructureGroup.group.hashCode) {
            updateState.createNew = true
        }

        if (updateState.createNew) {
            updateState.createGeometry = true
            return
        }

        setUpdateState(updateState, newProps, currentProps, theme, currentTheme)

        if (!ColorTheme.areEqual(theme.color, currentTheme.color)) {
            // console.log('new colorTheme')
            updateState.updateColor = true
        }
        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) {
            // console.log('new unitKinds')
            updateState.createGeometry = true
        }

        if (newStructureGroup.group.transformHash !== currentStructureGroup.group.transformHash) {
            // console.log('new transformHash')
            if (newStructureGroup.group.units.length !== currentStructureGroup.group.units.length || updateState.updateColor) {
                updateState.updateTransform = true
            } else {
                updateState.updateMatrix = true
            }
        }

        // check if the conformation of unit.model has changed
        if (Unit.conformationId(newStructureGroup.group.units[0]) !== Unit.conformationId(currentStructureGroup.group.units[0])) {
            // console.log('new conformation')
            updateState.createGeometry = true
        }

        if (updateState.updateTransform) {
            updateState.updateColor = true
            updateState.updateMatrix = true
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true
        }
    }

    function update(newGeometry?: Geometry) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newStructureGroup.group)
            if (newGeometry) {
                renderObject = createRenderObject(newStructureGroup.group, newGeometry, locationIt, newTheme, newProps)
            } else {
                throw new Error('expected geometry to be given')
            }
        } else {
            if (!renderObject) {
                throw new Error('expected renderObject to be available')
            }

            locationIt.reset()

            if (updateState.updateTransform) {
                // console.log('update transform')
                locationIt = createLocationIterator(newStructureGroup.group)
                const { instanceCount, groupCount } = locationIt
                createMarkers(instanceCount * groupCount, renderObject.values)
            }

            if (updateState.updateMatrix) {
                // console.log('update matrix')
                createUnitsTransform(newStructureGroup.group, renderObject.values)
            }

            if (updateState.createGeometry) {
                // console.log('update geometry')
                if (newGeometry) {
                    ValueCell.update(renderObject.values.drawCount, Geometry.getDrawCount(newGeometry))
                } else {
                    throw new Error('expected geometry to be given')
                }
            }

            if (updateState.updateTransform || updateState.createGeometry) {
                updateBoundingSphere(renderObject.values, newGeometry || geometry)
            }

            if (updateState.updateSize) {
                // not all geometries have size data, so check here
                if ('uSize' in renderObject.values) {
                    // console.log('update size')
                    createSizes(locationIt, newTheme.size, renderObject.values)
                }
            }

            if (updateState.updateColor) {
                // console.log('update color')
                createColors(locationIt, newTheme.color, renderObject.values)
            }

            updateValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)
        }

        currentProps = newProps
        currentTheme = newTheme
        currentStructureGroup = newStructureGroup
        if (newGeometry) geometry = newGeometry
    }

    function _createGeometry(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: Geometry) {
        return includesUnitKind(props.unitKinds, unit)
                ? createGeometry(ctx, unit, structure, theme, props, geometry)
                : createEmptyGeometry(geometry)
    }

    return {
        get groupCount() { return locationIt ? locationIt.count : 0 },
        get renderObject () { return locationIt && locationIt.count ? renderObject : undefined },
        createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, structureGroup?: StructureGroup) {
            prepareUpdate(theme, props, structureGroup || currentStructureGroup)
            if (updateState.createGeometry) {
                const newGeometry = _createGeometry(ctx, newStructureGroup.group.units[0], newStructureGroup.structure, newTheme, newProps, geometry)
                return newGeometry instanceof Promise ? newGeometry.then(update) : update(newGeometry)
            } else {
                update()
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentStructureGroup, renderObject.id) : EmptyLoci
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
            if (isEveryLoci(loci) || (Structure.isLoci(loci) && loci.structure === currentStructureGroup.structure)) {
                changed = apply(Interval.ofBounds(0, groupCount * instanceCount))
            } else {
                changed = mark(loci, currentStructureGroup, apply)
            }
            if (changed) {
                ValueCell.update(tMarker, tMarker.ref.value)
            }
            return changed
        },
        setVisibility(value: boolean) {
            if (renderObject) renderObject.state.visible = value
        },
        setPickable(value: boolean) {
            if (renderObject) renderObject.state.pickable = value
        },
        setTransform(value: Mat4) {
            if (renderObject) setTransform(value, renderObject.values)
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
export type UnitsMeshParams = typeof UnitsMeshParams
export interface UnitsMeshVisualBuilder<P extends UnitsMeshParams> extends UnitsVisualBuilder<P, Mesh> { }

export function UnitsMeshVisual<P extends UnitsMeshParams>(builder: UnitsMeshVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<StructureMeshParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        createEmptyGeometry: Mesh.createEmpty,
        createRenderObject: createUnitsMeshRenderObject,
        updateValues: Mesh.updateValues,
        updateBoundingSphere: Mesh.updateBoundingSphere,
        updateRenderableState: Geometry.updateRenderableState
    })
}

// points

export const UnitsPointsParams = {
    ...StructurePointsParams,
    ...UnitsParams,
}
export type UnitsPointsParams = typeof UnitsPointsParams
export interface UnitsPointVisualBuilder<P extends UnitsPointsParams> extends UnitsVisualBuilder<P, Points> { }

export function UnitsPointsVisual<P extends UnitsPointsParams>(builder: UnitsPointVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<StructurePointsParams & UnitsParams>({
        ...builder,
        createEmptyGeometry: Points.createEmpty,
        createRenderObject: createUnitsPointsRenderObject,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true
        },
        updateValues: Points.updateValues,
        updateBoundingSphere: Points.updateBoundingSphere,
        updateRenderableState: Points.updateRenderableState
    })
}

// lines

export const UnitsLinesParams = {
    ...StructureLinesParams,
    ...UnitsParams,
}
export type UnitsLinesParams = typeof UnitsLinesParams
export interface UnitsLinesVisualBuilder<P extends UnitsLinesParams> extends UnitsVisualBuilder<P, Lines> { }

export function UnitsLinesVisual<P extends UnitsLinesParams>(builder: UnitsLinesVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<StructureLinesParams & UnitsParams>({
        ...builder,
        createEmptyGeometry: Lines.createEmpty,
        createRenderObject: createUnitsLinesRenderObject,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true
        },
        updateValues: Lines.updateValues,
        updateBoundingSphere: Lines.updateBoundingSphere,
        updateRenderableState: Geometry.updateRenderableState
    })
}

// direct-volume

export const UnitsDirectVolumeParams = {
    ...StructureDirectVolumeParams,
    ...UnitsParams,
}
export type UnitsDirectVolumeParams = typeof UnitsDirectVolumeParams
export interface UnitsDirectVolumeVisualBuilder<P extends UnitsDirectVolumeParams> extends UnitsVisualBuilder<P, DirectVolume> { }

export function UnitsDirectVolumeVisual<P extends UnitsDirectVolumeParams>(builder: UnitsDirectVolumeVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<StructureDirectVolumeParams & UnitsParams>({
        ...builder,
        createEmptyGeometry: DirectVolume.createEmpty,
        createRenderObject: createUnitsDirectVolumeRenderObject,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        updateValues: DirectVolume.updateValues,
        updateBoundingSphere: DirectVolume.updateBoundingSphere,
        updateRenderableState: DirectVolume.updateRenderableState
    })
}