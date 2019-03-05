/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { RepresentationProps } from '../representation';
import { Visual, VisualContext } from '../visual';
import { StructureMeshParams, StructurePointsParams, StructureLinesParams, StructureDirectVolumeParams, StructureSpheresParams } from './representation';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { GraphicsRenderObject, createRenderObject } from 'mol-gl/render-object';
import { deepEqual, ValueCell } from 'mol-util';
import { Interval } from 'mol-data/int';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Geometry, GeometryUtils } from 'mol-geo/geometry/geometry';
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
import { Mat4 } from 'mol-math/linear-algebra';
import { Spheres } from 'mol-geo/geometry/spheres/spheres';
import { createUnitsTransform, includesUnitKind } from './visual/util/common';

export type StructureGroup = { structure: Structure, group: Unit.SymmetryGroup }

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> { }

function createUnitsRenderObject<G extends Geometry>(group: Unit.SymmetryGroup, geometry: G, locationIt: LocationIterator, theme: Theme, props: PD.Values<Geometry.Params<G>>) {
    const { createValues, createRenderableState } = Geometry.getUtils(geometry)
    const transform = createUnitsTransform(group)
    const values = createValues(geometry, transform, locationIt, theme, props)
    const state = createRenderableState(props)
    return createRenderObject(geometry.kind, values, state)
}

interface UnitsVisualBuilder<P extends UnitsParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(group: Unit.SymmetryGroup): LocationIterator
    getLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number): Loci
    mark(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme): void
}

interface UnitsVisualGeometryBuilder<P extends UnitsParams, G extends Geometry> extends UnitsVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function UnitsVisual<G extends Geometry, P extends UnitsParams & Geometry.Params<G>>(builder: UnitsVisualGeometryBuilder<P, G>): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const { createEmpty: createEmptyGeometry, updateValues, updateBoundingSphere, updateRenderableState } = builder.geometryUtils
    const updateState = VisualUpdateState.create()

    let renderObject: GraphicsRenderObject | undefined

    let newProps: PD.Values<P> = Object.assign({}, defaultProps)
    let newTheme: Theme = createEmptyTheme()
    let newStructureGroup: StructureGroup

    let currentProps: PD.Values<P>
    let currentTheme: Theme
    let currentStructureGroup: StructureGroup

    let geometry: G
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
        } else if (!currentStructureGroup || !Unit.SymmetryGroup.areInvariantElementsEqual(newStructureGroup.group, currentStructureGroup.group)) {
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
            updateState.updateSize = true
            updateState.updateMatrix = true
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true
            updateState.updateSize = true
        }
    }

    function update(newGeometry?: G) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newStructureGroup.group)
            if (newGeometry) {
                renderObject = createUnitsRenderObject(newStructureGroup.group, newGeometry, locationIt, newTheme, newProps)
            } else {
                throw new Error('expected geometry to be given')
            }
        } else {
            if (!renderObject) {
                throw new Error('expected renderObject to be available')
            }

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
                // console.log('UnitsVisual.updateBoundingSphere')
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

    function _createGeometry(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G) {
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
                return newGeometry instanceof Promise ? newGeometry.then(update) : update(newGeometry as G)
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
        setVisibility(visible: boolean) {
            Visual.setVisibility(renderObject, visible)
        },
        setPickable(pickable: boolean) {
            Visual.setPickable(renderObject, pickable)
        },
        setTransform(matrix?: Mat4, instanceMatrices?: Float32Array | null) {
            Visual.setTransform(renderObject, matrix, instanceMatrices)
        },
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}

// mesh

export const UnitsMeshParams = { ...StructureMeshParams, ...UnitsParams }
export type UnitsMeshParams = typeof UnitsMeshParams
export interface UnitsMeshVisualBuilder<P extends UnitsMeshParams> extends UnitsVisualBuilder<P, Mesh> { }

export function UnitsMeshVisual<P extends UnitsMeshParams>(builder: UnitsMeshVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<Mesh, StructureMeshParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        geometryUtils: Mesh.Utils
    })
}

// spheres

export const UnitsSpheresParams = { ...StructureSpheresParams, ...UnitsParams }
export type UnitsSpheresParams = typeof UnitsSpheresParams
export interface UnitsSpheresVisualBuilder<P extends UnitsSpheresParams> extends UnitsVisualBuilder<P, Spheres> { }

export function UnitsSpheresVisual<P extends UnitsSpheresParams>(builder: UnitsSpheresVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<Spheres, StructureSpheresParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true
        },
        geometryUtils: Spheres.Utils
    })
}

// points

export const UnitsPointsParams = { ...StructurePointsParams, ...UnitsParams }
export type UnitsPointsParams = typeof UnitsPointsParams
export interface UnitsPointVisualBuilder<P extends UnitsPointsParams> extends UnitsVisualBuilder<P, Points> { }

export function UnitsPointsVisual<P extends UnitsPointsParams>(builder: UnitsPointVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<Points, StructurePointsParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true
        },
        geometryUtils: Points.Utils
    })
}

// lines

export const UnitsLinesParams = { ...StructureLinesParams, ...UnitsParams }
export type UnitsLinesParams = typeof UnitsLinesParams
export interface UnitsLinesVisualBuilder<P extends UnitsLinesParams> extends UnitsVisualBuilder<P, Lines> { }

export function UnitsLinesVisual<P extends UnitsLinesParams>(builder: UnitsLinesVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<Lines, StructureLinesParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true
        },
        geometryUtils: Lines.Utils
    })
}

// direct-volume

export const UnitsDirectVolumeParams = { ...StructureDirectVolumeParams, ...UnitsParams }
export type UnitsDirectVolumeParams = typeof UnitsDirectVolumeParams
export interface UnitsDirectVolumeVisualBuilder<P extends UnitsDirectVolumeParams> extends UnitsVisualGeometryBuilder<P, DirectVolume> { }

export function UnitsDirectVolumeVisual<P extends UnitsDirectVolumeParams>(builder: UnitsDirectVolumeVisualBuilder<P>): UnitsVisual<P> {
    return UnitsVisual<DirectVolume, StructureDirectVolumeParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        geometryUtils: DirectVolume.Utils
    })
}