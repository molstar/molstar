/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { Visual, VisualContext } from '../representation';
import { MeshRenderObject, LinesRenderObject, PointsRenderObject, DirectVolumeRenderObject } from 'mol-gl/render-object';
import { createComplexMeshRenderObject, UnitKind, UnitKindOptions, createComplexDirectVolumeRenderObject } from './visual/util/common';
import { StructureMeshParams, StructureParams, StructureDirectVolumeParams } from './representation';
import { deepEqual, ValueCell } from 'mol-util';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { Interval } from 'mol-data/int';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { RenderableValues } from 'mol-gl/renderable/schema';
import { createSizes } from 'mol-geo/geometry/size-data';
import { Geometry } from 'mol-geo/geometry/geometry';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { createColors } from 'mol-geo/geometry/color-data';
import { MarkerAction, applyMarkerAction } from 'mol-geo/geometry/marker-data';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { VisualUpdateState } from 'mol-repr/util';
import { Theme, createEmptyTheme } from 'mol-theme/theme';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { RenderableState } from 'mol-gl/renderable';
import { UnitsParams } from './units-representation';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';
import { Mat4 } from 'mol-math/linear-algebra';
import { setTransform } from 'mol-geo/geometry/transform-data';

export interface  ComplexVisual<P extends StructureParams> extends Visual<Structure, P> { }

const ComplexParams = {
    ...StructureParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic', 'spheres'], UnitKindOptions),
}
type ComplexParams = typeof ComplexParams

type ComplexRenderObject = MeshRenderObject | LinesRenderObject | PointsRenderObject | DirectVolumeRenderObject

interface ComplexVisualBuilder<P extends ComplexParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(structure: Structure): LocationIterator
    getLoci(pickingId: PickingId, structure: Structure, id: number): Loci
    mark(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean): boolean,
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme): void
}

interface ComplexVisualGeometryBuilder<P extends ComplexParams, G extends Geometry> extends ComplexVisualBuilder<P, G> {
    createEmptyGeometry(geometry?: G): G
    createRenderObject(structure: Structure, geometry: Geometry, locationIt: LocationIterator, theme: Theme, currentProps: PD.Values<P>): ComplexRenderObject
    updateValues(values: RenderableValues, newProps: PD.Values<P>): void,
    updateBoundingSphere(values: RenderableValues, geometry: Geometry): void
    updateRenderableState(state: RenderableState, props: PD.Values<P>): void
}

export function ComplexVisual<P extends ComplexParams>(builder: ComplexVisualGeometryBuilder<P, Geometry>): ComplexVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const { createRenderObject, updateValues, updateBoundingSphere, updateRenderableState } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: ComplexRenderObject | undefined

    let newProps: PD.Values<P>
    let newTheme: Theme
    let newStructure: Structure

    let currentProps: PD.Values<P> = Object.assign({}, defaultProps)
    let currentTheme: Theme = createEmptyTheme()
    let currentStructure: Structure

    let geometry: Geometry
    let locationIt: LocationIterator

    function prepareUpdate(theme: Theme, props: Partial<PD.Values<P>>, structure: Structure) {
        if (!structure && !currentStructure) {
            throw new Error('missing structure')
        }

        newProps = Object.assign({}, currentProps, props)
        newTheme = theme
        newStructure = structure

        VisualUpdateState.reset(updateState)

        if (!renderObject) {
            updateState.createNew = true
        } else if (!currentStructure || !Structure.areEquivalent(newStructure, currentStructure)) {
            updateState.createNew = true
        }

        if (updateState.createNew) {
            updateState.createGeometry = true
            return
        }

        setUpdateState(updateState, newProps, currentProps, newTheme, currentTheme)

        if (Structure.conformationHash(newStructure) !== Structure.conformationHash(currentStructure)) {
            updateState.createGeometry = true
        }

        if (!ColorTheme.areEqual(theme.color, currentTheme.color)) updateState.updateColor = true
        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) updateState.createGeometry = true

        if (updateState.createGeometry) {
            updateState.updateColor = true
        }
    }

    function update(newGeometry?: Geometry) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newStructure)
            if (newGeometry) {
                renderObject = createRenderObject(newStructure, newGeometry, locationIt, newTheme, newProps)
            } else {
                throw new Error('expected geometry to be given')
            }
        } else {
            if (!renderObject) {
                throw new Error('expected renderObject to be available')
            }

            locationIt.reset()

            if (updateState.createGeometry) {
                if (newGeometry) {
                    ValueCell.update(renderObject.values.drawCount, Geometry.getDrawCount(newGeometry))
                    updateBoundingSphere(renderObject.values, newGeometry)
                } else {
                    throw new Error('expected geometry to be given')
                }
            }

            if (updateState.updateSize) {
                // not all geometries have size data, so check here
                if ('uSize' in renderObject.values) {
                    createSizes(locationIt, newTheme.size, renderObject.values)
                }
            }

            if (updateState.updateColor) {
                createColors(locationIt, newTheme.color, renderObject.values)
            }

            updateValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)
        }

        currentProps = newProps
        currentTheme = newTheme
        currentStructure = newStructure
        if (newGeometry) geometry = newGeometry
    }

    return {
        get groupCount() { return locationIt ? locationIt.count : 0 },
        get renderObject () { return locationIt && locationIt.count ? renderObject : undefined },
        createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, structure?: Structure) {
            prepareUpdate(theme, props, structure || currentStructure)
            if (updateState.createGeometry) {
                const newGeometry = createGeometry(ctx, newStructure, newTheme, newProps, geometry)
                return newGeometry instanceof Promise ? newGeometry.then(update) : update(newGeometry)
            } else {
                update()
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentStructure, renderObject.id) : EmptyLoci
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
            if (isEveryLoci(loci) || (Structure.isLoci(loci) && loci.structure === currentStructure)) {
                changed = apply(Interval.ofBounds(0, groupCount * instanceCount))
            } else {
                changed = mark(loci, currentStructure, apply)
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

export const ComplexMeshParams = {
    ...StructureMeshParams,
    unitKinds: PD.MultiSelect<UnitKind>([ 'atomic', 'spheres' ], UnitKindOptions),
}
export type ComplexMeshParams = typeof ComplexMeshParams

export interface ComplexMeshVisualBuilder<P extends ComplexMeshParams> extends ComplexVisualBuilder<P, Mesh> { }

export function ComplexMeshVisual<P extends ComplexMeshParams>(builder: ComplexMeshVisualBuilder<P>): ComplexVisual<P> {
    return ComplexVisual<StructureMeshParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        createEmptyGeometry: Mesh.createEmpty,
        createRenderObject: createComplexMeshRenderObject,
        updateValues: Mesh.updateValues,
        updateBoundingSphere: Mesh.updateBoundingSphere,
        updateRenderableState: Geometry.updateRenderableState
    })
}

// direct-volume

export const ComplexDirectVolumeParams = {
    ...StructureDirectVolumeParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic', 'spheres', 'gaussians'], UnitKindOptions),
}
export type ComplexDirectVolumeParams = typeof ComplexDirectVolumeParams

export interface ComplexDirectVolumeVisualBuilder<P extends ComplexDirectVolumeParams> extends ComplexVisualBuilder<P, DirectVolume> { }

export function ComplexDirectVolumeVisual<P extends ComplexDirectVolumeParams>(builder: ComplexDirectVolumeVisualBuilder<P>): ComplexVisual<P> {
    return ComplexVisual<StructureDirectVolumeParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        createEmptyGeometry: DirectVolume.createEmpty,
        createRenderObject: createComplexDirectVolumeRenderObject,
        updateValues: DirectVolume.updateValues,
        updateBoundingSphere: DirectVolume.updateBoundingSphere,
        updateRenderableState: DirectVolume.updateRenderableState
    })
}