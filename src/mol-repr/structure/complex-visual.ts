/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureParams, StructureMeshParams, StructureDirectVolumeParams } from './representation';
import { Visual, VisualContext } from '../visual';
import { Structure } from '../../mol-model/structure';
import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { Theme, createEmptyTheme } from '../../mol-theme/theme';
import { createIdentityTransform } from '../../mol-geo/geometry/transform-data';
import { createRenderObject, RenderObjectKindType, RenderObjectValuesType } from '../../mol-gl/render-object';
import { UnitKind, UnitKindOptions } from './visual/util/common';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci, isEveryLoci, EmptyLoci } from '../../mol-model/loci';
import { Interval } from '../../mol-data/int';
import { VisualUpdateState } from '../util';
import { UnitsParams } from './units-representation';
import { ColorTheme } from '../../mol-theme/color';
import { ValueCell, deepEqual } from '../../mol-util';
import { createSizes } from '../../mol-geo/geometry/size-data';
import { createColors } from '../../mol-geo/geometry/color-data';
import { MarkerAction } from '../../mol-geo/geometry/marker-data';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { SizeTheme } from '../../mol-theme/size';
import { DirectVolume } from '../../mol-geo/geometry/direct-volume/direct-volume';

export interface  ComplexVisual<P extends StructureParams> extends Visual<Structure, P> { }

function createComplexRenderObject<G extends Geometry>(structure: Structure, geometry: G, locationIt: LocationIterator, theme: Theme, props: PD.Values<Geometry.Params<G>>, materialId: number) {
    const { createValues, createRenderableState } = Geometry.getUtils(geometry)
    const transform = createIdentityTransform()
    const values = createValues(geometry, transform, locationIt, theme, props)
    const state = createRenderableState(props)
    return createRenderObject(geometry.kind, values, state, materialId)
}

const ComplexParams = {
    ...StructureParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic', 'spheres'], UnitKindOptions),
}
type ComplexParams = typeof ComplexParams

interface ComplexVisualBuilder<P extends ComplexParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(structure: Structure): LocationIterator
    getLoci(pickingId: PickingId, structure: Structure, id: number): Loci
    eachLocation(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean): boolean,
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme): void
}

interface ComplexVisualGeometryBuilder<P extends UnitsParams, G extends Geometry> extends ComplexVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function ComplexVisual<G extends Geometry, P extends ComplexParams & Geometry.Params<G>>(builder: ComplexVisualGeometryBuilder<P, G>, materialId: number): ComplexVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, eachLocation, setUpdateState } = builder
    const { updateValues, updateBoundingSphere, updateRenderableState } = builder.geometryUtils
    const updateState = VisualUpdateState.create()

    let renderObject: RenderObjectKindType[G['kind']] | undefined

    let newProps: PD.Values<P>
    let newTheme: Theme
    let newStructure: Structure

    let currentProps: PD.Values<P> = Object.assign({}, defaultProps)
    let currentTheme: Theme = createEmptyTheme()
    let currentStructure: Structure

    let geometry: G
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

        if (!ColorTheme.areEqual(theme.color, currentTheme.color)) {
            updateState.updateColor = true
        }

        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) {
            updateState.createGeometry = true
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true
        }
    }

    function update(newGeometry?: G) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newStructure)
            if (newGeometry) {
                renderObject = createComplexRenderObject(newStructure, newGeometry, locationIt, newTheme, newProps, materialId)
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
                    updateBoundingSphere(renderObject.values as RenderObjectValuesType[G['kind']], newGeometry)
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

            updateValues(renderObject.values as RenderObjectValuesType[G['kind']], newProps)
            updateRenderableState(renderObject.state, newProps)
        }

        currentProps = newProps
        currentTheme = newTheme
        currentStructure = newStructure
        if (newGeometry) geometry = newGeometry
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean) {
        if (isEveryLoci(loci) || (Structure.isLoci(loci) && Structure.areEquivalent(loci.structure, currentStructure))) {
            return apply(Interval.ofBounds(0, locationIt.groupCount * locationIt.instanceCount))
        } else {
            return eachLocation(loci, currentStructure, apply)
        }
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
            return Visual.mark(renderObject, loci, action, lociApply)
        },
        setVisibility(visible: boolean) {
            Visual.setVisibility(renderObject, visible)
        },
        setAlphaFactor(alphaFactor: number) {
            Visual.setAlphaFactor(renderObject, alphaFactor)
        },
        setPickable(pickable: boolean) {
            Visual.setPickable(renderObject, pickable)
        },
        setTransform(matrix?: Mat4, instanceMatrices?: Float32Array | null) {
            Visual.setTransform(renderObject, matrix, instanceMatrices)
        },
        setOverpaint(overpaint: Overpaint) {
            return Visual.setOverpaint(renderObject, overpaint, lociApply, true)
        },
        setTransparency(transparency: Transparency) {
            return Visual.setTransparency(renderObject, transparency, lociApply, true)
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

export function ComplexMeshVisual<P extends ComplexMeshParams>(builder: ComplexMeshVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Mesh, StructureMeshParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        geometryUtils: Mesh.Utils
    }, materialId)
}

// direct-volume

export const ComplexDirectVolumeParams = {
    ...StructureDirectVolumeParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic', 'spheres', 'gaussians'], UnitKindOptions),
}
export type ComplexDirectVolumeParams = typeof ComplexDirectVolumeParams

export interface ComplexDirectVolumeVisualBuilder<P extends ComplexDirectVolumeParams> extends ComplexVisualBuilder<P, DirectVolume> { }

export function ComplexDirectVolumeVisual<P extends ComplexDirectVolumeParams>(builder: ComplexDirectVolumeVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<DirectVolume, StructureDirectVolumeParams & UnitsParams>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme)
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true
        },
        geometryUtils: DirectVolume.Utils
    }, materialId)
}