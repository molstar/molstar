/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { Visual, VisualContext } from '..';
import { MeshRenderObject, LinesRenderObject, PointsRenderObject, DirectVolumeRenderObject } from 'mol-gl/render-object';
import { createComplexMeshRenderObject, UnitKind, UnitKindOptions } from './visual/util/common';
import { StructureProps, StructureMeshParams, StructureParams } from './index';
import { deepEqual, ValueCell } from 'mol-util';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { Interval } from 'mol-data/int';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { RenderableValues } from 'mol-gl/renderable/schema';
import { createSizes } from 'mol-geo/geometry/size-data';
import { Geometry, updateRenderableState, Theme } from 'mol-geo/geometry/geometry';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { createColors } from 'mol-geo/geometry/color-data';
import { MarkerAction, applyMarkerAction } from 'mol-geo/geometry/marker-data';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { VisualUpdateState, colorChanged, sizeChanged } from 'mol-repr/util';

export interface  ComplexVisual<P extends StructureProps> extends Visual<Structure, P> { }

const ComplexParams = {
    ...StructureParams,
    unitKinds: PD.MultiSelect<UnitKind>('Unit Kind', '', ['atomic', 'spheres'], UnitKindOptions),
}
const DefaultComplexProps = PD.getDefaultValues(ComplexParams)
type ComplexProps = typeof DefaultComplexProps

type ComplexRenderObject = MeshRenderObject | LinesRenderObject | PointsRenderObject | DirectVolumeRenderObject

interface ComplexVisualBuilder<P extends ComplexProps, G extends Geometry> {
    defaultProps: P
    createGeometry(ctx: VisualContext, structure: Structure, theme: Theme, props: P, geometry?: G): Promise<G>
    createLocationIterator(structure: Structure): LocationIterator
    getLoci(pickingId: PickingId, structure: Structure, id: number): Loci
    mark(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean): boolean,
    setUpdateState(state: VisualUpdateState, newProps: P, currentProps: P): void
}

interface ComplexVisualGeometryBuilder<P extends ComplexProps, G extends Geometry> extends ComplexVisualBuilder<P, G> {
    createEmptyGeometry(geometry?: G): G
    createRenderObject(ctx: VisualContext, structure: Structure, geometry: Geometry, locationIt: LocationIterator, theme: Theme, currentProps: P): Promise<ComplexRenderObject>
    updateValues(values: RenderableValues, newProps: P): void
}

export function ComplexVisual<P extends ComplexMeshProps>(builder: ComplexVisualGeometryBuilder<P, Geometry>): ComplexVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const { createRenderObject, updateValues } = builder
    const updateState = VisualUpdateState.create()

    let renderObject: ComplexRenderObject | undefined
    let currentProps: P
    let geometry: Geometry
    let currentStructure: Structure
    let locationIt: LocationIterator
    let conformationHash: number

    async function create(ctx: VisualContext, structure: Structure, theme: Theme, props: Partial<P> = {}) {
        currentProps = Object.assign({}, defaultProps, props)
        currentStructure = structure

        conformationHash = Structure.conformationHash(currentStructure)
        geometry = await createGeometry(ctx, currentStructure, theme, currentProps, geometry)

        locationIt = createLocationIterator(structure)
        renderObject = await createRenderObject(ctx, structure, geometry, locationIt, theme, currentProps)
    }

    async function update(ctx: VisualContext, theme: Theme, props: Partial<P>) {
        const newProps = Object.assign({}, currentProps, props, { structure: currentStructure })

        if (!renderObject) return false

        locationIt.reset()
        VisualUpdateState.reset(updateState)
        setUpdateState(updateState, newProps, currentProps)

        const newConformationHash = Structure.conformationHash(currentStructure)
        if (newConformationHash !== conformationHash) {
            conformationHash = newConformationHash
            updateState.createGeometry = true
        }

        if (colorChanged(currentProps, newProps)) updateState.updateColor = true
        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) updateState.createGeometry = true

        //

        if (updateState.createGeometry) {
            geometry = await createGeometry(ctx, currentStructure, theme, newProps, geometry)
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
        return true
    }

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<P> = {}, structure?: Structure) {
            if (!structure && !currentStructure) {
                throw new Error('missing structure')
            } else if (structure && (!currentStructure || !renderObject)) {
                await create(ctx, structure, theme, props)
            } else if (structure && structure.hashCode !== currentStructure.hashCode) {
                await create(ctx, structure, theme, props)
            } else {
                if (structure && Structure.conformationHash(structure) !== Structure.conformationHash(currentStructure)) {
                    currentStructure = structure
                }
                await update(ctx, theme, props)
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
            if (isEveryLoci(loci)) {
                changed = apply(Interval.ofBounds(0, groupCount * instanceCount))
            } else {
                changed = mark(loci, currentStructure, apply)
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

export const ComplexMeshParams = {
    ...StructureMeshParams,
    unitKinds: PD.MultiSelect<UnitKind>('Unit Kind', '', [ 'atomic', 'spheres' ], UnitKindOptions),
}
export const DefaultComplexMeshProps = PD.getDefaultValues(ComplexMeshParams)
export type ComplexMeshProps = typeof DefaultComplexMeshProps

export interface ComplexMeshVisualBuilder<P extends ComplexMeshProps> extends ComplexVisualBuilder<P, Mesh> { }

export function ComplexMeshVisual<P extends ComplexMeshProps>(builder: ComplexMeshVisualBuilder<P>): ComplexVisual<P> {
    return ComplexVisual({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: P, currentProps: P) => {
            builder.setUpdateState(state, newProps, currentProps)
            if (sizeChanged(currentProps, newProps)) state.createGeometry = true
        },
        createEmptyGeometry: Mesh.createEmpty,
        createRenderObject: createComplexMeshRenderObject,
        updateValues: Mesh.updateValues
    })
}