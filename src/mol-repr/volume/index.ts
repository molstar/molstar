/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RepresentationProps, Representation, Visual, RepresentationContext, VisualContext } from '..';
import { VolumeData, VolumeIsoValue } from 'mol-model/volume';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { paramDefaultValues, RangeParam } from 'mol-util/parameter';
import { Geometry, updateRenderableState } from 'mol-geo/geometry/geometry';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction, applyMarkerAction } from 'mol-geo/geometry/marker-data';
import { DirectVolumeRenderObject, PointsRenderObject, LinesRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { Interval } from 'mol-data/int';
import { RenderableValues } from 'mol-gl/renderable/schema';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { NullLocation } from 'mol-model/location';
import { VisualUpdateState } from 'mol-repr/util';
import { ValueCell } from 'mol-util';

export interface VolumeVisual<P extends RepresentationProps = {}> extends Visual<VolumeData, P> { }

type VolumeRenderObject = MeshRenderObject | LinesRenderObject | PointsRenderObject | DirectVolumeRenderObject

interface VolumeVisualBuilder<P extends VolumeProps, G extends Geometry> {
    defaultProps: P
    createGeometry(ctx: VisualContext, volumeData: VolumeData, props: P, geometry?: G): Promise<G>
    getLoci(pickingId: PickingId, id: number): Loci
    mark(loci: Loci, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: P, currentProps: P): void
}

interface VolumeVisualGeometryBuilder<P extends VolumeProps, G extends Geometry> extends VolumeVisualBuilder<P, G> {
    createRenderObject(ctx: VisualContext, geometry: G, locationIt: LocationIterator, currentProps: P): Promise<VolumeRenderObject>
    updateValues(values: RenderableValues, newProps: P): void
}

export function VolumeVisual<P extends VolumeProps>(builder: VolumeVisualGeometryBuilder<P, Geometry>) {
    const { defaultProps, createGeometry, getLoci, mark, setUpdateState } = builder
    const { createRenderObject, updateValues } = builder
    const updateState = VisualUpdateState.create()

    let currentProps: P
    let renderObject: VolumeRenderObject | undefined
    let currentVolume: VolumeData
    let geometry: Geometry
    let locationIt: LocationIterator

    async function create(ctx: VisualContext, volume: VolumeData, props: Partial<VolumeProps> = {}) {
        currentProps = Object.assign({}, defaultProps, props)
        if (props.isoValueRelative) {
            currentProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
            console.log('create props.isoValueRelative', props.isoValueRelative, currentProps.isoValueAbsolute, currentVolume.dataStats)
        }

        geometry = await createGeometry(ctx, volume, currentProps, geometry)
        locationIt = LocationIterator(1, 1, () => NullLocation)
        renderObject = await createRenderObject(ctx, geometry, locationIt, currentProps)
    }

    async function update(ctx: VisualContext, props: Partial<VolumeProps> = {}) {
        if (!renderObject) return
        const newProps = Object.assign({}, currentProps, props)

        if (props.isoValueRelative) {
            newProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
            console.log('update props.isoValueRelative', props.isoValueRelative, newProps.isoValueAbsolute, currentVolume.dataStats)
        }

        VisualUpdateState.reset(updateState)
        setUpdateState(updateState, newProps, currentProps)

        if (updateState.createGeometry) {
            geometry = await createGeometry(ctx, currentVolume, currentProps, geometry)
            ValueCell.update(renderObject.values.drawCount, Geometry.getDrawCount(geometry))
        }

        updateValues(renderObject.values, newProps)
        updateRenderableState(renderObject.state, newProps)

        currentProps = newProps
    }

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: VisualContext, props: Partial<VolumeProps> = {}, volume?: VolumeData) {
            if (!volume && !currentVolume) {
                throw new Error('missing volume')
            } else if (volume && (!currentVolume || !renderObject)) {
                currentVolume = volume
                await create(ctx, volume, props)
            } else if (volume && volume !== currentVolume) {
                currentVolume = volume
                await create(ctx, volume, props)
            } else {
                await update(ctx, props)
            }

            currentProps = Object.assign({}, defaultProps, props)
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, renderObject.id) : EmptyLoci
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
                changed = mark(loci, apply)
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

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> { }

export const VolumeParams = {
    ...Geometry.Params,
    isoValueAbsolute: RangeParam('Iso Value Absolute', '', 0.22, -1, 1, 0.01),
    isoValueRelative: RangeParam('Iso Value Relative', '', 2, -10, 10, 0.1),
}
export const DefaultVolumeProps = paramDefaultValues(VolumeParams)
export type VolumeProps = typeof DefaultVolumeProps

export function VolumeRepresentation<P extends VolumeProps>(visualCtor: (volumeData: VolumeData) => VolumeVisual<P>): VolumeRepresentation<P> {
    let visual: VolumeVisual<any>
    let _props: P
    let busy = false

    function createOrUpdate(ctx: RepresentationContext, props: Partial<P> = {}, volumeData?: VolumeData) {
        _props = Object.assign({}, DefaultVolumeProps, _props, props)
        return Task.create('VolumeRepresentation.create', async runtime => {
            // TODO queue it somehow
            if (busy) return

            if (!visual && !volumeData) {
                throw new Error('volumeData missing')
            } else if (volumeData && !visual) {
                busy = true
                visual = visualCtor(volumeData)
                await visual.createOrUpdate({ ...ctx, runtime } , props, volumeData)
                busy = false
            } else {
                busy = true
                await visual.createOrUpdate({ ...ctx, runtime }, props, volumeData)
                busy = false
            }
        });
    }

    return {
        label: 'Volume',
        params: VolumeParams,
        get renderObjects() {
            return visual && visual.renderObject ? [ visual.renderObject ] : []
        },
        get props () { return _props },
        createOrUpdate,
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
            return false
        },
        destroy() {
            // TODO
        }
    }
}