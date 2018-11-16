/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { Representation, Visual, RepresentationContext, VisualContext, RepresentationProvider, RepresentationParamsGetter } from '../representation';
import { VolumeData } from 'mol-model/volume';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { Geometry, updateRenderableState } from 'mol-geo/geometry/geometry';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction, applyMarkerAction } from 'mol-geo/geometry/marker-data';
import { DirectVolumeRenderObject, PointsRenderObject, LinesRenderObject, MeshRenderObject } from 'mol-gl/render-object';
import { Interval } from 'mol-data/int';
import { RenderableValues } from 'mol-gl/renderable/schema';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { NullLocation } from 'mol-model/location';
import { VisualUpdateState } from 'mol-repr/util';
import { ValueCell } from 'mol-util';
import { Theme, createTheme } from 'mol-theme/theme';
import { BehaviorSubject } from 'rxjs';

export interface VolumeVisual<P extends VolumeParams> extends Visual<VolumeData, P> { }

type VolumeRenderObject = MeshRenderObject | LinesRenderObject | PointsRenderObject | DirectVolumeRenderObject

interface VolumeVisualBuilder<P extends VolumeParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, volumeData: VolumeData, props: PD.Values<P>, geometry?: G): Promise<G>
    getLoci(pickingId: PickingId, id: number): Loci
    mark(loci: Loci, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>): void
}

interface VolumeVisualGeometryBuilder<P extends VolumeParams, G extends Geometry> extends VolumeVisualBuilder<P, G> {
    createRenderObject(ctx: VisualContext, geometry: G, locationIt: LocationIterator, theme: Theme, currentProps: PD.Values<P>): Promise<VolumeRenderObject>
    updateValues(values: RenderableValues, newProps: PD.Values<P>): void
}

export function VolumeVisual<P extends VolumeParams>(builder: VolumeVisualGeometryBuilder<P, Geometry>) {
    const { defaultProps, createGeometry, getLoci, mark, setUpdateState } = builder
    const { createRenderObject, updateValues } = builder
    const updateState = VisualUpdateState.create()

    let currentProps: PD.Values<P>
    let renderObject: VolumeRenderObject | undefined
    let currentVolume: VolumeData
    let geometry: Geometry
    let locationIt: LocationIterator

    async function create(ctx: VisualContext, volume: VolumeData, theme: Theme, props: Partial<PD.Values<P>> = {}) {
        currentProps = Object.assign({}, defaultProps, props)
        geometry = await createGeometry(ctx, volume, currentProps, geometry)
        locationIt = LocationIterator(1, 1, () => NullLocation)
        renderObject = await createRenderObject(ctx, geometry, locationIt, theme, currentProps)
    }

    async function update(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}) {
        if (!renderObject) return
        const newProps = Object.assign({}, currentProps, props)

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
        async createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, volume?: VolumeData) {
            if (!volume && !currentVolume) {
                throw new Error('missing volume')
            } else if (volume && (!currentVolume || !renderObject)) {
                currentVolume = volume
                await create(ctx, volume, theme, props)
            } else if (volume && volume !== currentVolume) {
                currentVolume = volume
                await create(ctx, volume, theme, props)
            } else {
                await update(ctx, theme, props)
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
        setVisibility(value: boolean) {
            if (renderObject) renderObject.state.visible = value
        },
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}

export interface VolumeRepresentation<P extends VolumeParams> extends Representation<VolumeData, P> { }

export type VolumeRepresentationProvider<P extends VolumeParams> = RepresentationProvider<VolumeData, P>

//

export const VolumeParams = {
    ...Geometry.Params,
    isoValue: PD.Numeric(0.22, { min: -1, max: 1, step: 0.01 }),
}
export type VolumeParams = typeof VolumeParams

export function VolumeRepresentation<P extends VolumeParams>(label: string, getParams: RepresentationParamsGetter<VolumeData, P>, visualCtor: (volume: VolumeData) => VolumeVisual<P>): VolumeRepresentation<P> {
    const updated = new BehaviorSubject(0)
    let visual: VolumeVisual<P>

    let _volume: VolumeData
    let _props: PD.Values<P>
    let _params: P
    let _theme: Theme
    let busy = false

    function createOrUpdate(ctx: RepresentationContext, props: Partial<PD.Values<P>> = {}, volume?: VolumeData) {
        if (volume && volume !== _volume) {
            _params = getParams(ctx, volume)
            _volume = volume
            if (!_props) _props = PD.getDefaultValues(_params)
        }
        _props = Object.assign({}, _props, props)
        _theme = createTheme(ctx, _props, {}, _theme)

        return Task.create('VolumeRepresentation.create', async runtime => {
            // TODO queue it somehow
            if (busy) return

            if (!visual && !volume) {
                throw new Error('volume data missing')
            } else if (volume && !visual) {
                busy = true
                visual = visualCtor(volume)
                await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, volume)
                busy = false
            } else {
                busy = true
                await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, volume)
                busy = false
            }
            updated.next(updated.getValue() + 1)
        });
    }

    function getLoci(pickingId: PickingId) {
        return visual ? visual.getLoci(pickingId) : EmptyLoci
    }

    function mark(loci: Loci, action: MarkerAction) {
        return visual ? visual.mark(loci, action) : false
    }

    function destroy() {
        if (visual) visual.destroy()
    }

    function setVisibility(value: boolean) {
        if (visual) visual.setVisibility(value)
    }

    return {
        label,
        get renderObjects() {
            return visual && visual.renderObject ? [ visual.renderObject ] : []
        },
        get props () { return _props },
        get params() { return _params },
        get updated() { return updated },
        createOrUpdate,
        getLoci,
        mark,
        setVisibility,
        destroy
    }
}