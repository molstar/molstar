/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RepresentationProps, Representation, Visual, RepresentationContext, VisualContext, RepresentationProvider, RepresentationParamsGetter } from '../representation';
import { VolumeData, VolumeIsoValue } from 'mol-model/volume';
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
import { ThemeProps, Theme, createTheme } from 'mol-theme/theme';

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
    createRenderObject(ctx: VisualContext, geometry: G, locationIt: LocationIterator, theme: Theme, currentProps: P): Promise<VolumeRenderObject>
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

    async function create(ctx: VisualContext, volume: VolumeData, theme: Theme, props: Partial<VolumeProps> = {}) {
        currentProps = Object.assign({}, defaultProps, props)
        if (props.isoValueRelative) {
            currentProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
            // console.log('create props.isoValueRelative', props.isoValueRelative, currentProps.isoValueAbsolute, currentVolume.dataStats)
        }

        geometry = await createGeometry(ctx, volume, currentProps, geometry)
        locationIt = LocationIterator(1, 1, () => NullLocation)
        renderObject = await createRenderObject(ctx, geometry, locationIt, theme, currentProps)
    }

    async function update(ctx: VisualContext, theme: Theme, props: Partial<VolumeProps> = {}) {
        if (!renderObject) return
        const newProps = Object.assign({}, currentProps, props)

        if (props.isoValueRelative) {
            newProps.isoValueAbsolute = VolumeIsoValue.calcAbsolute(currentVolume.dataStats, props.isoValueRelative)
            // console.log('update props.isoValueRelative', props.isoValueRelative, newProps.isoValueAbsolute, currentVolume.dataStats)
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
        async createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<VolumeProps> = {}, volume?: VolumeData) {
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
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> { }

export type VolumeRepresentationProvider<P extends PD.Params> = RepresentationProvider<VolumeData, P>

//

export const VolumeParams = {
    ...Geometry.Params,
    isoValueAbsolute: PD.Range('Iso Value Absolute', '', 0.22, -1, 1, 0.01),
    isoValueRelative: PD.Range('Iso Value Relative', '', 2, -10, 10, 0.1),
}
export const DefaultVolumeProps = PD.getDefaultValues(VolumeParams)
export type VolumeProps = typeof DefaultVolumeProps

export function VolumeRepresentation<P extends VolumeProps>(label: string, getParams: RepresentationParamsGetter<VolumeData>, visualCtor: (volume: VolumeData) => VolumeVisual<P>): VolumeRepresentation<P> {
    let visual: VolumeVisual<P>

    let _volume: VolumeData
    let _props: P
    let _params: PD.Params
    let _theme: Theme
    let busy = false

    function createOrUpdate(ctx: RepresentationContext, props: Partial<P> = {}, themeProps: ThemeProps = {}, volume?: VolumeData) {
        if (volume && volume !== _volume) {
            _params = getParams(ctx, volume)
            _volume = volume
            if (!_props) _props = PD.getDefaultValues(_params) as P
        }
        _props = Object.assign({}, DefaultVolumeProps, _props, props)
        _theme = createTheme(ctx, _props, themeProps, {}, _theme)

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

    return {
        label,
        get renderObjects() {
            return visual && visual.renderObject ? [ visual.renderObject ] : []
        },
        get props () { return _props },
        get params() { return _params },
        createOrUpdate,
        getLoci,
        mark,
        destroy
    }
}