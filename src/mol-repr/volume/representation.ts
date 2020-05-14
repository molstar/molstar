/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Visual, VisualContext } from '../visual';
import { Volume } from '../../mol-model/volume';
import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { Theme } from '../../mol-theme/theme';
import { createIdentityTransform } from '../../mol-geo/geometry/transform-data';
import { createRenderObject, RenderObjectValues, getNextMaterialId, GraphicsRenderObject } from '../../mol-gl/render-object';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci, isEveryLoci, EmptyLoci } from '../../mol-model/loci';
import { Interval } from '../../mol-data/int';
import { VisualUpdateState } from '../util';
import { ColorTheme } from '../../mol-theme/color';
import { ValueCell } from '../../mol-util';
import { createSizes } from '../../mol-geo/geometry/size-data';
import { createColors } from '../../mol-geo/geometry/color-data';
import { MarkerAction } from '../../mol-util/marker-action';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { Representation, RepresentationProvider, RepresentationContext, RepresentationParamsGetter } from '../representation';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { Subject } from 'rxjs';
import { Task } from '../../mol-task';
import { SizeValues } from '../../mol-gl/renderable/schema';
import { Clipping } from '../../mol-theme/clipping';

export interface VolumeVisual<P extends VolumeParams> extends Visual<Volume, P> { }

function createVolumeRenderObject<G extends Geometry>(volume: Volume, geometry: G, locationIt: LocationIterator, theme: Theme, props: PD.Values<Geometry.Params<G>>, materialId: number) {
    const { createValues, createRenderableState } = Geometry.getUtils(geometry);
    const transform = createIdentityTransform();
    const values = createValues(geometry, transform, locationIt, theme, props);
    const state = createRenderableState(props);
    return createRenderObject(geometry.kind, values, state, materialId);
}

interface VolumeVisualBuilder<P extends VolumeParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, volume: Volume, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(volume: Volume): LocationIterator
    getLoci(pickingId: PickingId, volume: Volume, props: PD.Values<P>, id: number): Loci
    eachLocation(loci: Loci, volume: Volume, props: PD.Values<P>, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, volume: Volume, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme): void
}

interface VolumeVisualGeometryBuilder<P extends VolumeParams, G extends Geometry> extends VolumeVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function VolumeVisual<G extends Geometry, P extends VolumeParams & Geometry.Params<G>>(builder: VolumeVisualGeometryBuilder<P, G>, materialId: number): VolumeVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, eachLocation, setUpdateState } = builder;
    const { updateValues, updateBoundingSphere, updateRenderableState } = builder.geometryUtils;
    const updateState = VisualUpdateState.create();

    let renderObject: GraphicsRenderObject<G['kind']> | undefined;

    let newProps: PD.Values<P>;
    let newTheme: Theme;
    let newVolume: Volume;

    let currentProps: PD.Values<P> = Object.assign({}, defaultProps);
    let currentTheme: Theme = Theme.createEmpty();
    let currentVolume: Volume;

    let geometry: G;
    let locationIt: LocationIterator;

    function prepareUpdate(theme: Theme, props: Partial<PD.Values<P>>, volume: Volume) {
        if (!volume && !currentVolume) {
            throw new Error('missing volume');
        }

        newProps = Object.assign({}, currentProps, props);
        newTheme = theme;
        newVolume = volume;

        VisualUpdateState.reset(updateState);

        if (!renderObject) {
            updateState.createNew = true;
        } else if (!currentVolume || !Volume.areEquivalent(newVolume, currentVolume)) {
            updateState.createNew = true;
        }

        if (updateState.createNew) {
            updateState.createGeometry = true;
            return;
        }

        setUpdateState(updateState, volume, newProps, currentProps, newTheme, currentTheme);

        if (!ColorTheme.areEqual(theme.color, currentTheme.color)) updateState.updateColor = true;

        if (updateState.createGeometry) {
            updateState.updateColor = true;
        }
    }

    function update(newGeometry?: G) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newVolume);
            if (newGeometry) {
                renderObject = createVolumeRenderObject(newVolume, newGeometry, locationIt, newTheme, newProps, materialId);
            } else {
                throw new Error('expected geometry to be given');
            }
        } else {
            if (!renderObject) {
                throw new Error('expected renderObject to be available');
            }

            locationIt.reset();

            if (updateState.createGeometry) {
                if (newGeometry) {
                    ValueCell.update(renderObject.values.drawCount, Geometry.getDrawCount(newGeometry));
                    updateBoundingSphere(renderObject.values as RenderObjectValues<G['kind']>, newGeometry);
                } else {
                    throw new Error('expected geometry to be given');
                }
            }

            if (updateState.updateSize) {
                // not all geometries have size data, so check here
                if ('uSize' in renderObject.values) {
                    createSizes(locationIt, newTheme.size, renderObject.values as SizeValues);
                }
            }

            if (updateState.updateColor) {
                createColors(locationIt, newTheme.color, renderObject.values);
            }

            updateValues(renderObject.values, newProps);
            updateRenderableState(renderObject.state, newProps);
        }

        currentProps = newProps;
        currentTheme = newTheme;
        currentVolume = newVolume;
        if (newGeometry) geometry = newGeometry;
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean) {
        if (isEveryLoci(loci)) {
            return apply(Interval.ofBounds(0, locationIt.groupCount * locationIt.instanceCount));
        } else {
            return eachLocation(loci, currentVolume, currentProps, apply);
        }
    }

    return {
        get groupCount() { return locationIt ? locationIt.count : 0; },
        get renderObject () { return renderObject; },
        async createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, volume?: Volume) {
            prepareUpdate(theme, props, volume || currentVolume);
            if (updateState.createGeometry) {
                const newGeometry = createGeometry(ctx, newVolume, newTheme, newProps, geometry);
                return newGeometry instanceof Promise ? newGeometry.then(update) : update(newGeometry);
            } else {
                update();
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentVolume, currentProps, renderObject.id) : EmptyLoci;
        },
        mark(loci: Loci, action: MarkerAction) {
            return Visual.mark(renderObject, loci, action, lociApply);
        },
        setVisibility(visible: boolean) {
            Visual.setVisibility(renderObject, visible);
        },
        setAlphaFactor(alphaFactor: number) {
            Visual.setAlphaFactor(renderObject, alphaFactor);
        },
        setPickable(pickable: boolean) {
            Visual.setPickable(renderObject, pickable);
        },
        setTransform(matrix?: Mat4, instanceMatrices?: Float32Array | null) {
            Visual.setTransform(renderObject, matrix, instanceMatrices);
        },
        setOverpaint(overpaint: Overpaint) {
            return Visual.setOverpaint(renderObject, overpaint, lociApply, true);
        },
        setTransparency(transparency: Transparency) {
            return Visual.setTransparency(renderObject, transparency, lociApply, true);
        },
        setClipping(clipping: Clipping) {
            return Visual.setClipping(renderObject, clipping, lociApply, true);
        },
        destroy() {
            // TODO
            renderObject = undefined;
        }
    };
}

export interface VolumeRepresentation<P extends VolumeParams> extends Representation<Volume, P> { }

export type VolumeRepresentationProvider<P extends VolumeParams, Id extends string = string> = RepresentationProvider<Volume, P, Representation.State, Id>
export function VolumeRepresentationProvider<P extends VolumeParams, Id extends string>(p: VolumeRepresentationProvider<P, Id>): VolumeRepresentationProvider<P, Id> { return p; }

//

export const VolumeParams = {
    ...BaseGeometry.Params,
};
export type VolumeParams = typeof VolumeParams

export function VolumeRepresentation<P extends VolumeParams>(label: string, ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, P>, visualCtor: (materialId: number) => VolumeVisual<P>, getLoci: (volume: Volume, props: PD.Values<P>) => Loci): VolumeRepresentation<P> {
    let version = 0;
    const updated = new Subject<number>();
    const materialId = getNextMaterialId();
    const renderObjects: GraphicsRenderObject[] = [];
    const _state = Representation.createState();
    let visual: VolumeVisual<P> | undefined;

    let _volume: Volume;
    let _params: P;
    let _props: PD.Values<P>;
    let _theme = Theme.createEmpty();

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, volume?: Volume) {
        if (volume && volume !== _volume) {
            _params = getParams(ctx, volume);
            _volume = volume;
            if (!_props) _props = PD.getDefaultValues(_params);
        }
        _props = Object.assign({}, _props, props);

        return Task.create('Creating or updating VolumeRepresentation', async runtime => {
            if (!visual) visual = visualCtor(materialId);
            const promise = visual.createOrUpdate({ webgl: ctx.webgl, runtime }, _theme, _props, volume);
            if (promise) await promise;
            // update list of renderObjects
            renderObjects.length = 0;
            if (visual && visual.renderObject) renderObjects.push(visual.renderObject);
            // increment version
            updated.next(version++);
        });
    }

    function mark(loci: Loci, action: MarkerAction) {
        return visual ? visual.mark(loci, action) : false;
    }

    function setState(state: Partial<Representation.State>) {
        if (state.visible !== undefined && visual) visual.setVisibility(state.visible);
        if (state.alphaFactor !== undefined && visual) visual.setAlphaFactor(state.alphaFactor);
        if (state.pickable !== undefined && visual) visual.setPickable(state.pickable);
        if (state.overpaint !== undefined && visual) visual.setOverpaint(state.overpaint);
        if (state.transparency !== undefined && visual) visual.setTransparency(state.transparency);
        if (state.transform !== undefined && visual) visual.setTransform(state.transform);

        Representation.updateState(_state, state);
    }

    function setTheme(theme: Theme) {
        _theme = theme;
    }

    function destroy() {
        if (visual) visual.destroy();
    }

    return {
        label,
        get groupCount() {
            return visual ? visual.groupCount : 0;
        },
        get props () { return _props; },
        get params() { return _params; },
        get state() { return _state; },
        get theme() { return _theme; },
        renderObjects,
        updated,
        createOrUpdate,
        setState,
        setTheme,
        getLoci: (pickingId?: PickingId): Loci => {
            if (pickingId === undefined) return getLoci(_volume, _props);
            return visual ? visual.getLoci(pickingId) : EmptyLoci;
        },
        mark,
        destroy
    };
}