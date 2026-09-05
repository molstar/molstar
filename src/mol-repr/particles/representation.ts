/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParticleList } from '../../mol-model/particles/particle-list';
import { Theme } from '../../mol-theme/theme';
import { getNextMaterialId, GraphicsRenderObject } from '../../mol-gl/render-object';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci as ModelLoci, EmptyLoci } from '../../mol-model/loci';
import { getQualityProps, LocationCallback } from '../util';
import { MarkerAction } from '../../mol-util/marker-action';
import { EPSILON, Mat4 } from '../../mol-math/linear-algebra';
import { Representation, RepresentationProvider, RepresentationContext, RepresentationParamsGetter } from '../representation';
import { Subject } from 'rxjs';
import { Task } from '../../mol-task';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { ParticleVisual, ParticleParams } from './visual';

export interface ParticleRepresentation<P extends PD.Params = PD.Params> extends Representation<ParticleList, P> { }

export type ParticleRepresentationProvider<P extends PD.Params, Id extends string = string> = RepresentationProvider<ParticleList, P, Representation.State, Id>
export function ParticleRepresentationProvider<P extends PD.Params, Id extends string>(p: ParticleRepresentationProvider<P, Id>): ParticleRepresentationProvider<P, Id> { return p; }

//

export function ParticleRepresentation<P extends ParticleParams>(label: string, ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, P>, visualCtor: (materialId: number, particles: ParticleList, props: PD.Values<P>, webgl?: WebGLContext) => ParticleVisual<P>, getLoci: (particles: ParticleList, props: PD.Values<P>) => ModelLoci): ParticleRepresentation<P> {
    let version = 0;
    const { webgl } = ctx;
    const updated = new Subject<number>();
    const geometryState = new Representation.GeometryState();
    const materialId = getNextMaterialId();
    const renderObjects: GraphicsRenderObject[] = [];
    const _state = Representation.createState();
    let _theme = Theme.createEmpty();

    let visual: ParticleVisual<P> | undefined;
    let _particles: ParticleList;
    let _params: P;
    let _props: PD.Values<P>;

    function setVisualState(state: Partial<Representation.State>) {
        if (!visual) return;
        if (state.visible !== undefined) visual.setVisibility(state.visible);
        if (state.alphaFactor !== undefined) visual.setAlphaFactor(state.alphaFactor);
        if (state.pickable !== undefined) visual.setPickable(state.pickable);
        if (state.colorOnly !== undefined) visual.setColorOnly(state.colorOnly);
        if (state.overpaint !== undefined) visual.setOverpaint(state.overpaint);
        if (state.transparency !== undefined) visual.setTransparency(state.transparency);
        if (state.emissive !== undefined) visual.setEmissive(state.emissive);
        if (state.substance !== undefined) visual.setSubstance(state.substance);
        if (state.clipping !== undefined) visual.setClipping(state.clipping);
        if (state.wiggle !== undefined) visual.setWiggle(state.wiggle);
        if (state.transform !== undefined) visual.setTransform(state.transform);
        if (state.themeStrength !== undefined) visual.setThemeStrength(state.themeStrength);
    }

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, particles?: ParticleList) {
        if (particles && particles !== _particles) {
            _params = getParams(ctx, particles);
            _particles = particles;
            if (!_props) _props = PD.getDefaultValues(_params);
        }
        const qualityProps = getQualityProps(Object.assign({}, _props, props), _particles);
        Object.assign(_props, props, qualityProps);

        return Task.create(`Creating or updating '${label}' representation`, async runtime => {
            if (!visual) {
                visual = visualCtor(materialId, _particles, _props, webgl);
            } else if (visual.mustRecreate?.({ particles: _particles }, _props, webgl)) {
                visual.destroy();
                visual = visualCtor(materialId, _particles, _props, webgl);
            }
            await visual.createOrUpdate({ webgl, runtime }, _theme, _props, { particles: _particles });
            renderObjects.length = 0;
            if (visual.renderObject) {
                renderObjects.push(visual.renderObject);
                geometryState.add(visual.renderObject.id, visual.geometryVersion);
            }
            geometryState.snapshot();
            updated.next(version++);
        });
    }

    function setState(state: Partial<Representation.State>) {
        const newState: Partial<Representation.State> = {};
        if (state.visible !== undefined) newState.visible = state.visible;
        if (state.alphaFactor !== undefined) newState.alphaFactor = state.alphaFactor;
        if (state.pickable !== undefined) newState.pickable = state.pickable;
        if (state.colorOnly !== undefined) newState.colorOnly = state.colorOnly;
        if (state.overpaint !== undefined) newState.overpaint = state.overpaint;
        if (state.transparency !== undefined) newState.transparency = state.transparency;
        if (state.emissive !== undefined) newState.emissive = state.emissive;
        if (state.substance !== undefined) newState.substance = state.substance;
        if (state.clipping !== undefined) newState.clipping = state.clipping;
        if (state.wiggle !== undefined) newState.wiggle = state.wiggle;
        if (state.themeStrength !== undefined) newState.themeStrength = state.themeStrength;
        if (state.transform !== undefined && !Mat4.areEqual(state.transform, _state.transform, EPSILON)) {
            newState.transform = state.transform;
        }
        if (state.syncManually !== undefined) newState.syncManually = state.syncManually;
        if (state.markerActions !== undefined) newState.markerActions = state.markerActions;

        setVisualState(newState);
        Representation.updateState(_state, state);
    }

    function setTheme(theme: Theme) {
        _theme = theme;
    }

    function destroy() {
        visual?.destroy();
        visual = undefined;
        renderObjects.length = 0;
    }

    return {
        label,
        updated,
        get groupCount() { return visual?.groupCount ?? 0; },
        renderObjects,
        get geometryVersion() { return geometryState.version; },
        get props() { return _props; },
        get params() { return _params; },
        get state() { return _state; },
        get theme() { return _theme; },
        createOrUpdate,
        setState,
        setTheme,
        getLoci: (pickingId: PickingId) => {
            return visual?.getLoci(pickingId) ?? EmptyLoci;
        },
        getAllLoci: () => {
            if (!_particles) return [];
            return [getLoci(_particles, _props)];
        },
        eachLocation: (cb: LocationCallback) => {
            visual?.eachLocation(cb);
        },
        mark: (loci: ModelLoci, action: MarkerAction) => {
            return visual?.mark(loci, action) ?? false;
        },
        destroy,
    };
}
