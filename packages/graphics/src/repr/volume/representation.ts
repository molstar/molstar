/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Volume } from '@molstar/model/model/volume';
import { Theme } from '@molstar/graphics/theme/theme';
import { getNextMaterialId, type GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { type Loci, EmptyLoci, isEmptyLoci } from '@molstar/model/model/loci';
import { getQualityProps, type LocationCallback } from '../util.js';
import type { MarkerAction } from '@molstar/core/util/marker-action';
import { EPSILON, Mat4 } from '@molstar/core/math/linear-algebra';
import {
  Representation,
  type RepresentationProvider,
  type RepresentationContext,
  type RepresentationParamsGetter,
} from '../representation.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { Subject } from 'rxjs';
import { type RuntimeContext, Task } from '@molstar/core/task';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { VolumeKey, VolumeVisual } from './visual.js';
import { getQualitySignatures, type SignatureContext } from '../signatures.js';
import { createVisual, getFactorySignatures, type VisualFactory } from '../visual-factory.js';

export interface VolumeRepresentation<P extends VolumeParams> extends Representation<Volume, P> {}

export type VolumeRepresentationProvider<P extends VolumeParams, Id extends string = string> = RepresentationProvider<
  Volume,
  P,
  Representation.State,
  Id
>;
export function VolumeRepresentationProvider<P extends VolumeParams, Id extends string>(
  p: VolumeRepresentationProvider<P, Id>,
): VolumeRepresentationProvider<P, Id> {
  return p;
}

//

export const VolumeParams = {
  ...BaseGeometry.Params,
};
export type VolumeParams = typeof VolumeParams;

export function VolumeRepresentation<P extends VolumeParams>(
  label: string,
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Volume, P>,
  visualCtor: VisualFactory<VolumeKey, P>,
  getLoci: (volume: Volume, props: PD.Values<P>) => Loci,
  getKeys: (props: PD.Values<P>, volume: Volume) => ArrayLike<number> = () => [-1],
): VolumeRepresentation<P> {
  let version = 0;
  let { graphics } = ctx;
  const updated = new Subject<number>();
  const geometryState = new Representation.GeometryState();
  const materialId = getNextMaterialId();
  const renderObjects: GraphicsRenderObject[] = [];
  const _state = Representation.createState();
  const visuals = new Map<number, VolumeVisual<P>>();

  let _volume: Volume;
  let _keys: ArrayLike<number>;
  let _params: P;
  let _props: PD.Values<P>;
  let _theme = Theme.createEmpty();

  async function visual(runtime: RuntimeContext, key: number) {
    let visual = visuals.get(key);
    if (!visual) {
      visual = createVisual(visualCtor, materialId, { volume: _volume, key }, _props, graphics);
      visuals.set(key, visual);
    } else if (visual.mustRecreate?.({ volume: _volume, key }, _props, graphics)) {
      visual.destroy();
      visual = createVisual(visualCtor, materialId, { volume: _volume, key }, _props, graphics);
      visuals.set(key, visual);
    }
    const renderObject = visual.renderObject;
    await visual.createOrUpdate({ graphics, runtime }, _theme, _props, { volume: _volume, key });
    // a new render object (new visual, `mustRecreate`, data that is not equivalent) starts from the defaults: give it
    // the current state, transform included (perf-audit-2 07/B1)
    if (visual.renderObject && visual.renderObject !== renderObject) setVisualState(visual, _state);
  }

  /** updates one after another, see `Representation.createUpdateQueue` */
  const serialize = Representation.createUpdateQueue();
  function createOrUpdate(props: Partial<PD.Values<P>> = {}, volume?: Volume) {
    if (volume && volume !== _volume) {
      _params = getParams(ctx, volume);
      _volume = volume;
      if (!_props) _props = PD.getDefaultValues(_params);
    }
    const qualityProps = getQualityProps(Object.assign({}, _props, props), _volume);
    Object.assign(_props, props, qualityProps);
    _keys = getKeys(_props, _volume);

    const task = Task.create('Creating or updating VolumeRepresentation', async (runtime) => {
      const toDelete = new Set(visuals.keys());
      for (let i = 0, il = _keys.length; i < il; ++i) {
        const key = _keys[i];
        toDelete.delete(key);
        const promise = visual(runtime, key);
        if (promise) await promise;
      }
      toDelete.forEach((key) => {
        visuals.get(key)?.destroy();
        visuals.delete(key);
      });
      // update list of renderObjects
      renderObjects.length = 0;
      visuals.forEach((visual) => {
        if (visual.renderObject) {
          renderObjects.push(visual.renderObject);
          geometryState.add(visual.renderObject.id, visual.geometryVersion);
        }
      });
      geometryState.snapshot();
      // increment version
      updated.next(version++);
    });
    return serialize(task);
  }

  function mark(loci: Loci, action: MarkerAction) {
    let changed = false;
    visuals.forEach((visual) => {
      changed = visual.mark(loci, action) || changed;
    });
    return changed;
  }

  function setVisualState(visual: VolumeVisual<P>, state: Partial<Representation.State>) {
    if (state.visible !== undefined && visual) visual.setVisibility(state.visible);
    if (state.pickable !== undefined && visual) visual.setPickable(state.pickable);
    if (state.overpaint !== undefined && visual) visual.setOverpaint(state.overpaint);
    if (state.transparency !== undefined && visual) visual.setTransparency(state.transparency);
    if (state.emissive !== undefined && visual) visual.setEmissive(state.emissive);
    if (state.substance !== undefined && visual) visual.setSubstance(state.substance);
    if (state.clipping !== undefined && visual) visual.setClipping(state.clipping);
    if (state.wiggle !== undefined && visual) visual.setWiggle(state.wiggle);
    if (state.transform !== undefined && visual) visual.setTransform(state.transform);
    if (state.themeStrength !== undefined && visual) visual.setThemeStrength(state.themeStrength);
  }

  function setState(state: Partial<Representation.State>) {
    const {
      visible,
      pickable,
      overpaint,
      transparency,
      emissive,
      substance,
      clipping,
      wiggle,
      transform,
      themeStrength,
      syncManually,
      markerActions,
    } = state;
    const newState: Partial<Representation.State> = {};

    if (visible !== undefined) newState.visible = visible;
    if (pickable !== undefined) newState.pickable = pickable;
    if (overpaint !== undefined) newState.overpaint = overpaint;
    if (transparency !== undefined) newState.transparency = transparency;
    if (emissive !== undefined) newState.emissive = emissive;
    if (substance !== undefined) newState.substance = substance;
    if (clipping !== undefined) newState.clipping = clipping;
    if (wiggle !== undefined) newState.wiggle = wiggle;
    if (themeStrength !== undefined) newState.themeStrength = themeStrength;
    if (transform !== undefined && !Mat4.areEqual(transform, _state.transform, EPSILON)) {
      newState.transform = transform;
    }
    if (syncManually !== undefined) newState.syncManually = syncManually;
    if (markerActions !== undefined) newState.markerActions = markerActions;

    visuals.forEach((visual) => setVisualState(visual, newState));

    Representation.updateState(_state, state);
  }

  function setTheme(theme: Theme) {
    _theme = theme;
  }

  function destroy() {
    visuals.forEach((visual) => visual.destroy());
    visuals.clear();
  }

  return {
    label,
    get groupCount() {
      let groupCount = 0;
      visuals.forEach((visual) => {
        if (visual.renderObject) groupCount += visual.groupCount;
      });
      return groupCount;
    },
    get props() {
      return _props;
    },
    get params() {
      return _params;
    },
    get state() {
      return _state;
    },
    get theme() {
      return _theme;
    },
    get geometryVersion() {
      return geometryState.version;
    },
    renderObjects,
    updated,
    createOrUpdate,
    setState,
    setTheme,
    getLoci: (pickingId: PickingId): Loci => {
      let loci: Loci = EmptyLoci;
      for (const visual of visuals.values()) {
        const _loci = visual.getLoci(pickingId);
        if (!isEmptyLoci(_loci)) {
          loci = _loci;
          break;
        }
      }
      return loci;
    },
    getAllLoci: (): Loci[] => {
      return [getLoci(_volume, _props)];
    },
    eachLocation: (cb: LocationCallback) => {
      visuals.forEach((visual) => {
        visual.eachLocation(cb);
      });
    },
    mark,
    setGraphics: (value: GraphicsContext) => {
      graphics = value;
    },
    getSignatures: (props: PD.Values<P>, ctx: SignatureContext) =>
      getQualitySignatures(props, (props) => getFactorySignatures(visualCtor, props, ctx)),
    destroy,
  };
}
