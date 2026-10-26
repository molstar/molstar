/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RuntimeContext } from '@molstar/core/task';
import { createForeignResourcesCheck, type GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { Loci, isEmptyLoci, isEveryLoci, EveryLoci } from '@molstar/model/model/loci';
import {
  MarkerAction,
  applyMarkerAction,
  getMarkerInfo,
  setMarkerValue,
  getPartialMarkerAverage,
  MarkerActions,
  type MarkerInfo,
} from '@molstar/core/util/marker-action';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { updateTransformData } from '@molstar/graphics/geo/geometry/transform-data';
import { calculateTransformBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import { ValueCell } from '@molstar/core/util';
import { TextureChanges } from '@molstar/graphics/gl/renderable/texture-changes';
import type { Overpaint } from '@molstar/graphics/theme/overpaint';
import { createOverpaint, clearOverpaint, applyOverpaintColor } from '@molstar/graphics/geo/geometry/overpaint-data';
import { Interval } from '@molstar/core/data/int/interval';
import type { Transparency } from '@molstar/graphics/theme/transparency';
import {
  createTransparency,
  clearTransparency,
  applyTransparencyValue,
  getTransparencyAverage,
  getTransparencyMin,
} from '@molstar/graphics/geo/geometry/transparency-data';
import type { Clipping } from '@molstar/graphics/theme/clipping';
import { createClipping, applyClippingGroups, clearClipping } from '@molstar/graphics/geo/geometry/clipping-data';
import { getMarkersAverage } from '@molstar/graphics/geo/geometry/marker-data';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { Geometry, GeometryUtils } from '@molstar/graphics/geo/geometry/geometry';
import { GeometryRequest } from '@molstar/graphics/geo/geometry/request';
import { isPromiseLike } from '@molstar/core/util/type-helpers';
import { getColorSmoothingProps, hasColorSmoothingProp } from '@molstar/graphics/geo/geometry/base';
import {
  applyMeshEmissiveSmoothing,
  applyMeshOverpaintSmoothing,
  applyMeshSubstanceSmoothing,
  applyMeshTransparencySmoothing,
} from '@molstar/graphics/geo/geometry/mesh/color-smoothing';
import {
  applyTextureMeshEmissiveSmoothing,
  applyTextureMeshOverpaintSmoothing,
  applyTextureMeshSubstanceSmoothing,
  applyTextureMeshTransparencySmoothing,
} from '@molstar/graphics/geo/geometry/texture-mesh/color-smoothing';
import type { Substance } from '@molstar/graphics/theme/substance';
import { applySubstanceMaterial, clearSubstance, createSubstance } from '@molstar/graphics/geo/geometry/substance-data';
import type { LocationCallback } from './util.js';
import type { SignatureContext } from './signatures.js';
import type { ObjectSignature } from '@molstar/graphics/gl/precompile';
import type { Emissive } from '@molstar/graphics/theme/emissive';
import {
  applyEmissiveValue,
  clearEmissive,
  createEmissive,
  getEmissiveAverage,
} from '@molstar/graphics/geo/geometry/emissive-data';
import type { Wiggle } from '@molstar/graphics/theme/wiggle';
import {
  applyWiggleValue,
  clearWiggle,
  createWiggle,
  getWiggleAverage,
} from '@molstar/graphics/geo/geometry/wiggle-data';

export interface VisualContext {
  readonly runtime: RuntimeContext;
  readonly graphics?: GraphicsContext;
}

export { Visual };
interface Visual<D, P extends PD.Params> {
  /** Number of addressable groups in all instances of the visual */
  readonly groupCount: number;
  readonly renderObject: GraphicsRenderObject | undefined;
  readonly geometryVersion: number;
  createOrUpdate: (ctx: VisualContext, theme: Theme, props: PD.Values<P>, data?: D) => Promise<void> | void;
  getLoci: (pickingId: PickingId) => Loci;
  eachLocation: (cb: LocationCallback) => void;
  mark: (loci: Loci, action: MarkerAction) => boolean;
  setVisibility: (visible: boolean) => void;
  setPickable: (pickable: boolean) => void;
  setColorOnly: (colorOnly: boolean) => void;
  setTransform: (matrix?: Mat4, instanceMatrices?: Float32Array | null) => void;
  /**
   * The layer setters smooth the layer where the visual's props ask for it, unless `smooth` is false (members of a
   * merged render object: the merged object has one layer type for all members, and smoothing changes it)
   */
  setOverpaint: (overpaint: Overpaint, graphics?: GraphicsContext, smooth?: boolean) => void;
  setTransparency: (transparency: Transparency, graphics?: GraphicsContext, smooth?: boolean) => void;
  setEmissive: (emissive: Emissive, graphics?: GraphicsContext, smooth?: boolean) => void;
  setSubstance: (substance: Substance, graphics?: GraphicsContext, smooth?: boolean) => void;
  setClipping: (clipping: Clipping) => void;
  setWiggle: (wiggle: Wiggle, graphics?: GraphicsContext) => void;
  setThemeStrength: (strength: {
    overpaint: number;
    transparency: number;
    emissive: number;
    substance: number;
    wiggle: number;
  }) => void;
  destroy: () => void;
  mustRecreate?: (data: D, props: PD.Values<P>, graphics?: GraphicsContext) => boolean;
  /**
   * The signatures of the objects the visual may create with `props` in `ctx`, without data (see
   * `getGeometrySignatures`)
   */
  getSignatures?(props: PD.Values<P>, ctx: SignatureContext): ObjectSignature[];
}
namespace Visual {
  /**
   * `createForeignResourcesCheck` that also knows the render object the visual created with a context: its resources
   * are that context's, so the first update after a build checks nothing (perf-audit-2 07/13/P4).
   */
  export function foreignResourcesCheck() {
    const check = createForeignResourcesCheck();
    let createdObject: GraphicsRenderObject | undefined;
    let createdGraphics: GraphicsContext | undefined;
    return {
      created(o: GraphicsRenderObject, graphics: GraphicsContext | undefined) {
        createdObject = graphics ? o : undefined;
        createdGraphics = graphics;
      },
      hasForeign(o: GraphicsRenderObject, graphics: GraphicsContext) {
        if (o === createdObject && graphics === createdGraphics) return false;
        return check(o, graphics);
      },
    };
  }

  /**
   * Runs `update` with `geometry` (created, possibly asynchronously, or kept) once it is prepared for `theme`
   * (`GeometryUtils.prepare`, e.g. GPU data read back for the position iterator); synchronous where creating and
   * preparing are
   */
  export function withPreparedGeometry<G extends Geometry>(
    geometry: G | Promise<G>,
    utils: GeometryUtils<G>,
    theme: Theme,
    graphics: GraphicsContext | undefined,
    update: (geometry: G) => void,
  ): Promise<void> | void {
    const prepared = (g: G): Promise<void> | void => {
      const preparing = utils.prepare?.(g, GeometryRequest.fromTheme(theme), graphics);
      if (preparing) return preparing.then(() => update(g));
      update(g);
    };
    return isPromiseLike(geometry) ? geometry.then(prepared) : prepared(geometry);
  }

  export type LociApply = (loci: Loci, apply: (interval: Interval) => boolean, isMarking: boolean) => boolean;

  export function setVisibility(renderObject: GraphicsRenderObject | undefined, visible: boolean) {
    if (renderObject) renderObject.state.visible = visible;
  }

  export function setPickable(renderObject: GraphicsRenderObject | undefined, pickable: boolean) {
    if (renderObject) renderObject.state.pickable = pickable;
  }

  export function setColorOnly(renderObject: GraphicsRenderObject | undefined, colorOnly: boolean) {
    if (renderObject) renderObject.state.colorOnly = colorOnly;
  }

  export type PreviousMark = { loci: Loci; action: MarkerAction; status: MarkerInfo['status'] };

  export function mark(
    renderObject: GraphicsRenderObject | undefined,
    loci: Loci,
    action: MarkerAction,
    lociApply: LociApply,
    previous?: PreviousMark,
  ) {
    if (!renderObject || isEmptyLoci(loci)) return false;

    const {
      tMarker,
      uMarker,
      markerAverage,
      markerStatus,
      uGroupCount,
      instanceCount,
      uInstanceGranularity: instanceGranularity,
    } = renderObject.values;
    const count = instanceGranularity.ref.value
      ? instanceCount.ref.value
      : uGroupCount.ref.value * instanceCount.ref.value;
    const { array } = tMarker.ref.value;
    const currentStatus = markerStatus.ref.value as MarkerInfo['status'];
    // all clear: an action that only removes changes nothing; skips the loci's lookup and a marking that changed nothing
    if (
      currentStatus === 0 &&
      (action === MarkerAction.RemoveHighlight || action === MarkerAction.Deselect || action === MarkerAction.Clear)
    ) {
      return false;
    }

    let changed = false;
    let average = -1;
    let status: MarkerInfo['status'] = -1;
    /** the element ranges changed, uploaded instead of the layer (`TextureChanges`); undefined: all of it */
    let ranges: number[] | undefined;
    if (isEveryLoci(loci)) {
      const info = getMarkerInfo(action, currentStatus);
      if (info.status !== -1) {
        changed = currentStatus !== info.status;
        if (changed) setMarkerValue(array, info.status, count);
      } else {
        changed = applyMarkerAction(array, Interval.ofLength(count), action);
      }
      average = info.average;
      status = info.status;
    } else {
      // one pass that applies and sums the covered size, intervals assumed non-overlapping (perf-audit 07 P5)
      let size = 0;
      ranges = [];
      lociApply(
        loci,
        (interval) => {
          size += Interval.size(interval);
          if (applyMarkerAction(array, interval, action)) changed = true;
          ranges!.push(Interval.start(interval), Interval.end(interval));
          return true;
        },
        true,
      );
      if (size === 0) return false;
      if (size === count) {
        // all covered: the outcome is that of `EveryLoci`
        loci = EveryLoci;
        const info = getMarkerInfo(action, currentStatus);
        if (info.status !== -1) {
          changed = currentStatus !== info.status;
          if (changed) setMarkerValue(array, info.status, count);
        }
        average = info.average;
        status = info.status;
      } else if (changed) {
        average = getPartialMarkerAverage(action, currentStatus);
        if (
          previous &&
          previous.status !== -1 &&
          average === -1 &&
          MarkerActions.isReverse(previous.action, action) &&
          Loci.areEqual(loci, previous.loci)
        ) {
          status = previous.status;
          average = status === 0 ? 0 : 0.5;
        }
      }
    }
    if (changed) {
      if (average === -1) average = getMarkersAverage(array, count);
      // nothing marked: the status says so, whichever way the average was found (a partial average can be 0)
      if (average === 0) status = 0;
      if (previous) {
        previous.action = action;
        previous.loci = loci;
        previous.status = currentStatus;
      }
      ValueCell.updateIfChanged(uMarker, status);
      // with a status the uniform stands for the layer (`uMarker`): the image changed all the same, uploaded with the
      // next update of the layer (only its changed ranges are)
      if (status === -1) TextureChanges.update(tMarker, ranges);
      else TextureChanges.defer(tMarker, ranges);
      ValueCell.updateIfChanged(markerAverage, average);
      ValueCell.updateIfChanged(markerStatus, status);
    }
    return changed;
  }

  type KnownClear = { readonly count: number; readonly type: string; readonly version: number };
  type LayerName = 'overpaint' | 'transparency' | 'emissive' | 'substance' | 'clipping' | 'wiggle';

  /**
   * Layers of a render object that the last layer change left all clear (it reached none of their groups), with the
   * texture cell version and size they were cleared at: the next change that reaches none of them either needs no
   * clearing and no texture upload, so a layer change costs the visuals it reaches (perf-audit 07 P1)
   */
  const KnownClearLayers = new WeakMap<GraphicsRenderObject, { [K in LayerName]?: KnownClear }>();

  function isKnownClear(
    renderObject: GraphicsRenderObject,
    name: LayerName,
    texture: ValueCell<unknown>,
    count: number,
    type: string,
  ) {
    const known = KnownClearLayers.get(renderObject)?.[name];
    return !!known && known.count === count && known.type === type && known.version === texture.ref.version;
  }

  function setKnownClear(
    renderObject: GraphicsRenderObject,
    name: LayerName,
    texture: ValueCell<unknown>,
    count: number,
    type: string,
    isClear: boolean,
  ) {
    let layers = KnownClearLayers.get(renderObject);
    if (!isClear) {
      if (layers) layers[name] = undefined;
      return;
    }
    if (!layers) {
      layers = {};
      KnownClearLayers.set(renderObject, layers);
    }
    layers[name] = { count, type, version: texture.ref.version };
  }

  /** whether any of `layers` reaches a group, applying nothing */
  function reachesAny(layers: ReadonlyArray<{ readonly loci: Loci }>, lociApply: LociApply) {
    let reached = false;
    const probe = () => {
      reached = true;
      return true;
    };
    for (let i = 0, il = layers.length; i < il && !reached; ++i) lociApply(layers[i].loci, probe, false);
    return reached;
  }

  type SurfaceMeta = {
    resolution?: number;
    overpaintTexture?: Texture;
    transparencyTexture?: Texture;
    emissiveTexture?: Texture;
    substanceTexture?: Texture;
  };

  type SmoothingContext = {
    geometry: Geometry;
    props: PD.Values<any>;
    graphics?: GraphicsContext;
  };

  export function setOverpaint(
    renderObject: GraphicsRenderObject | undefined,
    overpaint: Overpaint,
    lociApply: LociApply,
    clear: boolean,
    smoothing?: SmoothingContext,
  ) {
    if (!renderObject) return;

    const {
      tOverpaint,
      dOverpaint,
      uGroupCount,
      instanceCount,
      uInstanceGranularity: instanceGranularity,
    } = renderObject.values;
    const count = instanceGranularity.ref.value
      ? instanceCount.ref.value
      : uGroupCount.ref.value * instanceCount.ref.value;

    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    if (
      clear &&
      overpaint.layers.length > 0 &&
      dOverpaint.ref.value &&
      isKnownClear(renderObject, 'overpaint', tOverpaint, count, type) &&
      !reachesAny(overpaint.layers, lociApply)
    ) {
      return; // stays all clear
    }

    // ensure texture has right size and type
    createOverpaint(overpaint.layers.length ? count : 0, type, renderObject.values);
    const { array } = tOverpaint.ref.value;

    // clear all if requested
    if (clear) clearOverpaint(array, 0, count);
    const ranges: number[] = [];

    let touched = false;
    for (let i = 0, il = overpaint.layers.length; i < il; ++i) {
      const { loci, color, clear } = overpaint.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        ranges.push(start, end);
        return clear ? clearOverpaint(array, start, end) : applyOverpaintColor(array, start, end, color);
      };
      lociApply(loci, apply, false);
    }
    // cleared: all of it changed
    TextureChanges.update(tOverpaint, clear ? undefined : ranges);
    ValueCell.updateIfChanged(dOverpaint, overpaint.layers.length > 0);

    setKnownClear(renderObject, 'overpaint', tOverpaint, count, type, clear && !touched && overpaint.layers.length > 0);

    if (overpaint.layers.length === 0) return;
    if (type === 'instance') return;
    if (clear && !touched) return; // all clear: nothing to smooth (perf-audit 07 P1)

    if (smoothing && hasColorSmoothingProp(smoothing.props)) {
      const { geometry, props, graphics } = smoothing;
      if (geometry.kind === 'mesh') {
        const { resolution, overpaintTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyMeshOverpaintSmoothing(renderObject.values as any, csp, graphics, overpaintTexture);
          (geometry.meta as SurfaceMeta).overpaintTexture = renderObject.values.tOverpaintGrid.ref.value;
        }
      } else if (graphics && geometry.kind === 'texture-mesh') {
        const { resolution, overpaintTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyTextureMeshOverpaintSmoothing(renderObject.values as any, csp, graphics, overpaintTexture);
          (geometry.meta as SurfaceMeta).overpaintTexture = renderObject.values.tOverpaintGrid.ref.value;
        }
      }
    }
  }

  export function setTransparency(
    renderObject: GraphicsRenderObject | undefined,
    transparency: Transparency,
    lociApply: LociApply,
    clear: boolean,
    smoothing?: SmoothingContext,
  ) {
    if (!renderObject) return;

    const {
      tTransparency,
      transparencyAverage,
      transparencyMin,
      dTransparency,
      uGroupCount,
      instanceCount,
      uInstanceGranularity: instanceGranularity,
    } = renderObject.values;
    const count = instanceGranularity.ref.value
      ? instanceCount.ref.value
      : uGroupCount.ref.value * instanceCount.ref.value;

    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    if (
      clear &&
      transparency.layers.length > 0 &&
      dTransparency.ref.value &&
      isKnownClear(renderObject, 'transparency', tTransparency, count, type) &&
      !reachesAny(transparency.layers, lociApply)
    ) {
      return; // stays all clear
    }

    // ensure texture has right size and type
    createTransparency(transparency.layers.length ? count : 0, type, renderObject.values);
    const { array } = tTransparency.ref.value;

    // clear if requested
    if (clear) clearTransparency(array, 0, count);
    const ranges: number[] = [];

    let touched = false;
    for (let i = 0, il = transparency.layers.length; i < il; ++i) {
      const { loci, value } = transparency.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        ranges.push(start, end);
        return applyTransparencyValue(array, start, end, value);
      };
      lociApply(loci, apply, false);
    }
    // cleared: all of it changed
    TextureChanges.update(tTransparency, clear ? undefined : ranges);
    ValueCell.updateIfChanged(transparencyAverage, getTransparencyAverage(array, count));
    ValueCell.updateIfChanged(transparencyMin, getTransparencyMin(array, count));
    ValueCell.updateIfChanged(dTransparency, transparency.layers.length > 0);

    setKnownClear(
      renderObject,
      'transparency',
      tTransparency,
      count,
      type,
      clear && !touched && transparency.layers.length > 0,
    );

    if (transparency.layers.length === 0) return;
    if (type === 'instance') return;
    if (clear && !touched) return; // all clear: nothing to smooth (perf-audit 07 P1)

    if (smoothing && hasColorSmoothingProp(smoothing.props)) {
      const { geometry, props, graphics } = smoothing;
      if (geometry.kind === 'mesh') {
        const { resolution, transparencyTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyMeshTransparencySmoothing(renderObject.values as any, csp, graphics, transparencyTexture);
          (geometry.meta as SurfaceMeta).transparencyTexture = renderObject.values.tTransparencyGrid.ref.value;
        }
      } else if (graphics && geometry.kind === 'texture-mesh') {
        const { resolution, transparencyTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyTextureMeshTransparencySmoothing(renderObject.values as any, csp, graphics, transparencyTexture);
          (geometry.meta as SurfaceMeta).transparencyTexture = renderObject.values.tTransparencyGrid.ref.value;
        }
      }
    }
  }

  export function setEmissive(
    renderObject: GraphicsRenderObject | undefined,
    emissive: Emissive,
    lociApply: LociApply,
    clear: boolean,
    smoothing?: SmoothingContext,
  ) {
    if (!renderObject) return;

    const {
      tEmissive,
      emissiveAverage,
      dEmissive,
      uGroupCount,
      instanceCount,
      uInstanceGranularity: instanceGranularity,
    } = renderObject.values;
    const count = instanceGranularity.ref.value
      ? instanceCount.ref.value
      : uGroupCount.ref.value * instanceCount.ref.value;

    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    if (
      clear &&
      emissive.layers.length > 0 &&
      dEmissive.ref.value &&
      isKnownClear(renderObject, 'emissive', tEmissive, count, type) &&
      !reachesAny(emissive.layers, lociApply)
    ) {
      return; // stays all clear
    }

    // ensure texture has right size and type
    createEmissive(emissive.layers.length ? count : 0, type, renderObject.values);
    const { array } = tEmissive.ref.value;

    // clear if requested
    if (clear) clearEmissive(array, 0, count);
    const ranges: number[] = [];

    let touched = false;
    for (let i = 0, il = emissive.layers.length; i < il; ++i) {
      const { loci, value } = emissive.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        ranges.push(start, end);
        return applyEmissiveValue(array, start, end, value);
      };
      lociApply(loci, apply, false);
    }
    // cleared: all of it changed
    TextureChanges.update(tEmissive, clear ? undefined : ranges);
    ValueCell.updateIfChanged(emissiveAverage, getEmissiveAverage(array, count));
    ValueCell.updateIfChanged(dEmissive, emissive.layers.length > 0);

    setKnownClear(renderObject, 'emissive', tEmissive, count, type, clear && !touched && emissive.layers.length > 0);

    if (emissive.layers.length === 0) return;
    if (type === 'instance') return;
    if (clear && !touched) return; // all clear: nothing to smooth (perf-audit 07 P1)

    if (smoothing && hasColorSmoothingProp(smoothing.props)) {
      const { geometry, props, graphics } = smoothing;
      if (geometry.kind === 'mesh') {
        const { resolution, emissiveTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyMeshEmissiveSmoothing(renderObject.values as any, csp, graphics, emissiveTexture);
          (geometry.meta as SurfaceMeta).emissiveTexture = renderObject.values.tEmissiveGrid.ref.value;
        }
      } else if (graphics && geometry.kind === 'texture-mesh') {
        const { resolution, emissiveTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyTextureMeshEmissiveSmoothing(renderObject.values as any, csp, graphics, emissiveTexture);
          (geometry.meta as SurfaceMeta).emissiveTexture = renderObject.values.tEmissiveGrid.ref.value;
        }
      }
    }
  }

  export function setSubstance(
    renderObject: GraphicsRenderObject | undefined,
    substance: Substance,
    lociApply: LociApply,
    clear: boolean,
    smoothing?: SmoothingContext,
  ) {
    if (!renderObject) return;

    const {
      tSubstance,
      dSubstance,
      uGroupCount,
      instanceCount,
      uInstanceGranularity: instanceGranularity,
    } = renderObject.values;
    const count = instanceGranularity.ref.value
      ? instanceCount.ref.value
      : uGroupCount.ref.value * instanceCount.ref.value;

    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    if (
      clear &&
      substance.layers.length > 0 &&
      dSubstance.ref.value &&
      isKnownClear(renderObject, 'substance', tSubstance, count, type) &&
      !reachesAny(substance.layers, lociApply)
    ) {
      return; // stays all clear
    }

    // ensure texture has right size and type
    createSubstance(substance.layers.length ? count : 0, type, renderObject.values);
    const { array } = tSubstance.ref.value;

    // clear all if requested
    if (clear) clearSubstance(array, 0, count);
    const ranges: number[] = [];

    let touched = false;
    for (let i = 0, il = substance.layers.length; i < il; ++i) {
      const { loci, material, clear } = substance.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        ranges.push(start, end);
        return clear ? clearSubstance(array, start, end) : applySubstanceMaterial(array, start, end, material);
      };
      lociApply(loci, apply, false);
    }
    // cleared: all of it changed
    TextureChanges.update(tSubstance, clear ? undefined : ranges);
    ValueCell.updateIfChanged(dSubstance, substance.layers.length > 0);

    setKnownClear(renderObject, 'substance', tSubstance, count, type, clear && !touched && substance.layers.length > 0);

    if (substance.layers.length === 0) return;
    if (type === 'instance') return;
    if (clear && !touched) return; // all clear: nothing to smooth (perf-audit 07 P1)

    if (smoothing && hasColorSmoothingProp(smoothing.props)) {
      const { geometry, props, graphics } = smoothing;
      if (geometry.kind === 'mesh') {
        const { resolution, substanceTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyMeshSubstanceSmoothing(renderObject.values as any, csp, graphics, substanceTexture);
          (geometry.meta as SurfaceMeta).substanceTexture = renderObject.values.tSubstanceGrid.ref.value;
        }
      } else if (graphics && geometry.kind === 'texture-mesh') {
        const { resolution, substanceTexture } = geometry.meta as SurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, true, resolution);
        if (csp) {
          applyTextureMeshSubstanceSmoothing(renderObject.values as any, csp, graphics, substanceTexture);
          (geometry.meta as SurfaceMeta).substanceTexture = renderObject.values.tSubstanceGrid.ref.value;
        }
      }
    }
  }

  export function setClipping(
    renderObject: GraphicsRenderObject | undefined,
    clipping: Clipping,
    lociApply: LociApply,
    clear: boolean,
  ) {
    if (!renderObject) return;

    const {
      tClipping,
      dClipping,
      uGroupCount,
      instanceCount,
      uInstanceGranularity: instanceGranularity,
    } = renderObject.values;
    const count = instanceGranularity.ref.value
      ? instanceCount.ref.value
      : uGroupCount.ref.value * instanceCount.ref.value;
    const { layers } = clipping;

    // as the four layers above: a visual the change does not reach stays clear (perf-audit-2 07/P14)
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    if (
      clear &&
      layers.length > 0 &&
      dClipping.ref.value &&
      isKnownClear(renderObject, 'clipping', tClipping, count, type) &&
      !reachesAny(layers, lociApply)
    ) {
      return; // stays all clear
    }

    // ensure texture has right size and type
    createClipping(layers.length ? count : 0, type, renderObject.values);
    const { array } = tClipping.ref.value;

    // clear if requested
    if (clear) clearClipping(array, 0, count);
    const ranges: number[] = [];

    let touched = false;
    for (let i = 0, il = clipping.layers.length; i < il; ++i) {
      const { loci, groups } = clipping.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        ranges.push(start, end);
        return applyClippingGroups(array, start, end, groups);
      };
      lociApply(loci, apply, false);
    }
    // cleared: all of it changed
    TextureChanges.update(tClipping, clear ? undefined : ranges);
    ValueCell.updateIfChanged(dClipping, clipping.layers.length > 0);

    setKnownClear(renderObject, 'clipping', tClipping, count, type, clear && !touched && layers.length > 0);
  }

  export function setWiggle(
    renderObject: GraphicsRenderObject | undefined,
    wiggle: Wiggle,
    lociApply: LociApply,
    clear: boolean,
  ) {
    if (!renderObject) return;

    const {
      tWiggle,
      wiggleAverage,
      dWiggle,
      uGroupCount,
      instanceCount,
      uInstanceGranularity: instanceGranularity,
    } = renderObject.values;
    const count = instanceGranularity.ref.value
      ? instanceCount.ref.value
      : uGroupCount.ref.value * instanceCount.ref.value;

    // as the four layers above: a visual the change does not reach stays clear (perf-audit-2 07/P14)
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    if (
      clear &&
      wiggle.layers.length > 0 &&
      dWiggle.ref.value &&
      isKnownClear(renderObject, 'wiggle', tWiggle, count, type) &&
      !reachesAny(wiggle.layers, lociApply)
    ) {
      return; // stays all clear
    }

    // ensure texture has right size and type
    createWiggle(wiggle.layers.length ? count : 0, type, renderObject.values);
    const { array } = tWiggle.ref.value;

    // clear if requested
    if (clear) clearWiggle(array, 0, count);

    let touched = false;
    for (let i = 0, il = wiggle.layers.length; i < il; ++i) {
      const { loci, value } = wiggle.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        return applyWiggleValue(array, start, end, value);
      };
      lociApply(loci, apply, false);
    }
    ValueCell.update(tWiggle, tWiggle.ref.value);
    ValueCell.updateIfChanged(wiggleAverage, getWiggleAverage(array, count));
    ValueCell.updateIfChanged(dWiggle, wiggle.layers.length > 0);

    setKnownClear(renderObject, 'wiggle', tWiggle, count, type, clear && !touched && wiggle.layers.length > 0);
  }

  export function setThemeStrength(
    renderObject: GraphicsRenderObject | undefined,
    strength: { overpaint: number; transparency: number; emissive: number; substance: number; wiggle: number },
  ) {
    if (renderObject) {
      ValueCell.updateIfChanged(renderObject.values.uOverpaintStrength, strength.overpaint);
      ValueCell.updateIfChanged(renderObject.values.uTransparencyStrength, strength.transparency);
      ValueCell.updateIfChanged(renderObject.values.uEmissiveStrength, strength.emissive);
      ValueCell.updateIfChanged(renderObject.values.uSubstanceStrength, strength.substance);
      ValueCell.updateIfChanged(renderObject.values.uWiggleStrength, strength.wiggle);
    }
  }

  export function setTransform(
    renderObject: GraphicsRenderObject | undefined,
    transform?: Mat4,
    instanceTransforms?: Float32Array | null,
  ) {
    if (!renderObject || (!transform && !instanceTransforms)) return;

    const { values } = renderObject;
    if (transform) {
      Mat4.copy(values.matrix.ref.value, transform);
      ValueCell.update(values.matrix, values.matrix.ref.value);
    }
    if (instanceTransforms) {
      values.extraTransform.ref.value.set(instanceTransforms);
      ValueCell.update(values.extraTransform, values.extraTransform.ref.value);
      ValueCell.updateIfChanged(values.hasExtraTransform, true);
    } else if (instanceTransforms === null) {
      ValueCell.updateIfChanged(values.hasExtraTransform, false);
    }
    updateTransformData(
      values,
      values.invariantBoundingSphere.ref.value,
      values.instanceGrid.ref.value.cellSize,
      values.instanceGrid.ref.value.batchSize,
    );
    // from the transforms as drawn (`aTransform`: the matrix and any extra transform applied), as the geometries do
    const boundingSphere = calculateTransformBoundingSphere(
      values.invariantBoundingSphere.ref.value,
      values.aTransform.ref.value,
      values.instanceCount.ref.value,
      0,
    );
    ValueCell.update(values.boundingSphere, boundingSphere);
  }
}
