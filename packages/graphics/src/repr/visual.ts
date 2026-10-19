/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RuntimeContext } from '@molstar/core/task';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { Loci, isEmptyLoci, isEveryLoci, EveryLoci } from '@molstar/model/model/loci';
import {
  type MarkerAction,
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
  setAlphaFactor: (alphaFactor: number) => void;
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

  export function setAlphaFactor(renderObject: GraphicsRenderObject | undefined, alphaFactor: number) {
    if (renderObject) renderObject.state.alphaFactor = alphaFactor;
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

    let changed = false;
    let average = -1;
    let status: MarkerInfo['status'] = -1;
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
      lociApply(
        loci,
        (interval) => {
          size += Interval.size(interval);
          if (applyMarkerAction(array, interval, action)) changed = true;
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
      if (average === -1) {
        average = getMarkersAverage(array, count);
        if (average === 0) status = 0;
      }
      if (previous) {
        previous.action = action;
        previous.loci = loci;
        previous.status = currentStatus;
      }
      ValueCell.updateIfChanged(uMarker, status);
      if (status === -1) ValueCell.update(tMarker, tMarker.ref.value);
      ValueCell.updateIfChanged(markerAverage, average);
      ValueCell.updateIfChanged(markerStatus, status);
    }
    return changed;
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

    // ensure texture has right size and type
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    createOverpaint(overpaint.layers.length ? count : 0, type, renderObject.values);
    const { array } = tOverpaint.ref.value;

    // clear all if requested
    if (clear) clearOverpaint(array, 0, count);

    let touched = false;
    for (let i = 0, il = overpaint.layers.length; i < il; ++i) {
      const { loci, color, clear } = overpaint.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        return clear ? clearOverpaint(array, start, end) : applyOverpaintColor(array, start, end, color);
      };
      lociApply(loci, apply, false);
    }
    ValueCell.update(tOverpaint, tOverpaint.ref.value);
    ValueCell.updateIfChanged(dOverpaint, overpaint.layers.length > 0);

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

    // ensure texture has right size and type
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    createTransparency(transparency.layers.length ? count : 0, type, renderObject.values);
    const { array } = tTransparency.ref.value;

    // clear if requested
    if (clear) clearTransparency(array, 0, count);

    let touched = false;
    for (let i = 0, il = transparency.layers.length; i < il; ++i) {
      const { loci, value } = transparency.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        return applyTransparencyValue(array, start, end, value);
      };
      lociApply(loci, apply, false);
    }
    ValueCell.update(tTransparency, tTransparency.ref.value);
    ValueCell.updateIfChanged(transparencyAverage, getTransparencyAverage(array, count));
    ValueCell.updateIfChanged(transparencyMin, getTransparencyMin(array, count));
    ValueCell.updateIfChanged(dTransparency, transparency.layers.length > 0);

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

    // ensure texture has right size and type
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    createEmissive(emissive.layers.length ? count : 0, type, renderObject.values);
    const { array } = tEmissive.ref.value;

    // clear if requested
    if (clear) clearEmissive(array, 0, count);

    let touched = false;
    for (let i = 0, il = emissive.layers.length; i < il; ++i) {
      const { loci, value } = emissive.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        return applyEmissiveValue(array, start, end, value);
      };
      lociApply(loci, apply, false);
    }
    ValueCell.update(tEmissive, tEmissive.ref.value);
    ValueCell.updateIfChanged(emissiveAverage, getEmissiveAverage(array, count));
    ValueCell.updateIfChanged(dEmissive, emissive.layers.length > 0);

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

    // ensure texture has right size and type
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    createSubstance(substance.layers.length ? count : 0, type, renderObject.values);
    const { array } = tSubstance.ref.value;

    // clear all if requested
    if (clear) clearSubstance(array, 0, count);

    let touched = false;
    for (let i = 0, il = substance.layers.length; i < il; ++i) {
      const { loci, material, clear } = substance.layers[i];
      const apply = (interval: Interval) => {
        touched = true;
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        return clear ? clearSubstance(array, start, end) : applySubstanceMaterial(array, start, end, material);
      };
      lociApply(loci, apply, false);
    }
    ValueCell.update(tSubstance, tSubstance.ref.value);
    ValueCell.updateIfChanged(dSubstance, substance.layers.length > 0);

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

    // ensure texture has right size and type
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    createClipping(layers.length ? count : 0, type, renderObject.values);
    const { array } = tClipping.ref.value;

    // clear if requested
    if (clear) clearClipping(array, 0, count);

    for (let i = 0, il = clipping.layers.length; i < il; ++i) {
      const { loci, groups } = clipping.layers[i];
      const apply = (interval: Interval) => {
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        return applyClippingGroups(array, start, end, groups);
      };
      lociApply(loci, apply, false);
    }
    ValueCell.update(tClipping, tClipping.ref.value);
    ValueCell.updateIfChanged(dClipping, clipping.layers.length > 0);
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

    // ensure texture has right size and type
    const type = instanceGranularity.ref.value ? 'instance' : 'groupInstance';
    createWiggle(wiggle.layers.length ? count : 0, type, renderObject.values);
    const { array } = tWiggle.ref.value;

    // clear if requested
    if (clear) clearWiggle(array, 0, count);

    for (let i = 0, il = wiggle.layers.length; i < il; ++i) {
      const { loci, value } = wiggle.layers[i];
      const apply = (interval: Interval) => {
        const start = Interval.start(interval);
        const end = Interval.end(interval);
        return applyWiggleValue(array, start, end, value);
      };
      lociApply(loci, apply, false);
    }
    ValueCell.update(tWiggle, tWiggle.ref.value);
    ValueCell.updateIfChanged(wiggleAverage, getWiggleAverage(array, count));
    ValueCell.updateIfChanged(dWiggle, wiggle.layers.length > 0);
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
    const boundingSphere = calculateTransformBoundingSphere(
      values.invariantBoundingSphere.ref.value,
      values.transform.ref.value,
      values.instanceCount.ref.value,
      0,
    );
    ValueCell.update(values.boundingSphere, boundingSphere);
  }
}
