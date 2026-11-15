/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The renderer's solid interior caps: where the near plane or a clip object cuts a solid object, its interior marked in
// the stencil buffer (front and back faces counted) and filled, per cut and per segment it cuts of a merged object

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Clip } from '@molstar/core/util/clip';
import type { GraphicsContext } from '../context.js';
import { DrawFact, getDrawState } from '../draw-table.js';
import type { GlobalBindings, GraphicsRenderVariant } from '../render-item.js';
import type { RenderState } from '../render-state.js';
import type { GraphicsRenderable } from '../renderable.js';
import {
  type FactTest,
  Flag,
  getObjectFacts,
  getSolidInteriorFillStates,
  getSolidInteriorMarkState,
  type SolidInteriorMode,
} from './pass-draws.js';

/** What the caps read of the renderer: how it draws an object, and the view (getters: read when drawn) */
export interface SolidInteriorCapsHost {
  readonly ctx: GraphicsContext;
  readonly globals: GlobalBindings;
  /** the scaled bounding sphere of the object being drawn (`isVisible`, `setRowBoundingSphere`) */
  readonly boundingSphere: Sphere3D;
  readonly cameraPlane: Plane3D;
  readonly cameraPosition: Vec3;
  readonly modelScale: number;
  readonly near: number;
  /** objects whose second phase of GPU culling is to come */
  readonly cullSecondPending: number;
  /** draws `r` (visible, `boundingSphere` set) with `flag` in `passState`, its `segments` only where given */
  drawObject: (
    r: GraphicsRenderable,
    variant: GraphicsRenderVariant,
    flag: Flag,
    passState: RenderState,
    segments?: readonly number[],
  ) => void;
  /** whether `r` is drawn in `variant`, `boundingSphere` set */
  isVisible: (r: GraphicsRenderable, variant: GraphicsRenderVariant) => boolean;
}

export interface SolidInteriorCaps {
  /** caps drawn at all (fragment depth) */
  readonly supported: boolean;
  /** caps drawn in the open pass (its target with stencil) */
  canRender: () => boolean;
  /** of an object with the facts (`DrawFact.SolidInterior`) */
  has: (facts: number) => boolean;
  /** the caps `r` draws in this frame into `out` (-1 the near plane, else clip object indices) */
  get: (r: GraphicsRenderable, mode: SolidInteriorMode, out: number[]) => number[];
  /** per cap of the last `get` the segments it cuts, all where empty */
  readonly segments: readonly (readonly number[])[];
  /** draws the caps of `r`, visible with `boundingSphere` set */
  draw: (
    r: GraphicsRenderable,
    variant: GraphicsRenderVariant,
    mode: SolidInteriorMode,
    passState: RenderState,
  ) => void;
  /** draws the caps of the `renderables` passing `check`, those visible */
  drawAll: (
    renderables: ReadonlyArray<GraphicsRenderable>,
    check: FactTest,
    variant: GraphicsRenderVariant,
    mode: SolidInteriorMode,
    passState: RenderState,
  ) => void;
}

export function createSolidInteriorCaps(host: SolidInteriorCapsHost): SolidInteriorCaps {
  const { ctx, globals, boundingSphere, cameraPlane, cameraPosition, drawObject, isVisible } = host;

  const supported = ctx.capabilities.fragDepth;
  const canRenderSolidInteriorCaps = () => {
    return supported && ctx.hasRenderPassStencil();
  };
  /** of an object with the facts (`DrawFact.SolidInterior`) */
  const hasSolidInteriorCap = (facts: number) => {
    return supported && (facts & DrawFact.SolidInterior) !== 0;
  };

  const setSolidInteriorPass = (pass: number, clipIndex: number) => {
    globals.setSolidInteriorCap(pass, clipIndex);
  };

  const renderSolidInteriorMark = (
    r: GraphicsRenderable,
    variant: GraphicsRenderVariant,
    clipIndex: number,
    passState: RenderState,
    segments: readonly number[],
  ) => {
    ctx.clear({ stencil: 0 });
    setSolidInteriorPass(2, clipIndex);
    drawObject(r, variant, Flag.SolidInteriorMark, getSolidInteriorMarkState(passState), segments);
  };

  const renderSolidInteriorFill = (
    r: GraphicsRenderable,
    variant: GraphicsRenderVariant,
    mode: SolidInteriorMode,
    clipIndex: number,
    passState: RenderState,
    segments: readonly number[],
  ) => {
    setSolidInteriorPass(1, clipIndex);
    const { resolve, fill } = getSolidInteriorFillStates(passState, mode);
    if (resolve) {
      // the OIT depth attachment is otherwise unused: resolve the nearest back face, then only draw that one
      ctx.clear({ depth: 0 });
      drawObject(r, variant, Flag.SolidInteriorFill, resolve, segments);
    }
    drawObject(r, variant, Flag.SolidInteriorFill, fill, segments);
  };

  const solidInteriorEye = Vec3();
  const solidInteriorSphere = Sphere3D();

  /**
   * Whether near cuts the (scaled) sphere by more than a margin: the camera fits near to the front of the scene's
   * sphere, and an object spanning the scene has its sphere's front there too, its geometry behind near
   */
  const cutsNear = (s: Sphere3D, near: number) => {
    const d = Plane3D.distanceToPoint(cameraPlane, s.center);
    return d - s.radius < near * (1 - 1e-4) && d + s.radius >= near;
  };

  const cutsNearSegment = (s: Sphere3D) =>
    cutsNear(Sphere3D.scaleNX(solidInteriorSphere, s, host.modelScale), host.near);

  /**
   * Per cap of the last `getSolidInteriorCaps` (by position), the segments it cuts of an object drawing several (a
   * merged object, `Renderable.renderSegments`), drawn instead of every member the object draws; empty where it cuts
   * all or the object draws one
   */
  const capSegments: number[][] = [];

  /**
   * The segments `cuts` (by the segment's sphere, unscaled) into `capSegments[k]` (empty where all, or without
   * `spheres`); false where none
   */
  const setCapSegments = (
    spheres: readonly Sphere3D[] | undefined,
    k: number,
    cuts: (s: Sphere3D, objectIndex: number) => boolean,
    objectIndex: number,
  ) => {
    const segments = (capSegments[k] ??= []);
    segments.length = 0;
    if (!spheres) return true;
    for (let i = 0, il = spheres.length; i < il; ++i) if (cuts(spheres[i], objectIndex)) segments.push(i);
    if (segments.length === 0) return false;
    // whole also where it cuts all while a second phase of GPU culling is to come: its part draws no caps, and the
    // object's lists hold the first phase's instances only
    if (segments.length === spheres.length && host.cullSecondPending === 0) segments.length = 0;
    return true;
  };

  let clipObjects: Clip.Objects | undefined;
  const clipsSegment = (s: Sphere3D, objectIndex: number) => Clip.canIntersectSphere(clipObjects!, objectIndex, s);

  /**
   * The caps `r` draws in this frame into `out`: -1 for the near plane, the indices of the clip objects that cut it
   * (camera dependent), and per cap the segments it cuts (`capSegments`): a merged object's sphere spans its
   * members'; assumes `isVisible` returned true, which also set `boundingSphere`
   */
  const getSolidInteriorCaps = (r: GraphicsRenderable, mode: SolidInteriorMode, out: number[]) => {
    out.length = 0;
    const back = mode === 'back';
    const spheres = r.renderSegments ? r.getSegmentSpheres?.() : undefined;
    if (!back && cutsNear(boundingSphere, host.near)) {
      if (setCapSegments(spheres, 0, cutsNearSegment, -1)) out.push(-1);
    }

    const clipInfo = r.getClipInfo();
    if (clipInfo && clipInfo.capIndices.length > 0) {
      const { objects, planes, capIndices } = clipInfo;
      Vec3.scale(solidInteriorEye, cameraPosition, 1 / host.modelScale);
      clipObjects = objects;
      for (let j = 0, jl = capIndices.length; j < jl; ++j) {
        const i = capIndices[j];
        if (objects.type[i] === Clip.Type.plane) {
          const d = Plane3D.distanceToPoint(planes[i], solidInteriorEye) * host.modelScale;
          if (back ? d >= -1e-4 : d <= 1e-4) continue;
        }
        if (setCapSegments(spheres, out.length, clipsSegment, i)) out.push(i);
      }
    }
    return out;
  };

  const solidInteriorCapList: number[] = [];
  /** the caps of `r`; assumes it is visible, `boundingSphere` set (`isVisible`, `setRowBoundingSphere`) */
  const drawSolidInteriorCaps = (
    r: GraphicsRenderable,
    variant: GraphicsRenderVariant,
    mode: SolidInteriorMode,
    passState: RenderState,
  ) => {
    const caps = getSolidInteriorCaps(r, mode, solidInteriorCapList);
    for (let j = 0, jl = caps.length; j < jl; ++j) {
      renderSolidInteriorMark(r, variant, caps[j], passState, capSegments[j]);
      renderSolidInteriorFill(r, variant, mode, caps[j], passState, capSegments[j]);
    }
    setSolidInteriorPass(0, -1);
  };

  const renderSolidInteriorCaps = (
    renderables: ReadonlyArray<GraphicsRenderable>,
    check: FactTest,
    variant: GraphicsRenderVariant,
    mode: SolidInteriorMode,
    passState: RenderState,
  ) => {
    if (!supported) return;
    let hasCaps = false;
    for (let i = 0, il = renderables.length; i < il; ++i) {
      const r = renderables[i];
      const facts = getObjectFacts(r);
      if (!hasSolidInteriorCap(facts) || !check(facts, getDrawState(r.state))) continue;
      if (!hasCaps) {
        if (!canRenderSolidInteriorCaps()) return;
        hasCaps = true;
      }
      if (isVisible(r, variant)) drawSolidInteriorCaps(r, variant, mode, passState);
    }
  };

  return {
    supported,
    canRender: canRenderSolidInteriorCaps,
    has: hasSolidInteriorCap,
    get: getSolidInteriorCaps,
    segments: capSegments,
    draw: drawSolidInteriorCaps,
    drawAll: renderSolidInteriorCaps,
  };
}
