/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import type { ICamera } from '@molstar/graphics/canvas3d/camera';
import type { Scene } from './scene.js';
import type { GraphicsContext, PipelineRequest, TargetDescriptor } from './context.js';
import { Mat4, Vec3, Vec4, Vec2 } from '@molstar/core/math/linear-algebra';
import { type GraphicsRenderable, type Frame, createFrame } from './renderable.js';
import { Color } from '@molstar/core/util/color';
import { ValueCell, deepEqual } from '@molstar/core/util';
import { GlobalUniformSchema, type GlobalUniformValues, MaxLightCount } from './renderable/schema.js';
import type { GraphicsRenderVariant, ProgramItem } from './render-item.js';
import type { ObjectSignature } from './precompile.js';
import { getRenderableProgramSource } from './render-object.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { degToRad } from '@molstar/core/math/misc';
import type { Texture, Textures } from './texture.js';
import { arrayMapUpsert } from '@molstar/core/util/array';
import { clamp } from '@molstar/core/math/interpolate';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Clip } from '@molstar/core/util/clip';
import {
  BlendState,
  type CullMode,
  type DepthState,
  type FrontFace,
  RenderState,
  type StencilState,
} from './render-state.js';

export interface RendererStats {
  programCount: number;
  shaderCount: number;

  attributeCount: number;
  elementsCount: number;
  framebufferCount: number;
  renderbufferCount: number;
  textureCount: number;
  vertexArrayCount: number;

  drawCount: number;
  instanceCount: number;
  instancedDrawCount: number;
}

export enum PickType {
  None = 0,
  Object = 1,
  Instance = 2,
  Group = 3,
}

export enum MarkingType {
  None = 0,
  Depth = 1,
  Mask = 2,
}

interface Renderer {
  readonly stats: RendererStats;
  readonly props: Readonly<RendererProps>;
  readonly light: Readonly<Light>;
  readonly ambientColor: Vec3;

  clear: (toBackgroundColor: boolean, ignoreTransparentBackground?: boolean, forceToTransparency?: boolean) => void;
  clearDepth: (packed?: boolean) => void;
  update: (camera: ICamera, scene: Scene, frame: Frame) => void;
  setTime: (time: number) => void;

  renderPick: (group: Scene.Group, camera: ICamera, variant: 'pick' | 'depth', pickType: PickType) => void;
  renderDepth: (group: Scene.Group, camera: ICamera) => void;
  renderDepthOpaque: (group: Scene.Group, camera: ICamera) => void;
  renderDepthOpaqueBack: (group: Scene.Group, camera: ICamera) => void;
  renderDepthTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderMarkingDepth: (group: Scene.Group, camera: ICamera) => void;
  renderMarkingMask: (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => void;
  renderEmissiveOpaque: (
    group: Scene.Group,
    camera: ICamera,
    depthTexture: Texture,
    occludeWithOpaqueDepth: boolean,
  ) => void;
  renderEmissiveTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderTracing: (group: Scene.Group, camera: ICamera) => void;
  renderBlended: (group: Scene, camera: ICamera) => void;
  renderOpaque: (group: Scene.Group, camera: ICamera) => void;
  renderBlendedTransparent: (group: Scene.Group, camera: ICamera) => void;
  renderVolume: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderWboitTransparent: (group: Scene.Group, camera: ICamera, depthTexture: Texture) => void;
  renderDpoitTransparent: (
    group: Scene.Group,
    camera: ICamera,
    depthTexture: Texture,
    dpoitTextures: { depth: Texture; frontColor: Texture; backColor: Texture },
  ) => void;

  /**
   * The pipelines `pass` draws the objects of `group` with into a target like `target`, for
   * `resources.preparePipelines`: every object the pass would draw regardless of visibility and culling, and with `all`
   * also those that the objects' current state leaves out (e.g. unmarked objects for the marking passes)
   */
  collect: (
    pass: Renderer.Pass,
    group: Scene.Group,
    target: TargetDescriptor,
    options?: { all?: boolean },
  ) => PipelineRequest[];
  /**
   * The pipelines `pass` draws objects of the signatures with, through their program items
   * (`createSignatureProgramItem`): what drawing those objects would request, without them, for each way they may draw
   * (`ObjectSignature.opaque`, ...); `all` as for `collect`
   */
  collectSignatures: (
    pass: Renderer.Pass,
    items: readonly { readonly signature: ObjectSignature; readonly item: ProgramItem }[],
    target: TargetDescriptor,
    options?: { all?: boolean },
  ) => PipelineRequest[];

  setProps: (props: Partial<RendererProps>) => void;
  setViewport: (x: number, y: number, width: number, height: number) => void;
  setTransparentBackground: (value: boolean) => void;
  setDrawingBufferSize: (width: number, height: number) => void;
  setPixelRatio: (value: number) => void;
  setOcclusionTest: (f: ((s: Sphere3D) => boolean) | null) => void;

  dispose: () => void;
}

export const RendererParams = {
  backgroundColor: PD.Color(Color(0x000000), { description: 'Background color of the 3D canvas' }),

  pickingAlphaThreshold: PD.Numeric(
    0.5,
    { min: 0.0, max: 1.0, step: 0.01 },
    { description: 'The minimum opacity value needed for an object to be pickable.' },
  ),

  colorMarker: PD.Boolean(true, {
    description:
      'Tint marked and dim unmarked objects, by the marking pass when enabled, otherwise in the material color.',
  }),
  highlightColor: PD.Color(Color.fromNormalizedRgb(1.0, 0.4, 0.6)),
  selectColor: PD.Color(Color.fromNormalizedRgb(0.2, 1.0, 0.1)),
  dimColor: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
  highlightStrength: PD.Numeric(0.3, { min: 0.0, max: 1.0, step: 0.1 }),
  selectStrength: PD.Numeric(0.3, { min: 0.0, max: 1.0, step: 0.1 }),
  dimStrength: PD.Numeric(0.0, { min: 0.0, max: 1.0, step: 0.1 }),
  markerPriority: PD.Select(1, [
    [1, 'Highlight'],
    [2, 'Select'],
  ]),

  xrayEdgeFalloff: PD.Numeric(1, { min: 0.0, max: 3.0, step: 0.1 }),
  celSteps: PD.Numeric(5, { min: 2, max: 16, step: 1 }),
  exposure: PD.Numeric(1, { min: 0.0, max: 3.0, step: 0.01 }),

  light: PD.ObjectList(
    {
      inclination: PD.Numeric(150, { min: 0, max: 180, step: 1 }),
      azimuth: PD.Numeric(320, { min: 0, max: 360, step: 1 }),
      color: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
      intensity: PD.Numeric(0.6, { min: 0.0, max: 5.0, step: 0.01 }),
    },
    (o) => Color.toHexString(o.color),
    {
      defaultValue: [
        {
          inclination: 150,
          azimuth: 320,
          color: Color.fromNormalizedRgb(1.0, 1.0, 1.0),
          intensity: 0.6,
        },
      ],
    },
  ),
  ambientColor: PD.Color(Color.fromNormalizedRgb(1.0, 1.0, 1.0)),
  ambientIntensity: PD.Numeric(0.4, { min: 0.0, max: 2.0, step: 0.01 }),

  enableAnimation: PD.Boolean(true, { description: 'Enable time-based animations.' }),
};
export type RendererProps = PD.Values<typeof RendererParams>;

export type Light = {
  count: number;
  direction: number[];
  color: number[];
};

const tmpDir = Vec3();
const tmpColor = Vec3();
function getLight(props: RendererProps['light'], light?: Light): Light {
  const count = Math.min(props.length, MaxLightCount);
  const { direction, color } = light || {
    direction: new Array(count * 3).fill(0),
    color: new Array(count * 3).fill(0),
  };
  for (let i = 0; i < count; ++i) {
    const p = props[i];
    Vec3.directionFromSpherical(tmpDir, degToRad(p.inclination), degToRad(p.azimuth), 1);
    Vec3.toArray(tmpDir, direction, i * 3);
    Vec3.scale(tmpColor, Color.toVec3Normalized(tmpColor, p.color), p.intensity);
    Vec3.toArray(tmpColor, color, i * 3);
  }
  return { count, direction, color };
}

export function getTransformedLightDirection(light: Light, t: Mat4): Light['direction'] {
  const tld = new Array(light.count * 3);
  for (let i = 0, il = light.count; i < il; ++i) {
    Vec3.fromArray(tmpDir, light.direction, i * 3);
    Vec3.transformDirection(tmpDir, tmpDir, t);
    Vec3.toArray(tmpDir, tld, i * 3);
  }
  return tld;
}

namespace Renderer {
  const enum Flag {
    None = 0,
    BlendedFront = 1,
    BlendedBack = 2,
    DepthBack = 3,
    SolidInteriorMark = 4,
    SolidInteriorFill = 5,
  }

  type SolidInteriorMode = 'opaque' | 'back' | 'blended' | 'oit';

  const enum Mask {
    All = 0,
    Opaque = 1,
    Transparent = 2,
  }

  // depth is reversed (near 1, far 0, see `DepthRangeTag`): nearer is greater
  const DepthNearer: DepthState = { compare: 'greater', write: true };
  const DepthNearerReadOnly: DepthState = { compare: 'greater', write: false };

  const OpaqueState: RenderState = { depth: DepthNearer };
  const OpaqueBackState: RenderState = { depth: { compare: 'less', write: true } };
  /** test against caller-attached opaque depth, read-only + greater-equal (nearer or equal) so an emitter at its own depth still passes */
  const EmissiveOccludedState: RenderState = { depth: { compare: 'greater-equal', write: false } };
  /** max blending so overlapping faces don't accumulate */
  const EmissiveTransparentState: RenderState = { blend: BlendState.Max, depth: DepthNearerReadOnly };
  const BlendedTransparentState: RenderState = { blend: BlendState.Alpha, depth: DepthNearerReadOnly };
  const BlendedTransparentBackgroundState: RenderState = {
    blend: BlendState.Premultiplied,
    depth: DepthNearerReadOnly,
  };
  /** depth test done manually in shader against `depthTexture` */
  const VolumeState: RenderState = RenderState.Premultiplied;
  const WboitState: RenderState = {
    blend: BlendState(
      BlendState.Component('add', 'one', 'one'),
      BlendState.Component('add', 'zero', 'one-minus-src-alpha'),
    ),
  };
  const DpoitState: RenderState = { blend: BlendState.Max };

  const SolidInteriorMarkStencil: StencilState = {
    front: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'increment-wrap' },
    back: { compare: 'always', failOp: 'keep', depthFailOp: 'keep', passOp: 'decrement-wrap' },
    readMask: 0xff,
    writeMask: 0xff,
    reference: 0,
  };
  const SolidInteriorFillStencilKeep: StencilState = {
    front: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'keep' },
    back: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'keep' },
    readMask: 0xff,
    writeMask: 0xff,
    reference: 0,
  };
  const SolidInteriorFillStencilZero: StencilState = {
    front: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'zero' },
    back: { compare: 'not-equal', failOp: 'keep', depthFailOp: 'keep', passOp: 'zero' },
    readMask: 0xff,
    writeMask: 0xff,
    reference: 0,
  };
  const NoColorWrite = [false, false, false, false] as const;

  const CullModes: CullMode[] = ['none', 'front', 'back'];
  /** per-object cull and front face variants of the pass states */
  const FaceStates = new WeakMap<RenderState, RenderState[]>();
  function getFaceState(base: RenderState, cull: CullMode, frontFace: FrontFace) {
    let states = FaceStates.get(base);
    if (!states) {
      states = [];
      FaceStates.set(base, states);
    }
    const i = CullModes.indexOf(cull) * 2 + (frontFace === 'cw' ? 1 : 0);
    return (states[i] ??= { ...base, cull, frontFace });
  }

  /**
   * The render state of a draw of an object: the pass state with the object's cull and front face; `undefined` where
   * the object isn't drawn (direct volumes outside the `color` variant)
   */
  function getObjectState(
    values: GraphicsRenderable['values'],
    variant: GraphicsRenderVariant,
    flag: Flag,
    passState: RenderState,
  ): RenderState | undefined {
    const flipSided = !!values.dFlipSided?.ref.value;
    let cull: CullMode;
    let frontFace: FrontFace = flipSided ? 'cw' : 'ccw';
    if (values.dGeometryType.ref.value === 'directVolume') {
      if (variant !== 'color') {
        return; // only color supported
      }

      // culling done in fragment shader
      cull = 'none';
      frontFace = 'ccw';
    } else if (flag === Flag.BlendedFront) {
      cull = flipSided ? 'front' : 'back';
    } else if (
      flag === Flag.DepthBack ||
      flag === Flag.SolidInteriorMark ||
      (flag === Flag.SolidInteriorFill && values.hasReflection.ref.value)
    ) {
      // reflected instances flip the winding, so their exit faces are not back-facing
      cull = 'none';
    } else if (flag === Flag.BlendedBack || flag === Flag.SolidInteriorFill) {
      cull = flipSided ? 'back' : 'front';
    } else if (values.uDoubleSided && !values.uDoubleSided.ref.value && !values.hasReflection.ref.value) {
      cull = flipSided ? 'front' : 'back';
    } else {
      cull = 'none';
    }
    return getFaceState(passState, cull, frontFace);
  }

  /** marks the interior of a solid object in the stencil buffer, per pass state */
  const SolidInteriorMarkStates = new WeakMap<RenderState, RenderState>();
  function getSolidInteriorMarkState(passState: RenderState) {
    let state = SolidInteriorMarkStates.get(passState);
    if (!state) {
      state = { ...passState, colorWrite: NoColorWrite, depth: undefined, stencil: SolidInteriorMarkStencil };
      SolidInteriorMarkStates.set(passState, state);
    }
    return state;
  }

  type SolidInteriorFillStates = { readonly resolve?: RenderState; readonly fill: RenderState };
  /** fills the marked interior (`oit` resolves the nearest back face first), per pass state and mode */
  const SolidInteriorFillStatesCache = new WeakMap<RenderState, Map<SolidInteriorMode, SolidInteriorFillStates>>();
  function getSolidInteriorFillStates(passState: RenderState, mode: SolidInteriorMode): SolidInteriorFillStates {
    let cache = SolidInteriorFillStatesCache.get(passState);
    if (!cache) SolidInteriorFillStatesCache.set(passState, (cache = new Map()));
    let states = cache.get(mode);
    if (!states) {
      const writeDepth = mode === 'opaque' || mode === 'back';
      let resolve: RenderState | undefined;
      let depth: DepthState;
      if (mode === 'oit') {
        resolve = {
          ...passState,
          colorWrite: NoColorWrite,
          depth: { compare: 'greater-equal', write: true },
          stencil: SolidInteriorFillStencilKeep,
        };
        depth = { compare: 'equal', write: false };
      } else {
        depth = { compare: mode === 'back' ? 'less-equal' : 'greater-equal', write: writeDepth };
      }
      const fill: RenderState = {
        ...passState,
        depth,
        stencil: writeDepth ? SolidInteriorFillStencilKeep : SolidInteriorFillStencilZero,
      };
      states = { resolve, fill };
      cache.set(mode, states);
    }
    return states;
  }

  function getAlpha(r: GraphicsRenderable) {
    // uAlpha is updated in `r.render` so we need to recompute it here
    return clamp(r.values.alpha.ref.value * r.state.alphaFactor, 0, 1);
  }

  function isXrayShaded(r: GraphicsRenderable) {
    return r.values.dXrayShaded?.ref.value === 'on' || r.values.dXrayShaded?.ref.value === 'inverted';
  }

  function isDirectVolume(r: GraphicsRenderable) {
    return r.values.dGeometryType.ref.value === 'directVolume';
  }

  function checkOpaque(r: GraphicsRenderable) {
    return (
      (getAlpha(r) === 1 &&
        r.values.transparencyAverage.ref.value !== 1 &&
        !isDirectVolume(r) &&
        r.values.dPointStyle?.ref.value !== 'fuzzy' &&
        !isXrayShaded(r)) ||
      r.values.dTransparentBackfaces?.ref.value === 'opaque'
    );
  }

  function checkTransparent(r: GraphicsRenderable) {
    const alpha = getAlpha(r);
    return (
      alpha !== 0 &&
      (alpha < 1 ||
        r.values.transparencyAverage.ref.value > 0 ||
        isDirectVolume(r) ||
        r.values.dPointStyle?.ref.value === 'fuzzy' ||
        r.values.dGeometryType.ref.value === 'text' ||
        r.values.dGeometryType.ref.value === 'image' ||
        isXrayShaded(r))
    );
  }

  function checkPickable(r: GraphicsRenderable) {
    return !r.state.colorOnly;
  }

  function checkMarkingDepth(r: GraphicsRenderable) {
    return getAlpha(r) !== 0 && r.values.transparencyAverage.ref.value !== 1 && r.values.markerAverage.ref.value !== 1;
  }

  function checkMarkingMask(r: GraphicsRenderable) {
    return r.values.markerAverage.ref.value > 0;
  }

  function checkEmissive(r: GraphicsRenderable) {
    return r.values.emissiveAverage.ref.value + r.values.uEmissive.ref.value > 0;
  }

  /** caps are drawn from back faces, which transparent passes discard when `transparentBackfaces` is 'off' */
  function checkTransparentCap(r: GraphicsRenderable) {
    return checkTransparent(r) && r.values.dTransparentBackfaces?.ref.value !== 'off';
  }

  const checkAll = () => true;

  /** A draw of an object in a pass */
  type ObjectDraw = { readonly flag: Flag; readonly passState: RenderState };

  /** The draws of a renderer pass, shared by drawing (`renderPass`) and collecting its pipeline requests (`collect`) */
  interface PassDraws {
    readonly variant: GraphicsRenderVariant;
    /** objects the pass is for (by geometry), all if not given */
    readonly select?: (r: GraphicsRenderable) => boolean;
    /** objects the pass draws in their current state (opacity, markers, emission, ...) */
    readonly check: (r: GraphicsRenderable) => boolean;
    /** the draws of an object, in order */
    readonly draws: (r: GraphicsRenderable) => readonly ObjectDraw[];
    /** solid interior caps of the objects passing `check`, `inline` right after each object, else after all objects */
    readonly caps?: {
      readonly check: (r: GraphicsRenderable) => boolean;
      readonly mode: SolidInteriorMode;
      readonly passState: RenderState;
      readonly inline?: boolean;
    };
  }

  function singleDraw(passState: RenderState, flag = Flag.None): (r: GraphicsRenderable) => readonly ObjectDraw[] {
    const draws: readonly ObjectDraw[] = [{ flag, passState }];
    return () => draws;
  }

  /** draws objects of `check` with `passState`, caps of the same objects */
  function simplePassDraws(
    variant: GraphicsRenderVariant,
    check: (r: GraphicsRenderable) => boolean,
    passState: RenderState,
    caps: { check?: (r: GraphicsRenderable) => boolean; mode: SolidInteriorMode },
  ): PassDraws {
    return {
      variant,
      check,
      draws: singleDraw(passState),
      caps: { check: caps.check ?? check, mode: caps.mode, passState },
    };
  }

  function createEmissiveOpaquePassDraws(occludeWithOpaqueDepth: boolean): PassDraws {
    // without occluding opaque depth (packed depth): build own occluding depth from all opaque geometry
    const passState = occludeWithOpaqueDepth ? EmissiveOccludedState : OpaqueState;
    const check = (r: GraphicsRenderable) => checkOpaque(r) && (!occludeWithOpaqueDepth || checkEmissive(r));
    return simplePassDraws('emissive', check, passState, { mode: occludeWithOpaqueDepth ? 'blended' : 'opaque' });
  }

  function createBlendedTransparentPassDraws(transparentBackground: boolean): PassDraws {
    const passState = transparentBackground ? BlendedTransparentBackgroundState : BlendedTransparentState;
    const backAndFront: readonly ObjectDraw[] = [
      { flag: Flag.BlendedBack, passState },
      { flag: Flag.BlendedFront, passState },
    ];
    const front: readonly ObjectDraw[] = [{ flag: Flag.BlendedFront, passState }];
    const single: readonly ObjectDraw[] = [{ flag: Flag.None, passState }];
    return {
      variant: 'color',
      check: checkTransparent,
      // render frontfaces and backfaces separately to avoid artefacts
      draws: (r) => {
        if (!r.values.uDoubleSided?.ref.value) return single;
        return r.values.dTransparentBackfaces?.ref.value !== 'opaque' ? backAndFront : front;
      },
      caps: { check: checkTransparentCap, mode: 'blended', passState, inline: true },
    };
  }

  const PassDraws = {
    pick: simplePassDraws('pick', checkPickable, OpaqueState, { mode: 'opaque' }),
    pickDepth: simplePassDraws('depth', checkPickable, OpaqueState, { mode: 'opaque' }),
    depth: simplePassDraws('depth', checkAll, OpaqueState, { mode: 'opaque' }),
    depthOpaque: simplePassDraws('depth', checkOpaque, OpaqueState, { mode: 'opaque' }),
    depthOpaqueBack: {
      variant: 'depth',
      check: checkOpaque,
      draws: singleDraw(OpaqueBackState, Flag.DepthBack),
      caps: { check: checkOpaque, mode: 'back', passState: OpaqueBackState },
    },
    depthTransparent: simplePassDraws('depth', checkTransparent, OpaqueState, {
      check: checkTransparentCap,
      mode: 'opaque',
    }),
    markingDepth: simplePassDraws('marking', checkMarkingDepth, OpaqueState, { mode: 'opaque' }),
    markingMask: simplePassDraws('marking', checkMarkingMask, OpaqueState, { mode: 'opaque' }),
    emissiveOpaque: createEmissiveOpaquePassDraws(false),
    emissiveOpaqueOccluded: createEmissiveOpaquePassDraws(true),
    emissiveTransparent: {
      ...simplePassDraws('emissive', (r) => checkTransparent(r) && checkEmissive(r), EmissiveTransparentState, {
        check: (r) => checkTransparentCap(r) && checkEmissive(r),
        mode: 'blended',
      }),
      select: (r) => !isDirectVolume(r),
    },
    tracing: simplePassDraws('tracing', checkOpaque, OpaqueState, { mode: 'opaque' }),
    opaque: simplePassDraws('color', checkOpaque, OpaqueState, { mode: 'opaque' }),
    blendedTransparent: createBlendedTransparentPassDraws(false),
    blendedTransparentBackground: createBlendedTransparentPassDraws(true),
    volume: { variant: 'color', select: isDirectVolume, check: checkAll, draws: singleDraw(VolumeState) },
    wboitTransparent: simplePassDraws('color', checkTransparent, WboitState, {
      check: checkTransparentCap,
      mode: 'oit',
    }),
    dpoitTransparent: simplePassDraws('color', checkTransparent, DpoitState, {
      check: checkTransparentCap,
      mode: 'oit',
    }),
  } satisfies Record<string, PassDraws>;

  /** A pass of the renderer, as its `render*` methods draw it (`renderBlended` is `opaque` and `blendedTransparent`) */
  export type Pass =
    | { readonly name: 'pick'; readonly variant: 'pick' | 'depth' }
    | { readonly name: 'emissiveOpaque'; readonly occludeWithOpaqueDepth: boolean }
    | {
        readonly name:
          | 'depth'
          | 'depthOpaque'
          | 'depthOpaqueBack'
          | 'depthTransparent'
          | 'markingDepth'
          | 'markingMask'
          | 'emissiveTransparent'
          | 'tracing'
          | 'opaque'
          | 'blendedTransparent'
          | 'volume'
          | 'wboitTransparent'
          | 'dpoitTransparent';
      };

  function getPassDraws(pass: Pass, transparentBackground: boolean): PassDraws {
    switch (pass.name) {
      case 'pick':
        return pass.variant === 'pick' ? PassDraws.pick : PassDraws.pickDepth;
      case 'emissiveOpaque':
        return pass.occludeWithOpaqueDepth ? PassDraws.emissiveOpaqueOccluded : PassDraws.emissiveOpaque;
      case 'blendedTransparent':
        return transparentBackground ? PassDraws.blendedTransparentBackground : PassDraws.blendedTransparent;
      default:
        return PassDraws[pass.name];
    }
  }

  /** 1x1 depth texture cleared to the far plane */
  function createFarDepthTexture(ctx: GraphicsContext) {
    const texture = ctx.resources.texture('image-depth', 'depth', 'ushort', 'nearest');
    texture.define(1, 1);
    const framebuffer = ctx.resources.framebuffer();
    texture.attachFramebuffer(framebuffer, 'depth');
    ctx.beginRenderPass(framebuffer, { clear: { depth: 0 } });
    ctx.endRenderPass();
    texture.detachFramebuffer(framebuffer, 'depth');
    framebuffer.destroy();
    return texture;
  }

  export function create(ctx: GraphicsContext, props: Partial<RendererProps> = {}): Renderer {
    const { stats } = ctx;
    const p = PD.merge(RendererParams, PD.getDefaultValues(RendererParams), props);
    const light = getLight(p.light);

    const viewport = Viewport();
    const initialDrawingBufferSize = ctx.getDrawingBufferSize();
    const drawingBufferSize = Vec2.create(initialDrawingBufferSize.width, initialDrawingBufferSize.height);
    const bgColor = Color.toVec3Normalized(Vec3(), p.backgroundColor);

    let transparentBackground = false;
    let isOccluded: ((s: Sphere3D) => boolean) | null = null;

    // placeholders at the far plane (depth is reversed: 0, packed as zeros)
    const emptyDepthAlphaTexture = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    emptyDepthAlphaTexture.define(1, 1);
    emptyDepthAlphaTexture.load({ array: new Uint8Array([0, 0, 0, 0]), width: 1, height: 1 });
    const emptyDepthTexture = ctx.capabilities.depthTexture ? createFarDepthTexture(ctx) : emptyDepthAlphaTexture;
    const sharedTexturesList: Textures = [
      ['tDepth', emptyDepthTexture],
      ['tDepthAlpha', emptyDepthAlphaTexture],
    ];

    const model = Mat4();
    const view = Mat4();
    const invView = Mat4();
    const modelView = Mat4();
    const invModelView = Mat4();
    const invProjection = Mat4();
    const modelViewProjection = Mat4();
    const invModelViewProjection = Mat4();
    const invHeadRotation = Mat4();
    const modelViewEye = Mat4();
    const invModelViewEye = Mat4();

    const cameraDir = Vec3();
    const cameraPosition = Vec3();
    const cameraTarget = Vec3();
    const cameraPlane = Plane3D();
    const viewOffset = Vec2();
    const frustum = Frustum3D();

    let modelScale = 1;
    const boundingSphere = Sphere3D();
    let currentFrame: Frame = createFrame();

    const ambientColor = Vec3();
    Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);

    const globalUniforms: GlobalUniformValues = {
      uDrawId: ValueCell.create(0),

      uModel: ValueCell.create(Mat4.identity()),
      uView: ValueCell.create(view),
      uInvView: ValueCell.create(invView),
      uModelView: ValueCell.create(modelView),
      uInvModelView: ValueCell.create(invModelView),
      uInvProjection: ValueCell.create(invProjection),
      uProjection: ValueCell.create(Mat4()),
      uModelViewProjection: ValueCell.create(modelViewProjection),
      uInvModelViewProjection: ValueCell.create(invModelViewProjection),
      uHasHeadRotation: ValueCell.create(false),
      uInvHeadRotation: ValueCell.create(invHeadRotation),
      uHasEyeCamera: ValueCell.create(false),
      uModelViewEye: ValueCell.create(modelViewEye),
      uInvModelViewEye: ValueCell.create(invModelViewEye),
      uIsAsymmetricProjection: ValueCell.create(false),

      uIsOrtho: ValueCell.create(1),
      uViewOffset: ValueCell.create(viewOffset),
      uModelScale: ValueCell.create(1),

      uPixelRatio: ValueCell.create(ctx.pixelRatio),
      uViewport: ValueCell.create(Viewport.toVec4(Vec4(), viewport)),
      uDrawingBufferSize: ValueCell.create(drawingBufferSize),

      uCameraPosition: ValueCell.create(cameraPosition),
      uCameraDir: ValueCell.create(cameraDir),
      uCameraPlane: ValueCell.create(Plane3D.toArray(cameraPlane, Vec4(), 0)),
      uNear: ValueCell.create(1),
      uFar: ValueCell.create(10000),
      uFog: ValueCell.create(true),
      uFogNear: ValueCell.create(1),
      uFogFar: ValueCell.create(10000),
      uFogColor: ValueCell.create(bgColor),

      uRenderMask: ValueCell.create(0),
      uMarkingDepthTest: ValueCell.create(false),
      uDepthBack: ValueCell.create(false),
      uPickType: ValueCell.create(PickType.None),
      uMarkingType: ValueCell.create(MarkingType.None),
      uSolidInteriorPass: ValueCell.create(0),
      uSolidInteriorClip: ValueCell.create(-1),

      uTransparentBackground: ValueCell.create(false),

      uLightDirection: ValueCell.create(light.direction),
      uLightColor: ValueCell.create(light.color),
      uAmbientColor: ValueCell.create(ambientColor),

      uPickingAlphaThreshold: ValueCell.create(p.pickingAlphaThreshold),

      uHighlightColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.highlightColor)),
      uSelectColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.selectColor)),
      uDimColor: ValueCell.create(Color.toVec3Normalized(Vec3(), p.dimColor)),
      uHighlightStrength: ValueCell.create(p.highlightStrength),
      uSelectStrength: ValueCell.create(p.selectStrength),
      uDimStrength: ValueCell.create(p.dimStrength),
      uMarkerPriority: ValueCell.create(p.markerPriority),
      uMarkerAverage: ValueCell.create(0),

      uXrayEdgeFalloff: ValueCell.create(p.xrayEdgeFalloff),
      uCelSteps: ValueCell.create(p.celSteps),
      uExposure: ValueCell.create(p.exposure),

      uTime: ValueCell.create(0),
      uEnableAnimation: ValueCell.create(p.enableAnimation),
    };
    const globals = ctx.createGlobalBindings(GlobalUniformSchema, Object.entries(globalUniforms), sharedTexturesList);

    /** also sets `boundingSphere` to the scaled bounding sphere of `r` */
    const isVisible = (r: GraphicsRenderable, variant: GraphicsRenderVariant) => {
      if (r.state.disposed || !r.state.visible || (!r.state.pickable && variant === 'pick')) {
        return false;
      }

      if (!r.values.drawCount.ref.value) {
        return false;
      }

      Sphere3D.scaleNX(boundingSphere, r.values.boundingSphere.ref.value, modelScale);

      if (!Frustum3D.intersectsSphere3D(frustum, boundingSphere)) {
        return false;
      }

      const [minDistance, maxDistance] = r.values.uLod.ref.value;
      if (minDistance !== 0 || maxDistance !== 0) {
        const { center, radius } = boundingSphere;
        const d = Plane3D.distanceToPoint(cameraPlane, center);
        if (d + radius < minDistance * modelScale) return false;
        if (d - radius > maxDistance * modelScale) return false;
      }

      if (modelScale === 1 && isOccluded !== null && isOccluded(boundingSphere)) {
        return false;
      }

      return true;
    };

    /** assumes `isVisible` returned true, which also set `boundingSphere` */
    const drawObject = (r: GraphicsRenderable, variant: GraphicsRenderVariant, flag: Flag, passState: RenderState) => {
      const unscaled = modelScale === 1;
      if (unscaled) {
        const hasInstanceGrid = r.values.instanceGrid.ref.value.cellSize > 0;
        const hasMultipleInstances = r.values.uInstanceCount.ref.value > 1;
        if (hasInstanceGrid && (hasMultipleInstances || r.values.lodLevels)) {
          r.cull(cameraPlane, frustum, isOccluded, ctx.stats, currentFrame);
        } else {
          r.uncull();
        }
      } else {
        if (r.values.lodLevels) {
          const { center, radius } = boundingSphere;
          const d = Plane3D.distanceToPoint(cameraPlane, center);
          r.cullSimple(d, radius, modelScale);
        } else {
          r.uncull();
        }
      }

      const state = getObjectState(r.values, variant, flag, passState);
      if (!state) return;
      ctx.setRenderState(state);

      r.render(variant, globals);
    };

    const renderObject = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      flag: Flag,
      passState: RenderState,
    ) => {
      if (isVisible(r, variant)) drawObject(r, variant, flag, passState);
    };

    const solidInteriorCapSupported = ctx.capabilities.fragDepth;
    const canRenderSolidInteriorCaps = () => {
      return solidInteriorCapSupported && ctx.hasRenderPassStencil();
    };
    const hasSolidInteriorCap = (r: GraphicsRenderable) => {
      const geomType = r.values.dGeometryType.ref.value;
      return (
        solidInteriorCapSupported &&
        (geomType === 'mesh' || geomType === 'textureMesh') &&
        !!r.values.dSolidInterior?.ref.value
      );
    };

    const setSolidInteriorPass = (pass: number, clipIndex: number) => {
      ValueCell.updateIfChanged(globalUniforms.uSolidInteriorPass, pass);
      ValueCell.updateIfChanged(globalUniforms.uSolidInteriorClip, clipIndex);
      globals.invalidate();
    };

    const renderSolidInteriorMark = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      clipIndex: number,
      passState: RenderState,
    ) => {
      ctx.clear({ stencil: 0 });
      setSolidInteriorPass(2, clipIndex);
      drawObject(r, variant, Flag.SolidInteriorMark, getSolidInteriorMarkState(passState));
    };

    const renderSolidInteriorFill = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      mode: SolidInteriorMode,
      clipIndex: number,
      passState: RenderState,
    ) => {
      setSolidInteriorPass(1, clipIndex);
      const { resolve, fill } = getSolidInteriorFillStates(passState, mode);
      if (resolve) {
        // the OIT depth attachment is otherwise unused: resolve the nearest back face, then only draw that one
        ctx.clear({ depth: 0 });
        drawObject(r, variant, Flag.SolidInteriorFill, resolve);
      }
      drawObject(r, variant, Flag.SolidInteriorFill, fill);
    };

    const solidInteriorEye = Vec3();

    const renderSolidInteriorCap = (
      r: GraphicsRenderable,
      variant: GraphicsRenderVariant,
      mode: SolidInteriorMode,
      passState: RenderState,
    ) => {
      if (!isVisible(r, variant)) return;

      const back = mode === 'back';
      const near = globalUniforms.uNear.ref.value * 1.0001;
      if (
        !back &&
        Math.abs(Plane3D.distanceToPoint(cameraPlane, boundingSphere.center) - near) <= boundingSphere.radius
      ) {
        renderSolidInteriorMark(r, variant, -1, passState);
        renderSolidInteriorFill(r, variant, mode, -1, passState);
      }

      const clipInfo = r.getClipInfo();
      if (clipInfo && clipInfo.capIndices.length > 0) {
        const { objects, planes, capIndices } = clipInfo;
        Vec3.scale(solidInteriorEye, cameraPosition, 1 / modelScale);
        for (let j = 0, jl = capIndices.length; j < jl; ++j) {
          const i = capIndices[j];
          if (objects.type[i] === Clip.Type.plane) {
            const d = Plane3D.distanceToPoint(planes[i], solidInteriorEye) * modelScale;
            if (back ? d >= -1e-4 : d <= 1e-4) continue;
          }
          renderSolidInteriorMark(r, variant, i, passState);
          renderSolidInteriorFill(r, variant, mode, i, passState);
        }
      }
      setSolidInteriorPass(0, -1);
    };

    const renderSolidInteriorCaps = (
      renderables: ReadonlyArray<GraphicsRenderable>,
      check: (r: GraphicsRenderable) => boolean,
      variant: GraphicsRenderVariant,
      mode: SolidInteriorMode,
      passState: RenderState,
    ) => {
      if (!solidInteriorCapSupported) return;
      let hasCaps = false;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if (!hasSolidInteriorCap(r) || !check(r)) continue;
        if (!hasCaps) {
          if (!canRenderSolidInteriorCaps()) return;
          hasCaps = true;
        }
        renderSolidInteriorCap(r, variant, mode, passState);
      }
    };

    const update = (camera: ICamera, scene: Scene, frame: Frame) => {
      currentFrame = frame;

      ValueCell.update(globalUniforms.uView, camera.view);
      ValueCell.update(globalUniforms.uInvView, Mat4.invert(invView, camera.view));
      ValueCell.update(globalUniforms.uProjection, camera.projection);
      ValueCell.update(globalUniforms.uInvProjection, Mat4.invert(invProjection, camera.projection));

      ValueCell.updateIfChanged(globalUniforms.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
      ValueCell.update(
        globalUniforms.uViewOffset,
        camera.viewOffset.enabled
          ? Vec2.set(viewOffset, camera.viewOffset.offsetX * 16, camera.viewOffset.offsetY * 16)
          : Vec2.set(viewOffset, 0, 0),
      );
      ValueCell.updateIfChanged(globalUniforms.uModelScale, camera.scale);

      ValueCell.update(globalUniforms.uCameraPosition, Mat4.getTranslation(cameraPosition, invView));
      Vec3.scale(cameraTarget, camera.state.target, camera.scale);
      Vec3.normalize(cameraDir, Vec3.sub(cameraDir, cameraTarget, cameraPosition));
      ValueCell.update(globalUniforms.uCameraDir, cameraDir);

      ValueCell.updateIfChanged(globalUniforms.uFar, camera.far);
      ValueCell.updateIfChanged(globalUniforms.uNear, camera.near);
      ValueCell.updateIfChanged(globalUniforms.uFog, camera.state.fog > 0);
      ValueCell.updateIfChanged(globalUniforms.uFogFar, camera.fogFar);
      ValueCell.updateIfChanged(globalUniforms.uFogNear, camera.fogNear);
      ValueCell.updateIfChanged(globalUniforms.uTransparentBackground, transparentBackground);

      Frustum3D.fromProjectionMatrix(frustum, camera.projectionView);

      Plane3D.copy(cameraPlane, frustum[Frustum3D.PlaneIndex.Near]);
      cameraPlane.constant -= Plane3D.distanceToPoint(cameraPlane, cameraPosition);
      ValueCell.update(
        globalUniforms.uCameraPlane,
        Plane3D.toArray(cameraPlane, globalUniforms.uCameraPlane.ref.value, 0),
      );

      ValueCell.updateIfChanged(globalUniforms.uMarkerAverage, scene.markerAverage);

      const hasHeadRotation = !Mat4.isZero(camera.headRotation);
      if (hasHeadRotation) {
        ValueCell.updateIfChanged(globalUniforms.uHasHeadRotation, true);
        ValueCell.update(globalUniforms.uInvHeadRotation, Mat4.invert(invHeadRotation, camera.headRotation));
        ValueCell.update(globalUniforms.uLightDirection, getTransformedLightDirection(light, invHeadRotation));
      } else {
        ValueCell.updateIfChanged(globalUniforms.uHasHeadRotation, false);
        ValueCell.updateIfChanged(globalUniforms.uInvHeadRotation, Mat4.id);
        ValueCell.update(globalUniforms.uLightDirection, light.direction);
      }

      ValueCell.update(globalUniforms.uIsAsymmetricProjection, camera.isAsymmetricProjection);
    };

    const updateInternal = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture | null,
      depthAlphaTexture: Texture | null,
      renderMask: Mask,
      markingDepthTest: boolean,
    ) => {
      arrayMapUpsert(sharedTexturesList, 'tDepth', depthTexture || emptyDepthTexture);
      arrayMapUpsert(sharedTexturesList, 'tDepthAlpha', depthAlphaTexture || emptyDepthAlphaTexture);

      modelScale = camera.scale;

      ValueCell.update(globalUniforms.uModel, Mat4.scaleUniformly(model, group.view, modelScale));
      ValueCell.update(globalUniforms.uModelView, Mat4.mul(modelView, camera.view, model));
      ValueCell.update(globalUniforms.uInvModelView, Mat4.invert(invModelView, modelView));
      ValueCell.update(
        globalUniforms.uModelViewProjection,
        Mat4.mul(modelViewProjection, modelView, camera.projection),
      );
      ValueCell.update(
        globalUniforms.uInvModelViewProjection,
        Mat4.invert(invModelViewProjection, modelViewProjection),
      );

      ValueCell.updateIfChanged(globalUniforms.uRenderMask, renderMask);
      ValueCell.updateIfChanged(globalUniforms.uMarkingDepthTest, markingDepthTest);

      const hasEyeCamera = !Mat4.isZero(camera.viewEye);
      if (hasEyeCamera) {
        ValueCell.updateIfChanged(globalUniforms.uHasEyeCamera, true);
        ValueCell.update(globalUniforms.uModelViewEye, Mat4.mul(modelViewEye, camera.viewEye, model));
        ValueCell.update(globalUniforms.uInvModelViewEye, Mat4.invert(invModelViewEye, modelViewEye));
      } else {
        ValueCell.updateIfChanged(globalUniforms.uHasEyeCamera, false);
        ValueCell.updateIfChanged(globalUniforms.uModelViewEye, Mat4.id);
        ValueCell.updateIfChanged(globalUniforms.uInvModelViewEye, Mat4.id);
      }

      const { x, y, width, height } = viewport;
      ctx.setViewport(x, y, width, height);
      ctx.setScissor(x, y, width, height);

      globals.invalidate();
    };

    /** draws the objects of `group` as `pass` describes */
    const renderPass = (pass: PassDraws, group: Scene.Group) => {
      const { variant, select, check, draws, caps } = pass;
      const inlineCaps = !!caps?.inline && canRenderSolidInteriorCaps();
      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        if ((select && !select(r)) || !check(r)) continue;
        const objectDraws = draws(r);
        for (let j = 0, jl = objectDraws.length; j < jl; ++j) {
          renderObject(r, variant, objectDraws[j].flag, objectDraws[j].passState);
        }
        if (inlineCaps && hasSolidInteriorCap(r) && caps!.check(r)) {
          renderSolidInteriorCap(r, variant, caps!.mode, caps!.passState);
        }
      }
      if (caps && !caps.inline) renderSolidInteriorCaps(renderables, caps.check, variant, caps.mode, caps.passState);
    };

    const collect = (pass: Pass, group: Scene.Group, target: TargetDescriptor, options?: { all?: boolean }) => {
      return collectFrom(pass, group.renderables, target, !!options?.all);
    };

    /**
     * stand-ins of the signatures' objects, with the values the pass draws and object states read and those the passes
     * check (alpha, transparency, emission): one per way the object may draw (`ObjectSignature.opaque`, ...), so that the
     * passes' own checks select what to request
     */
    const signatureStandIns = new WeakMap<ProgramItem, Map<string, GraphicsRenderable[]>>();
    const getSignatureStandIns = (signature: ObjectSignature, item: ProgramItem) => {
      const { opaque, transparent, emissive } = signature;
      const flags = `${opaque}|${transparent}|${emissive}`;
      let cache = signatureStandIns.get(item);
      if (!cache) signatureStandIns.set(item, (cache = new Map()));
      let standIns = cache.get(flags);
      if (!standIns) {
        const { defines } = signature;
        const fixed: { [k: string]: unknown } = {
          dGeometryType: defines.dGeometryType,
          dSolidInterior: !!defines.dSolidInterior,
          dFlipSided: !!defines.dFlipSided,
          dTransparentBackfaces: defines.dTransparentBackfaces,
          dXrayShaded: defines.dXrayShaded,
          dPointStyle: defines.dPointStyle,
          uDoubleSided: signature.doubleSided,
          hasReflection: signature.hasReflection,
          transparencyAverage: 0,
          markerAverage: 0,
          emissiveAverage: 0,
          uEmissive: emissive ? 1 : 0,
        };
        // only the values of the kind: the object state reads a missing one apart from a false one (`uDoubleSided`)
        const schema = getRenderableProgramSource(signature)?.schema;
        const alphas = [...(opaque ? [1] : []), ...(transparent ? [0.5] : [])];
        standIns = (alphas.length ? alphas : [0]).map((alpha) => {
          const values: { [k: string]: ValueCell<unknown> } = { alpha: ValueCell.create(alpha) };
          for (const k of Object.keys(fixed)) if (!schema || schema[k]) values[k] = ValueCell.create(fixed[k]);
          return {
            values,
            renderItem: item,
            state: { disposed: false, alphaFactor: 1, colorOnly: false },
          } as unknown as GraphicsRenderable;
        });
        cache.set(flags, standIns);
      }
      return standIns;
    };

    const collectSignatures = (
      pass: Pass,
      items: readonly { readonly signature: ObjectSignature; readonly item: ProgramItem }[],
      target: TargetDescriptor,
      options?: { all?: boolean },
    ) => {
      const standIns = items.flatMap(({ signature, item }) => getSignatureStandIns(signature, item));
      return collectFrom(pass, standIns, target, !!options?.all);
    };

    const collectFrom = (
      pass: Pass,
      renderables: readonly GraphicsRenderable[],
      target: TargetDescriptor,
      all: boolean,
    ) => {
      const requests: PipelineRequest[] = [];
      const { variant, select, check, draws, caps } = getPassDraws(pass, transparentBackground);
      const capsDrawn = solidInteriorCapSupported && target.hasStencil;
      const request = (r: GraphicsRenderable, flag: Flag, passState: RenderState) => {
        const state = getObjectState(r.values, variant, flag, passState);
        if (state) requests.push({ item: r.renderItem, variant, state, target });
      };
      for (const r of renderables) {
        if (r.state.disposed || (select && !select(r)) || (!all && !check(r))) continue;
        for (const { flag, passState } of draws(r)) request(r, flag, passState);
        if (caps && capsDrawn && hasSolidInteriorCap(r) && (all || caps.check(r))) {
          request(r, Flag.SolidInteriorMark, getSolidInteriorMarkState(caps.passState));
          const { resolve, fill } = getSolidInteriorFillStates(caps.passState, caps.mode);
          if (resolve) request(r, Flag.SolidInteriorFill, resolve);
          request(r, Flag.SolidInteriorFill, fill);
        }
      }
      return requests;
    };

    const renderPick = (group: Scene.Group, camera: ICamera, variant: 'pick' | 'depth', pickType: PickType) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uPickType, pickType);
      renderPass(variant === 'pick' ? PassDraws.pick : PassDraws.pickDepth, group);
    };

    const renderDepth = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      renderPass(PassDraws.depth, group);
    };

    const renderDepthOpaque = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      renderPass(PassDraws.depthOpaque, group);
    };

    const renderDepthOpaqueBack = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, true);
      renderPass(PassDraws.depthOpaqueBack, group);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, false);
    };

    const renderDepthTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      renderPass(PassDraws.depthTransparent, group);
    };

    const renderMarkingDepth = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Depth);
      renderPass(PassDraws.markingDepth, group);
    };

    const renderMarkingMask = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, null, depthTexture, Mask.All, !!depthTexture);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Mask);
      renderPass(PassDraws.markingMask, group);
    };

    const renderEmissiveOpaque = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture,
      occludeWithOpaqueDepth: boolean,
    ) => {
      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Opaque, false);
      renderPass(occludeWithOpaqueDepth ? PassDraws.emissiveOpaqueOccluded : PassDraws.emissiveOpaque, group);
    };

    const renderEmissiveTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Transparent, false);
      renderPass(PassDraws.emissiveTransparent, group);
    };

    const renderTracing = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      renderPass(PassDraws.tracing, group);
    };

    const renderBlended = (scene: Scene, camera: ICamera) => {
      if (scene.hasOpaque) {
        renderOpaque(scene, camera);
      }
      if (scene.opacityAverage < 1) {
        renderBlendedTransparent(scene, camera);
      }
    };

    const renderOpaque = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      renderPass(PassDraws.opaque, group);
    };

    const renderBlendedTransparent = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Transparent, false);
      renderPass(transparentBackground ? PassDraws.blendedTransparentBackground : PassDraws.blendedTransparent, group);
    };

    const renderVolume = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      renderPass(PassDraws.volume, group);
    };

    const renderWboitTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      renderPass(PassDraws.wboitTransparent, group);
    };

    const renderDpoitTransparent = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture,
      dpoitTextures: { depth: Texture; frontColor: Texture; backColor: Texture },
    ) => {
      arrayMapUpsert(sharedTexturesList, 'tDpoitDepth', dpoitTextures.depth);
      arrayMapUpsert(sharedTexturesList, 'tDpoitFrontColor', dpoitTextures.frontColor);
      arrayMapUpsert(sharedTexturesList, 'tDpoitBackColor', dpoitTextures.backColor);

      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      renderPass(PassDraws.dpoitTransparent, group);
    };

    return {
      clear: (toBackgroundColor: boolean, ignoreTransparentBackground?: boolean, forceToTransparency?: boolean) => {
        ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
        if (forceToTransparency || (transparentBackground && !ignoreTransparentBackground)) {
          ctx.clear({ color: [0, 0, 0, 0], depth: 0 });
        } else if (toBackgroundColor) {
          ctx.clear({ color: [bgColor[0], bgColor[1], bgColor[2], 1], depth: 0 });
        } else {
          ctx.clear({ color: [1, 1, 1, 1], depth: 0 });
        }
      },
      clearDepth: (packed = false) => {
        ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
        // the far plane: depth 0, packed as zeros
        if (packed) {
          ctx.clear({ color: [0, 0, 0, 0], depth: 0 });
        } else {
          ctx.clear({ depth: 0 });
        }
      },
      update,

      renderPick,
      renderDepth,
      renderDepthOpaque,
      renderDepthOpaqueBack,
      renderDepthTransparent,
      renderMarkingDepth,
      renderMarkingMask,
      renderEmissiveOpaque,
      renderEmissiveTransparent,
      renderTracing,
      renderBlended,
      renderOpaque,
      renderBlendedTransparent,
      renderVolume,
      renderWboitTransparent,
      renderDpoitTransparent,

      collect,
      collectSignatures,

      setTime: (time: number) => {
        ValueCell.updateIfChanged(globalUniforms.uTime, time);
      },
      setProps: (props: Partial<RendererProps>) => {
        if (props.backgroundColor !== undefined && props.backgroundColor !== p.backgroundColor) {
          p.backgroundColor = props.backgroundColor;
          Color.toVec3Normalized(bgColor, p.backgroundColor);
          ValueCell.update(globalUniforms.uFogColor, Vec3.copy(globalUniforms.uFogColor.ref.value, bgColor));
        }

        if (props.pickingAlphaThreshold !== undefined && props.pickingAlphaThreshold !== p.pickingAlphaThreshold) {
          p.pickingAlphaThreshold = props.pickingAlphaThreshold;
          ValueCell.update(globalUniforms.uPickingAlphaThreshold, p.pickingAlphaThreshold);
        }

        if (props.colorMarker !== undefined && props.colorMarker !== p.colorMarker) {
          p.colorMarker = props.colorMarker;
        }
        if (props.highlightColor !== undefined && props.highlightColor !== p.highlightColor) {
          p.highlightColor = props.highlightColor;
          ValueCell.update(
            globalUniforms.uHighlightColor,
            Color.toVec3Normalized(globalUniforms.uHighlightColor.ref.value, p.highlightColor),
          );
        }
        if (props.selectColor !== undefined && props.selectColor !== p.selectColor) {
          p.selectColor = props.selectColor;
          ValueCell.update(
            globalUniforms.uSelectColor,
            Color.toVec3Normalized(globalUniforms.uSelectColor.ref.value, p.selectColor),
          );
        }
        if (props.dimColor !== undefined && props.dimColor !== p.dimColor) {
          p.dimColor = props.dimColor;
          ValueCell.update(
            globalUniforms.uDimColor,
            Color.toVec3Normalized(globalUniforms.uDimColor.ref.value, p.dimColor),
          );
        }
        if (props.highlightStrength !== undefined && props.highlightStrength !== p.highlightStrength) {
          p.highlightStrength = props.highlightStrength;
          ValueCell.update(globalUniforms.uHighlightStrength, p.highlightStrength);
        }
        if (props.selectStrength !== undefined && props.selectStrength !== p.selectStrength) {
          p.selectStrength = props.selectStrength;
          ValueCell.update(globalUniforms.uSelectStrength, p.selectStrength);
        }
        if (props.dimStrength !== undefined && props.dimStrength !== p.dimStrength) {
          p.dimStrength = props.dimStrength;
          ValueCell.update(globalUniforms.uDimStrength, p.dimStrength);
        }
        if (props.markerPriority !== undefined && props.markerPriority !== p.markerPriority) {
          p.markerPriority = props.markerPriority;
          ValueCell.update(globalUniforms.uMarkerPriority, p.markerPriority);
        }

        if (props.xrayEdgeFalloff !== undefined && props.xrayEdgeFalloff !== p.xrayEdgeFalloff) {
          p.xrayEdgeFalloff = props.xrayEdgeFalloff;
          ValueCell.update(globalUniforms.uXrayEdgeFalloff, p.xrayEdgeFalloff);
        }

        if (props.celSteps !== undefined && props.celSteps !== p.celSteps) {
          p.celSteps = props.celSteps;
          ValueCell.update(globalUniforms.uCelSteps, p.celSteps);
        }

        if (props.exposure !== undefined && props.exposure !== p.exposure) {
          p.exposure = props.exposure;
          ValueCell.update(globalUniforms.uExposure, p.exposure);
        }

        if (props.light !== undefined && !deepEqual(props.light, p.light)) {
          p.light = props.light;
          Object.assign(light, getLight(props.light, light));
          ValueCell.update(globalUniforms.uLightDirection, light.direction);
          ValueCell.update(globalUniforms.uLightColor, light.color);
        }
        if (props.ambientColor !== undefined && props.ambientColor !== p.ambientColor) {
          p.ambientColor = props.ambientColor;
          Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);
          ValueCell.update(globalUniforms.uAmbientColor, ambientColor);
        }
        if (props.ambientIntensity !== undefined && props.ambientIntensity !== p.ambientIntensity) {
          p.ambientIntensity = props.ambientIntensity;
          Vec3.scale(ambientColor, Color.toArrayNormalized(p.ambientColor, ambientColor, 0), p.ambientIntensity);
          ValueCell.update(globalUniforms.uAmbientColor, ambientColor);
        }

        if (props.enableAnimation !== undefined && props.enableAnimation !== p.enableAnimation) {
          p.enableAnimation = props.enableAnimation;
          ValueCell.update(globalUniforms.uEnableAnimation, p.enableAnimation);
        }
      },
      setViewport: (x: number, y: number, width: number, height: number) => {
        ctx.setViewport(x, y, width, height);
        ctx.setScissor(x, y, width, height);
        if (x !== viewport.x || y !== viewport.y || width !== viewport.width || height !== viewport.height) {
          Viewport.set(viewport, x, y, width, height);
          ValueCell.update(globalUniforms.uViewport, Vec4.set(globalUniforms.uViewport.ref.value, x, y, width, height));
        }
      },
      setTransparentBackground: (value: boolean) => {
        transparentBackground = value;
      },
      setDrawingBufferSize: (width: number, height: number) => {
        if (width !== drawingBufferSize[0] || height !== drawingBufferSize[1]) {
          ValueCell.update(globalUniforms.uDrawingBufferSize, Vec2.set(drawingBufferSize, width, height));
        }
      },
      setPixelRatio: (value: number) => {
        ValueCell.update(globalUniforms.uPixelRatio, value);
      },
      setOcclusionTest: (f: ((s: Sphere3D) => boolean) | null) => {
        isOccluded = f;
      },

      props: p,
      get stats(): RendererStats {
        return {
          programCount: ctx.stats.resourceCounts.program,
          shaderCount: ctx.stats.resourceCounts.shader,

          attributeCount: ctx.stats.resourceCounts.attribute,
          elementsCount: ctx.stats.resourceCounts.elements,
          framebufferCount: ctx.stats.resourceCounts.framebuffer,
          renderbufferCount: ctx.stats.resourceCounts.renderbuffer,
          textureCount: ctx.stats.resourceCounts.texture,
          vertexArrayCount: ctx.stats.resourceCounts.vertexArray,

          drawCount: stats.drawCount,
          instanceCount: stats.instanceCount,
          instancedDrawCount: stats.instancedDrawCount,
        };
      },
      get light(): Light {
        return light;
      },
      get ambientColor(): Vec3 {
        return globalUniforms.uAmbientColor.ref.value;
      },
      dispose: () => {
        globals.destroy();
        if (emptyDepthTexture !== emptyDepthAlphaTexture) emptyDepthTexture.destroy();
        emptyDepthAlphaTexture.destroy();
        // TODO
      },
    };
  }
}

export { Renderer };
