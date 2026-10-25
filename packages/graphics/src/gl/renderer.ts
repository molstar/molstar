/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { OcclusionPyramid } from './compute/cell-cull.js';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import type { ICamera } from '@molstar/graphics/canvas3d/camera';
import type { Scene } from './scene.js';
import type { GraphicsContext, PipelineRequest, TargetDescriptor } from './context.js';
import { Mat4, Vec3, Vec4, Vec2 } from '@molstar/core/math/linear-algebra';
import { type GraphicsRenderable, type Frame, createFrame, DrawKey, LodView } from './renderable.js';
import {
  DrawFact,
  DrawRow,
  DrawState,
  type DrawTable,
  getCheckFacts,
  getDrawFacts,
  getDrawState,
} from './draw-table.js';
import { Color } from '@molstar/core/util/color';
import { ValueCell, deepEqual } from '@molstar/core/util';
import { GlobalUniformSchema, type GlobalUniformValues, MaxLightCount } from './renderable/schema.js';
import type { DrawList, GraphicsRenderVariant, ProgramItem } from './render-item.js';
import type { ObjectSignature } from './precompile.js';
import { getRenderableProgramSource } from './render-object.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { degToRad } from '@molstar/core/math/misc';
import type { Texture, Textures } from './texture.js';
import { arrayMapUpsert } from '@molstar/core/util/array';
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
  /**
   * Culls on the GPU (merged render objects, `Renderable.prepareCull`) with the planes of the last `update`: outside of
   * passes (compute), before the passes of the frame that draw with it; others draw the CPU culling's lists
   */
  prepareCull: (scene: Scene.Group) => void;
  /** The objects the last `prepareCull` culled on the GPU */
  readonly gpuCulledCount: number;
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
    dpoitTextures: { depth: Texture; frontColor: Texture },
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
  /**
   * The depth pyramid GPU culling tests occlusion against (`prepareCull`, the previous frame's), null for none (e.g.
   * where passes need what is hidden, illumination)
   */
  setOcclusionPyramid: (get: (() => OcclusionPyramid | undefined) | null) => void;

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
  drawLists: PD.Boolean(true, {
    description:
      'Replay the passes from recordings while they draw the same objects (WebGPU render bundles); objects changing every frame are drawn directly.',
  }),
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
   * The render state of a draw of an object with the facts (`DrawFact`): the pass state with the object's cull and
   * front face; `undefined` where the object isn't drawn (direct volumes in the `emissive` and `tracing` variants)
   */
  function getObjectState(
    facts: number,
    variant: GraphicsRenderVariant,
    flag: Flag,
    passState: RenderState,
  ): RenderState | undefined {
    const flipSided = (facts & DrawFact.FlipSided) !== 0;
    const hasReflection = (facts & DrawFact.HasReflection) !== 0;
    let cull: CullMode;
    let frontFace: FrontFace = flipSided ? 'cw' : 'ccw';
    if ((facts & DrawFact.DirectVolume) !== 0) {
      if (variant === 'emissive' || variant === 'tracing') {
        return; // not supported
      }

      // culling done in fragment shader
      cull = 'none';
      frontFace = 'ccw';
    } else if (flag === Flag.BlendedFront) {
      cull = flipSided ? 'front' : 'back';
    } else if (
      flag === Flag.DepthBack ||
      flag === Flag.SolidInteriorMark ||
      (flag === Flag.SolidInteriorFill && hasReflection)
    ) {
      // reflected instances flip the winding, so their exit faces are not back-facing
      cull = 'none';
    } else if (flag === Flag.BlendedBack || flag === Flag.SolidInteriorFill) {
      cull = flipSided ? 'back' : 'front';
    } else if ((facts & DrawFact.SingleSided) !== 0 && !hasReflection) {
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

  /**
   * A check of a pass on an object's facts (`DrawFact`, `getCheckFacts`) and state (`DrawState`): bit tests, the same
   * for an object of a draw table (its row) and any other (`getObjectFacts`)
   */
  type FactTest = (facts: number, state: number) => boolean;

  const OpaqueBlockers =
    DrawFact.AlphaOne | DrawFact.TransparencyFull | DrawFact.DirectVolume | DrawFact.FuzzyPoints | DrawFact.Xray;
  const TransparentReasons =
    DrawFact.AlphaBelowOne |
    DrawFact.TransparencyAny |
    DrawFact.DirectVolume |
    DrawFact.FuzzyPoints |
    DrawFact.TextOrImage |
    DrawFact.Xray;

  /** alpha 1, not fully transparent, not a direct volume, fuzzy points or x-ray shaded; or back faces opaque */
  const checkOpaque: FactTest = (f) =>
    (f & OpaqueBlockers) === DrawFact.AlphaOne || (f & DrawFact.BackfacesOpaque) !== 0;

  /** not alpha 0, and alpha below 1, transparency, a direct volume, fuzzy points, text, an image or x-ray shaded */
  const checkTransparent: FactTest = (f) => (f & DrawFact.AlphaZero) === 0 && (f & TransparentReasons) !== 0;

  const checkPickable: FactTest = (_, s) => (s & DrawState.ColorOnly) === 0;

  const checkMarkingDepth: FactTest = (f) =>
    (f & (DrawFact.AlphaZero | DrawFact.TransparencyFull | DrawFact.MarkerFull)) === 0;

  const checkMarkingMask: FactTest = (f) => (f & DrawFact.MarkerAny) !== 0;

  const checkEmissive: FactTest = (f) => (f & DrawFact.Emissive) !== 0;

  const isDirectVolume: FactTest = (f) => (f & DrawFact.DirectVolume) !== 0;

  /** caps are drawn from back faces, which transparent passes discard when `transparentBackfaces` is 'off' */
  const checkTransparentCap: FactTest = (f, s) => checkTransparent(f, s) && (f & DrawFact.BackfacesOff) === 0;

  const checkAll: FactTest = () => true;

  /** The facts of an object outside a draw table, for the checks (`FactTest`; its state: `getDrawState`) */
  function getObjectFacts(r: GraphicsRenderable) {
    return getCheckFacts(r.values);
  }

  /** whether an object with `state` draws in `variant`: pick passes skip objects that are not pickable */
  function isVariantDrawn(variant: GraphicsRenderVariant, state: number) {
    return variant !== 'pick' || (state & DrawState.Pickable) !== 0;
  }

  /** A draw of an object in a pass */
  type ObjectDraw = { readonly flag: Flag; readonly passState: RenderState };

  /** The draws of a renderer pass, shared by drawing (`renderPass`) and collecting its pipeline requests (`collect`) */
  interface PassDraws {
    readonly variant: GraphicsRenderVariant;
    /** objects the pass is for (by geometry), all if not given */
    readonly select?: FactTest;
    /** objects the pass draws in their current state (opacity, markers, emission, ...) */
    readonly check: FactTest;
    /** the draws of an object with the facts, in order */
    readonly draws: (facts: number) => readonly ObjectDraw[];
    /** solid interior caps of the objects passing `check`, `inline` right after each object, else after all objects */
    readonly caps?: {
      readonly check: FactTest;
      readonly mode: SolidInteriorMode;
      readonly passState: RenderState;
      readonly inline?: boolean;
    };
  }

  function singleDraw(passState: RenderState, flag = Flag.None): (facts: number) => readonly ObjectDraw[] {
    const draws: readonly ObjectDraw[] = [{ flag, passState }];
    return () => draws;
  }

  /** draws objects of `check` with `passState`, caps of the same objects */
  function simplePassDraws(
    variant: GraphicsRenderVariant,
    check: FactTest,
    passState: RenderState,
    caps: { check?: FactTest; mode: SolidInteriorMode },
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
    const check: FactTest = (f, s) => checkOpaque(f, s) && (!occludeWithOpaqueDepth || checkEmissive(f, s));
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
      draws: (f) => {
        if ((f & DrawFact.DoubleSided) === 0) return single;
        return (f & DrawFact.BackfacesOpaque) === 0 ? backAndFront : front;
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
      ...simplePassDraws(
        'emissive',
        (f, s) => checkTransparent(f, s) && checkEmissive(f, s),
        EmissiveTransparentState,
        { check: (f, s) => checkTransparentCap(f, s) && checkEmissive(f, s), mode: 'blended' },
      ),
      select: (f, s) => !isDirectVolume(f, s),
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
    // every global texture from the start: the `dpoit` mode's programs sample the dpoit textures in all their color
    // draws, also those before the first dpoit pass sets them (the opaque pass), which would read whatever their units
    // hold, e.g. a texture attached to the target ('invalid operation' on WebGL)
    const sharedTexturesList: Textures = [
      ['tDepth', emptyDepthTexture],
      ['tDepthAlpha', emptyDepthAlphaTexture],
      ['tDpoitDepth', emptyDepthAlphaTexture],
      ['tDpoitFrontColor', emptyDepthAlphaTexture],
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
    const lodView = LodView();
    /** the culling view of the last update (camera plane, frustum, level scales), see `LodView.version` */
    const lastView = new Float64Array(30);
    let viewVersions = 0;
    /** whether the view was compared since the last update (`updateViewVersion`) */
    let viewChecked = false;
    /** compares the culling view with that of the last cull: a new `lodView.version` where it differs */
    const updateViewVersion = () => {
      let changed = lodView.version === 0;
      const set = (i: number, v: number) => {
        if (lastView[i] !== v) {
          lastView[i] = v;
          changed = true;
        }
      };
      const setPlane = (p: Plane3D, o: number) => {
        set(o, p.normal[0]);
        set(o + 1, p.normal[1]);
        set(o + 2, p.normal[2]);
        set(o + 3, p.constant);
      };
      setPlane(cameraPlane, 0);
      for (let i = 0; i < 6; ++i) setPlane(frustum[i as Frustum3D.PlaneIndex], (i + 1) * 4);
      set(28, lodView.depthScale);
      set(29, lodView.pixelSize);
      if (changed) lodView.version = ++viewVersions;
    };

    let modelScale = 1;
    /** of the object being drawn, with `drawFacts` (`isVisible`, `setRowBoundingSphere`) */
    const boundingSphere = Sphere3D();
    /** the facts of the object being drawn (`DrawFact`): its render states and culling read them */
    let drawFacts = 0;
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

    /** also sets `boundingSphere` to the scaled bounding sphere of `r`, and `drawFacts` */
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

      drawFacts = getDrawFacts(r.values);
      return true;
    };

    /**
     * sets `boundingSphere` to the scaled bounding sphere of a table's row, and `drawFacts` (as `isVisible` does for an
     * object): from the row, without reading the object's values
     */
    const setRowBoundingSphere = (table: DrawTable, slot: number) => {
      const o = slot * DrawRow.Stride;
      const { floats } = table;
      const s = modelScale;
      Vec3.set(
        boundingSphere.center,
        floats[o + DrawRow.CenterX] * s,
        floats[o + DrawRow.CenterY] * s,
        floats[o + DrawRow.CenterZ] * s,
      );
      boundingSphere.radius = floats[o + DrawRow.Radius] * s;
      drawFacts = table.bits[o + DrawRow.Facts];
    };

    /**
     * The planes deciding visibility, the frustum's six and the camera plane (normal, constant), as of the last `update`;
     * with the model scale and the occlusion test, the view what was derived from it (the tables' visibility, the draw
     * lists' snapshots) compares by value: passes of several cameras alternate in a frame (e.g. the camera helper's)
     */
    const viewPlanes = new Float64Array(28);
    const setViewPlane = (o: number, plane: Plane3D) => {
      const { normal, constant } = plane;
      viewPlanes[o] = normal[0];
      viewPlanes[o + 1] = normal[1];
      viewPlanes[o + 2] = normal[2];
      viewPlanes[o + 3] = constant;
    };
    const setViewPlanes = () => {
      for (let i = 0; i < 6; ++i) setViewPlane(i * 4, frustum[i as Frustum3D.PlaneIndex]);
      setViewPlane(24, cameraPlane);
    };
    const isSameView = (planes: Float64Array) => {
      for (let i = 0; i < 28; ++i) if (planes[i] !== viewPlanes[i]) return false;
      return true;
    };

    /**
     * The visibility of a table's rows for a view, by slot (1 visible): state (not disposed, visible), draw count,
     * frustum, lod, occlusion; as `isVisible` without the variant (`isVariantDrawn`). Its key: the table's and the view's
     * versions, the model scale, the occlusion test and with it the frame (Hi-Z changes between frames)
     */
    type RowVisibility = {
      tableVersion: number;
      readonly view: Float64Array;
      modelScale: number;
      isOccluded: ((s: Sphere3D) => boolean) | null;
      frame: Frame | undefined;
      visible: Uint8Array;
    };
    const rowVisibilities = new WeakMap<DrawTable, RowVisibility>();
    const occlusionSphere = Sphere3D();
    const getRowVisibility = (table: DrawTable) => {
      table.refresh();
      let v = rowVisibilities.get(table);
      if (!v) {
        v = {
          tableVersion: -1,
          view: new Float64Array(28).fill(NaN),
          modelScale: NaN,
          isOccluded: null,
          frame: undefined,
          visible: new Uint8Array(0),
        };
        rowVisibilities.set(table, v);
      }
      if (
        v.tableVersion === table.version &&
        isSameView(v.view) &&
        v.modelScale === modelScale &&
        v.isOccluded === isOccluded &&
        (isOccluded === null || v.frame === currentFrame)
      ) {
        return v.visible;
      }
      v.tableVersion = table.version;
      v.view.set(viewPlanes);
      v.modelScale = modelScale;
      v.isOccluded = isOccluded;
      v.frame = currentFrame;
      if (v.visible.length < table.capacity) v.visible = new Uint8Array(table.capacity);
      const { visible } = v;
      const { order, count, floats, bits } = table;
      const s = modelScale;
      const occlude = s === 1 ? isOccluded : null;
      const p = viewPlanes;
      // inlined plane tests: the loop runs over every row per view
      for (let k = 0; k < count; ++k) {
        const slot = order[k];
        const o = slot * DrawRow.Stride;
        visible[slot] = 0;
        if ((bits[o + DrawRow.State] & (DrawState.Disposed | DrawState.Visible)) !== DrawState.Visible) continue;
        if ((bits[o + DrawRow.Facts] & DrawFact.DrawCount) === 0) continue;
        const x = floats[o + DrawRow.CenterX] * s;
        const y = floats[o + DrawRow.CenterY] * s;
        const z = floats[o + DrawRow.CenterZ] * s;
        const r = floats[o + DrawRow.Radius] * s;
        let inside = true;
        for (let i = 0; i < 24; i += 4) {
          if (p[i] * x + p[i + 1] * y + p[i + 2] * z + p[i + 3] < -r) {
            inside = false;
            break;
          }
        }
        if (!inside) continue;
        const minDistance = floats[o + DrawRow.LodMin];
        const maxDistance = floats[o + DrawRow.LodMax];
        if (minDistance !== 0 || maxDistance !== 0) {
          const d = p[24] * x + p[25] * y + p[26] * z + p[27];
          if (d + r < minDistance * s) continue;
          if (d - r > maxDistance * s) continue;
        }
        if (occlude !== null) {
          Vec3.set(occlusionSphere.center, x, y, z);
          occlusionSphere.radius = r;
          if (occlude(occlusionSphere)) continue;
        }
        visible[slot] = 1;
      }
      return visible;
    };

    /** the instances and levels of detail of `r` to draw; assumes `isVisible` returned true (`boundingSphere`) */
    const cullObject = (r: GraphicsRenderable) => {
      const unscaled = modelScale === 1;
      const lodLevels = (drawFacts & DrawFact.LodLevels) !== 0;
      if (unscaled) {
        const hasInstanceGrid = (drawFacts & DrawFact.InstanceGrid) !== 0;
        const hasMultipleInstances = (drawFacts & DrawFact.MultipleInstances) !== 0;
        if (hasInstanceGrid && (hasMultipleInstances || lodLevels)) {
          if (!viewChecked) {
            updateViewVersion();
            viewChecked = true;
          }
          r.cull(cameraPlane, frustum, isOccluded, ctx.stats, currentFrame, lodView);
        } else {
          r.uncull(lodView);
        }
      } else {
        if (lodLevels) {
          const { center, radius } = boundingSphere;
          const d = Plane3D.distanceToPoint(cameraPlane, center);
          r.cullSimple(d, radius, modelScale, lodView);
        } else {
          r.uncull(lodView);
        }
      }
    };

    /** assumes `isVisible` returned true, which also set `boundingSphere` */
    const drawObject = (r: GraphicsRenderable, variant: GraphicsRenderVariant, flag: Flag, passState: RenderState) => {
      cullObject(r);
      const state = getObjectState(drawFacts, variant, flag, passState);
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
    /** of an object with the facts (`DrawFact.SolidInterior`) */
    const hasSolidInteriorCap = (facts: number) => {
      return solidInteriorCapSupported && (facts & DrawFact.SolidInterior) !== 0;
    };

    const setSolidInteriorPass = (pass: number, clipIndex: number) => {
      globals.setSolidInteriorCap(pass, clipIndex);
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

    /**
     * The caps `r` draws in this frame into `out`: -1 for the near plane, the indices of the clip objects that cut it
     * (camera dependent); assumes `isVisible` returned true, which also set `boundingSphere`
     */
    const getSolidInteriorCaps = (r: GraphicsRenderable, mode: SolidInteriorMode, out: number[]) => {
      out.length = 0;
      const back = mode === 'back';
      const near = globalUniforms.uNear.ref.value * 1.0001;
      if (
        !back &&
        Math.abs(Plane3D.distanceToPoint(cameraPlane, boundingSphere.center) - near) <= boundingSphere.radius
      ) {
        out.push(-1);
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
          out.push(i);
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
        renderSolidInteriorMark(r, variant, caps[j], passState);
        renderSolidInteriorFill(r, variant, mode, caps[j], passState);
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
      if (!solidInteriorCapSupported) return;
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

    /** the camera scale of the last `update`: culled on the GPU only unscaled, as `cull` tests (`cullObject`) */
    let updateScale = 1;
    let getOcclusionPyramid: (() => OcclusionPyramid | undefined) | null = null;
    let gpuCulledCount = 0;
    const prepareCull = (scene: Scene.Group) => {
      gpuCulledCount = 0;
      if (updateScale !== 1) return;
      const pyramid = getOcclusionPyramid?.();
      const { renderables } = scene;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        if (renderables[i].prepareCull?.(cameraPlane, frustum, currentFrame, lodView, pyramid)) gpuCulledCount += 1;
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
      // levels of detail by size on screen: in CSS pixels, of the camera's viewport (not a pass's target)
      LodView.set(
        lodView,
        camera.projection,
        camera.viewport.height / ctx.pixelRatio,
        camera.state.mode === 'orthographic',
      );

      Plane3D.copy(cameraPlane, frustum[Frustum3D.PlaneIndex.Near]);
      cameraPlane.constant -= Plane3D.distanceToPoint(cameraPlane, cameraPosition);
      // compared at the first cull after this update: updates of other cameras between (picking) don't cull
      viewChecked = false;
      setViewPlanes();
      ValueCell.update(
        globalUniforms.uCameraPlane,
        Plane3D.toArray(cameraPlane, globalUniforms.uCameraPlane.ref.value, 0),
      );
      updateScale = camera.scale;

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

    /** draws the objects of `group` as `pass` describes: from its table's rows where it has one */
    const renderPass = (pass: PassDraws, group: Scene.Group) => {
      if (group.table) renderPassRows(pass, group.table);
      else renderPassObjects(pass, group);
    };

    /**
     * A part of a pass over a table's order: the objects of `[begin, end)` with their inline caps (`objects`), the caps
     * drawn after all objects of the range (`caps`); a chunk of a draw list's recordings
     */
    type PassPart = { readonly begin: number; readonly end: number; readonly objects: boolean; readonly caps: boolean };
    const WholePass: PassPart = { begin: 0, end: Infinity, objects: true, caps: true };

    /** `renderPass` from a table's rows (of `part`): visibility per view, the checks on bits */
    const renderPassRows = (pass: PassDraws, table: DrawTable, part: PassPart = WholePass) => {
      const { variant, select, check, draws, caps } = pass;
      const visible = getRowVisibility(table);
      const inlineCaps = !!caps?.inline && canRenderSolidInteriorCaps();
      const { order, bits, objects } = table;
      const begin = part.begin;
      const end = Math.min(part.end, table.count);
      for (let k = part.objects ? begin : end; k < end; ++k) {
        const slot = order[k];
        if (!visible[slot]) continue;
        const o = slot * DrawRow.Stride;
        const facts = bits[o + DrawRow.Facts];
        const state = bits[o + DrawRow.State];
        if (!isVariantDrawn(variant, state) || (select && !select(facts, state)) || !check(facts, state)) continue;
        const r = objects[slot]!;
        setRowBoundingSphere(table, slot);
        const objectDraws = draws(facts);
        for (let j = 0, jl = objectDraws.length; j < jl; ++j) {
          drawObject(r, variant, objectDraws[j].flag, objectDraws[j].passState);
        }
        if (inlineCaps && hasSolidInteriorCap(facts) && caps!.check(facts, state)) {
          drawSolidInteriorCaps(r, variant, caps!.mode, caps!.passState);
        }
      }
      if (part.caps && caps && !caps.inline && canRenderSolidInteriorCaps()) {
        for (let k = begin; k < end; ++k) {
          const slot = order[k];
          if (!visible[slot]) continue;
          const o = slot * DrawRow.Stride;
          const facts = bits[o + DrawRow.Facts];
          const state = bits[o + DrawRow.State];
          if (!isVariantDrawn(variant, state) || !hasSolidInteriorCap(facts) || !caps.check(facts, state)) continue;
          const r = objects[slot]!;
          setRowBoundingSphere(table, slot);
          drawSolidInteriorCaps(r, variant, caps.mode, caps.passState);
        }
      }
    };

    /** `renderPass` for groups without a table (not a scene's) */
    const renderPassObjects = (pass: PassDraws, group: Scene.Group) => {
      const { variant, select, check, draws, caps } = pass;
      const inlineCaps = !!caps?.inline && canRenderSolidInteriorCaps();
      const { renderables } = group;
      for (let i = 0, il = renderables.length; i < il; ++i) {
        const r = renderables[i];
        const facts = getObjectFacts(r);
        const state = getDrawState(r.state);
        if ((select && !select(facts, state)) || !check(facts, state)) continue;
        const objectDraws = draws(facts);
        for (let j = 0, jl = objectDraws.length; j < jl; ++j) {
          renderObject(r, variant, objectDraws[j].flag, objectDraws[j].passState);
        }
        if (inlineCaps && hasSolidInteriorCap(facts) && caps!.check(facts, state) && isVisible(r, variant)) {
          drawSolidInteriorCaps(r, variant, caps!.mode, caps!.passState);
        }
      }
      if (caps && !caps.inline) renderSolidInteriorCaps(renderables, caps.check, variant, caps.mode, caps.passState);
    };

    /**
     * What a recording drew, by the rows of its table (`walkPass`): per draw its slot, object, render state and the
     * versions of its program, values (`DrawTable.valuesVersions`) and draws (`Renderable.drawVersion`), and what it
     * binds (`digests`, for data changes); the caps as a flat list (per object drawing caps: slot, program and values
     * versions, count, caps); and what the walk read besides the rows, which unchanged with the table means the same
     * draws (`isSnapshotCurrent`). Written in place, so that comparing a frame allocates nothing
     */
    type Snapshot = {
      valid: boolean;
      tableVersion: number;
      valuesEpoch: number;
      readonly view: Float64Array;
      modelScale: number;
      near: number;
      caps: boolean;
      /** `GraphicsContext.drawResourcesVersion` */
      drawResources: number;
      count: number;
      slots: Uint32Array;
      programVersions: Uint32Array;
      valuesVersions: Uint32Array;
      drawVersions: Uint32Array;
      readonly objects: GraphicsRenderable[];
      readonly states: RenderState[];
      readonly digests: (string | undefined)[];
      readonly capKey: DrawKey;
    };
    const createSnapshot = (): Snapshot => ({
      valid: false,
      tableVersion: -1,
      valuesEpoch: -1,
      view: new Float64Array(28).fill(NaN),
      modelScale: NaN,
      near: NaN,
      caps: false,
      drawResources: -1,
      count: 0,
      slots: new Uint32Array(0),
      programVersions: new Uint32Array(0),
      valuesVersions: new Uint32Array(0),
      drawVersions: new Uint32Array(0),
      objects: [],
      states: [],
      digests: [],
      capKey: DrawKey(),
    });
    const ensureSnapshotDraws = (snap: Snapshot, n: number) => {
      if (snap.slots.length > n) return;
      const grow = (a: Uint32Array) => {
        const b = new Uint32Array(Math.max(16, n * 2));
        b.set(a);
        return b;
      };
      snap.slots = grow(snap.slots);
      snap.programVersions = grow(snap.programVersions);
      snap.valuesVersions = grow(snap.valuesVersions);
      snap.drawVersions = grow(snap.drawVersions);
    };
    const isCapsDrawn = () => solidInteriorCapSupported && canRenderSolidInteriorCaps();
    const setSnapshotInputs = (snap: Snapshot, table: DrawTable) => {
      snap.tableVersion = table.version;
      snap.valuesEpoch = table.valuesEpoch;
      snap.view.set(viewPlanes);
      snap.modelScale = modelScale;
      snap.near = globalUniforms.uNear.ref.value;
      snap.caps = isCapsDrawn();
      snap.drawResources = ctx.drawResourcesVersion ?? 0;
    };
    /**
     * Whether nothing the walk read changed (the rows and order, values, programs, the view and model scale, the near
     * plane, whether caps are drawn): the snapshot holds without walking. Not with an occlusion test (Hi-Z changes
     * between frames)
     */
    const isSnapshotCurrent = (snap: Snapshot, table: DrawTable) => {
      if (!snap.valid || isOccluded !== null) return false;
      // volatile cells compared
      table.refresh();
      return snap.tableVersion === table.version && snap.valuesEpoch === table.valuesEpoch && isSameInputs(snap);
    };
    /** Whether what the walk read besides the table is unchanged (see `isSnapshotCurrent`) */
    const isSameInputs = (snap: Snapshot) =>
      isSameView(snap.view) &&
      snap.modelScale === modelScale &&
      snap.near === globalUniforms.uNear.ref.value &&
      snap.caps === isCapsDrawn() &&
      snap.drawResources === (ctx.drawResourcesVersion ?? 0);

    const walkCapList: number[] = [];
    /**
     * The caps `r` (in `slot`) draws, at `c` of the snapshot's cap list: written, or compared (-1 where they differ);
     * nothing where it draws none (one starting or stopping to draw caps changes the list all the same). Assumes `r`
     * is visible, `boundingSphere` set; returns the position after its entry
     */
    const walkCaps = (
      r: GraphicsRenderable,
      slot: number,
      valuesVersion: number,
      mode: SolidInteriorMode,
      key: DrawKey,
      c: number,
      write: boolean,
    ) => {
      const caps = getSolidInteriorCaps(r, mode, walkCapList);
      if (caps.length === 0) return c;
      if (write) {
        DrawKey.add(key, slot);
        DrawKey.add(key, r.programVersion);
        DrawKey.add(key, valuesVersion);
        DrawKey.add(key, caps.length);
        for (let k = 0, kl = caps.length; k < kl; ++k) DrawKey.add(key, caps[k]);
        return key.length;
      }
      const v = key.values;
      if (c + 4 + caps.length > key.length) return -1;
      if (v[c] !== slot || v[c + 1] !== r.programVersion || v[c + 2] !== valuesVersion || v[c + 3] !== caps.length) {
        return -1;
      }
      for (let k = 0, kl = caps.length; k < kl; ++k) {
        if (v[c + 4 + k] !== caps[k]) return -1;
      }
      return c + 4 + caps.length;
    };

    /**
     * The draws of `pass` as `renderPassRows` draws them (culled the same way) with the caps they draw (the near plane
     * and clip objects cutting objects, camera dependent), written into `snap` or compared with it: 'same'; 'data'
     * where only values of draws changed (their data may have, see `flushDataChanges`); 'different' at the first draw,
     * state, program, range or cap that differs
     */
    const walkPass = (
      pass: PassDraws,
      table: DrawTable,
      snap: Snapshot,
      write: boolean,
      part: PassPart = WholePass,
    ): 'same' | 'data' | 'different' => {
      const { variant, select, check, draws, caps } = pass;
      const withCaps = !!caps && isCapsDrawn();
      const visible = getRowVisibility(table);
      const { order, bits, objects, valuesVersions } = table;
      const begin = part.begin;
      const end = Math.min(part.end, table.count);
      if (write) snap.capKey.length = 0;
      let n = 0;
      let c = 0;
      let data = false;
      for (let k = part.objects ? begin : end; k < end; ++k) {
        const slot = order[k];
        if (!visible[slot]) continue;
        const o = slot * DrawRow.Stride;
        const facts = bits[o + DrawRow.Facts];
        const state = bits[o + DrawRow.State];
        if (!isVariantDrawn(variant, state) || (select && !select(facts, state)) || !check(facts, state)) continue;
        const r = objects[slot]!;
        setRowBoundingSphere(table, slot);
        const objectDraws = draws(facts);
        if (objectDraws.length > 0) {
          cullObject(r);
          const { programVersion, drawVersion } = r;
          const valuesVersion = valuesVersions[slot];
          for (let j = 0, jl = objectDraws.length; j < jl; ++j) {
            const renderState = getObjectState(facts, variant, objectDraws[j].flag, objectDraws[j].passState);
            if (!renderState) continue;
            if (write) {
              ensureSnapshotDraws(snap, n);
              snap.slots[n] = slot;
              snap.objects[n] = r;
              snap.states[n] = renderState;
              snap.programVersions[n] = programVersion;
              snap.valuesVersions[n] = valuesVersion;
              snap.drawVersions[n] = drawVersion;
              snap.digests[n] = undefined;
            } else {
              if (
                n >= snap.count ||
                snap.slots[n] !== slot ||
                snap.objects[n] !== r ||
                snap.states[n] !== renderState ||
                snap.programVersions[n] !== programVersion ||
                snap.drawVersions[n] !== drawVersion
              ) {
                return 'different';
              }
              if (snap.valuesVersions[n] !== valuesVersion) data = true;
            }
            n += 1;
          }
        }
        if (withCaps && caps!.inline && hasSolidInteriorCap(facts) && caps!.check(facts, state)) {
          c = walkCaps(r, slot, valuesVersions[slot], caps!.mode, snap.capKey, c, write);
          if (c < 0) return 'different';
        }
      }
      if (part.caps && withCaps && !caps!.inline) {
        for (let k = begin; k < end; ++k) {
          const slot = order[k];
          if (!visible[slot]) continue;
          const o = slot * DrawRow.Stride;
          const facts = bits[o + DrawRow.Facts];
          const state = bits[o + DrawRow.State];
          if (!isVariantDrawn(variant, state) || !hasSolidInteriorCap(facts) || !caps!.check(facts, state)) continue;
          const r = objects[slot]!;
          setRowBoundingSphere(table, slot);
          c = walkCaps(r, slot, valuesVersions[slot], caps!.mode, snap.capKey, c, write);
          if (c < 0) return 'different';
        }
      }
      if (write) {
        snap.count = n;
        // shorter only where draws were dropped: no references kept to their objects
        if (snap.objects.length > n) {
          snap.objects.length = n;
          snap.states.length = n;
          snap.digests.length = n;
        }
        return 'same';
      }
      if (n !== snap.count || c !== snap.capKey.length) return 'different';
      return data ? 'data' : 'same';
    };

    /** What the `i`th draw of `snap` binds, its data written as by drawing it (`GraphicsContext.captureDraws`) */
    const captureDraw = (snap: Snapshot, i: number, variant: GraphicsRenderVariant) => {
      return ctx.captureDraws!(() => {
        ctx.setRenderState(snap.states[i]);
        snap.objects[i].render(variant, globals);
      });
    };

    /**
     * The draws of `snap` whose values changed written as by drawing them, without encoding them: true where each binds
     * what it bound when recorded (only data changed, e.g. markers on hover), so the recording replays with it; false
     * where one binds otherwise (a texture or buffer replaced, a draw count changed): record again
     */
    const flushDataChanges = (snap: Snapshot, table: DrawTable, variant: GraphicsRenderVariant) => {
      if (!ctx.captureDraws) return false;
      const { valuesVersions } = table;
      for (let i = 0, il = snap.count; i < il; ++i) {
        const version = valuesVersions[snap.slots[i]];
        if (snap.valuesVersions[i] === version) continue;
        const digest = snap.digests[i];
        if (digest === undefined || captureDraw(snap, i, variant) !== digest) return false;
        snap.valuesVersions[i] = version;
      }
      return true;
    };

    /**
     * The draw lists of a group's pass, one per time it is drawn between two submissions (stereo eyes, sub-renders, the
     * pick types, marking samples): the n-th draw of the pass after a submission binds the n-th global slot each time
     * (`GraphicsContext.submitCount`), whether a frame or a marking redraw without a new frame; per pass, so that
     * submissions drawing other passes (a pick, a marking redraw) don't shift the lists of a pass
     */
    type PassDrawLists = {
      submit: number;
      next: number;
      lists: PassDrawList[];
    };
    /**
     * A recording of a chunk (`part`) of a pass: its draw list and what it drew (`snapshot`, `valid` once recorded);
     * for a chunk of objects, the rows its walk read (`rows`: per position of its range the slot, its row and values
     * versions; `invalidations` of the table): unchanged, it holds while other objects change (`isChunkCurrent`)
     */
    type Chunk = {
      readonly list: DrawList;
      readonly snapshot: Snapshot;
      readonly part: PassPart;
      rows: Uint32Array;
      rowCount: number;
      invalidations: number;
      /** the visibility of its rows per position (`getRowVisibility`) when drawn */
      readonly visible: Uint8Array;
      /** whether a visible row's draws follow the view besides its visibility (`isViewDependent`) */
      viewDependent: boolean;
      /** the submission its inputs last changed in (`drawChunk`: drawn directly while changing) */
      changedAt: number;
      /** how it was drawn last (`drawChunk`) */
      outcome: ChunkOutcome;
    };
    type ChunkOutcome = 'unplanned' | 'replayed' | 'flushed' | 'empty' | 'recorded' | 'direct';
    /**
     * The recordings of a pass in chunks of `ChunkSize` objects of the table's order, then one of the caps drawn after
     * all objects: a change of what the pass draws (an object marked, hidden, entering the view) records the chunks it
     * touches again, the others replay (a hover sweeping across many objects, a zoomed view rotating)
     */
    type PassDrawList = {
      readonly chunks: Chunk[];
      caps: Chunk | undefined;
      /** what the walks read when every chunk was last drawn: unchanged, all replay without a check of their own */
      readonly inputs: Snapshot;
    };
    const ChunkSize = 64;
    /**
     * Submissions a chunk's inputs stay unchanged before it is recorded again after changing in consecutive ones: drawn
     * directly meanwhile, its data or view changing every frame (a hover, a recolor, a rotation) would else be walked,
     * flushed or recorded each frame, which cost more than drawing it (merged hover 19-21%, unmerged recolor 17-21%)
     */
    const StableSubmissions = 3;
    /**
     * Rows of a chunk drawn directly as cheaply as walked and flushed (a merged object is one row): more rows changing
     * every frame (a hover sweeping across separate objects, a rotation) are flushed and replayed, cheaper than drawn
     */
    const FewChunkRows = 4;
    /** The caps drawn after all objects, a chunk of their own */
    const CapsPart: PassPart = { begin: 0, end: Infinity, objects: false, caps: true };
    let groupDrawLists = new WeakMap<Scene.Group, Map<PassDraws, PassDrawLists>>();
    const resetDrawLists = () => {
      groupDrawLists = new WeakMap();
    };
    const getDrawList = (pass: PassDraws, group: Scene.Group) => {
      if (!ctx.createDrawList) return;
      let passLists = groupDrawLists.get(group);
      if (!passLists) groupDrawLists.set(group, (passLists = new Map()));
      let lists = passLists.get(pass);
      if (!lists) passLists.set(pass, (lists = { submit: -1, next: 0, lists: [] }));
      const submit = ctx.submitCount ?? 0;
      if (lists.submit !== submit) {
        lists.submit = submit;
        lists.next = 0;
      }
      const i = lists.next++;
      return (lists.lists[i] ??= { chunks: [], caps: undefined, inputs: createSnapshot() });
    };
    const createChunk = (part: PassPart): Chunk => ({
      list: ctx.createDrawList!(),
      snapshot: createSnapshot(),
      part,
      rows: new Uint32Array(part.objects ? 3 * ChunkSize : 0),
      rowCount: -1,
      invalidations: -1,
      visible: new Uint8Array(part.objects ? ChunkSize : 0),
      viewDependent: true,
      changedAt: -Infinity,
      outcome: 'empty',
    });
    /**
     * Whether the draws of a row (`facts`) follow the view besides whether it is visible: instances culled
     * (`cullObject`), levels of detail, solid interior caps (cut by the near plane)
     */
    const isViewDependent = (facts: number) =>
      (facts & (DrawFact.LodLevels | DrawFact.SolidInterior)) !== 0 ||
      (facts & (DrawFact.InstanceGrid | DrawFact.MultipleInstances)) ===
        (DrawFact.InstanceGrid | DrawFact.MultipleInstances);
    /**
     * Whether nothing the walk of `chunk` read changed: for a chunk of objects its rows (slots in order, row and values
     * versions) and the table's invalidations, so that a change of other objects leaves it current; the caps chunk
     * reads all objects (`isSnapshotCurrent`). Assumes the table refreshed
     */
    const isChunkCurrent = (chunk: Chunk, table: DrawTable) =>
      chunk.part.objects
        ? chunk.snapshot.valid && isChunkInputsCurrent(chunk, table)
        : isSnapshotCurrent(chunk.snapshot, table);
    /** `isChunkCurrent` of a chunk of objects but its recording: whether what it would draw is as when last drawn */
    const isChunkInputsCurrent = (chunk: Chunk, table: DrawTable) => {
      const { snapshot, part, rows } = chunk;
      if (chunk.rowCount < 0 || isOccluded !== null) return false;
      const end = Math.min(part.end, table.count);
      if (snapshot.tableVersion !== table.version || snapshot.valuesEpoch !== table.valuesEpoch) {
        if (chunk.invalidations !== table.invalidations || chunk.rowCount !== end - part.begin) return false;
        const { order, rowVersions, valuesVersions } = table;
        for (let k = part.begin, i = 0; k < end; ++k, i += 3) {
          const slot = order[k];
          if (rows[i] !== slot || rows[i + 1] !== rowVersions[slot] || rows[i + 2] !== valuesVersions[slot]) {
            return false;
          }
        }
      }
      if (isSameInputs(snapshot)) return true;
      // a recording binding a replaced resource is not current, whatever the view did
      if (snapshot.drawResources !== (ctx.drawResourcesVersion ?? 0)) return false;
      // the view changed (a camera rotating): the same draws where its rows are as visible as they were and none of them
      // draws by the view otherwise
      if (chunk.viewDependent || snapshot.modelScale !== modelScale || snapshot.caps !== isCapsDrawn()) return false;
      const visible = getRowVisibility(table);
      const { order } = table;
      for (let k = part.begin, i = 0; k < end; ++k, ++i) {
        if (visible[order[k]] !== chunk.visible[i]) return false;
      }
      return true;
    };
    /** `chunk` drawn as the table is now (`isChunkCurrent`) */
    const setChunkInputs = (chunk: Chunk, table: DrawTable) => {
      const { snapshot, part, rows } = chunk;
      setSnapshotInputs(snapshot, table);
      if (!part.objects) return;
      const end = Math.min(part.end, table.count);
      const { order, rowVersions, valuesVersions } = table;
      for (let k = part.begin, i = 0; k < end; ++k, i += 3) {
        const slot = order[k];
        rows[i] = slot;
        rows[i + 1] = rowVersions[slot];
        rows[i + 2] = valuesVersions[slot];
      }
      chunk.rowCount = Math.max(0, end - part.begin);
      chunk.invalidations = table.invalidations;
      const visible = getRowVisibility(table);
      const { bits } = table;
      let viewDependent = false;
      for (let k = part.begin, i = 0; k < end; ++k, ++i) {
        const slot = order[k];
        chunk.visible[i] = visible[slot];
        if (visible[slot] && isViewDependent(bits[slot * DrawRow.Stride + DrawRow.Facts])) viewDependent = true;
      }
      chunk.viewDependent = viewDependent;
    };
    const isEmptySnapshot = (snapshot: Snapshot) => snapshot.count === 0 && snapshot.capKey.length === 0;

    /**
     * Draws `part` of `pass` from `chunk`: 'unplanned' replayed without walking (nothing the walk read changed),
     * 'replayed' or 'flushed' (its data written first) after walking, 'empty' drawing nothing, 'recorded', 'direct'
     * drawn without its recording while its inputs change in consecutive submissions (`StableSubmissions`)
     */
    const drawChunk = (pass: PassDraws, table: DrawTable, chunk: Chunk): ChunkOutcome => {
      const outcome = drawChunkOnce(pass, table, chunk);
      chunk.outcome = outcome;
      return outcome;
    };
    /**
     * Whether `chunk`, changed again in the submission after its last change, is drawn directly: where its recording
     * was made again (draws changing every frame, a recolor), or it was drawn directly, or it has few rows
     */
    const isDrawnDirectly = (chunk: Chunk) =>
      chunk.outcome === 'recorded' || chunk.outcome === 'direct' || chunk.rowCount <= FewChunkRows;
    const drawChunkOnce = (pass: PassDraws, table: DrawTable, chunk: Chunk): ChunkOutcome => {
      const { snapshot, list, part } = chunk;
      if (part.objects) {
        const submit = ctx.submitCount ?? 0;
        const current = isChunkInputsCurrent(chunk, table);
        if (current && snapshot.valid) {
          if (isEmptySnapshot(snapshot)) return 'empty';
          if (list.replay(globals)) return 'unplanned';
        } else if (
          current
            ? submit - chunk.changedAt < StableSubmissions && chunk.outcome === 'direct'
            : submit - chunk.changedAt <= 1 && isDrawnDirectly(chunk)
        ) {
          // changing every frame, or not stable for long enough since: drawn as without draw lists
          if (!current) chunk.changedAt = submit;
          renderPassRows(pass, table, part);
          // its recording no longer of the draws, released once
          if (snapshot.valid) list.reset();
          snapshot.valid = false;
          setChunkInputs(chunk, table);
          return 'direct';
        } else if (!current) {
          chunk.changedAt = submit;
        }
      } else if (isChunkCurrent(chunk, table)) {
        if (isEmptySnapshot(snapshot)) return 'empty';
        if (list.replay(globals)) return 'unplanned';
      }
      // a recording binding a replaced resource (`GraphicsContext.drawResourcesVersion`) does not replay, whatever the
      // walk finds
      const bound = snapshot.drawResources === (ctx.drawResourcesVersion ?? 0);
      const match = snapshot.valid && bound ? walkPass(pass, table, snapshot, false, part) : 'different';
      if (match === 'same' || (match === 'data' && flushDataChanges(snapshot, table, pass.variant))) {
        if (isEmptySnapshot(snapshot)) {
          setChunkInputs(chunk, table);
          return 'empty';
        }
        if (list.replay(globals)) {
          setChunkInputs(chunk, table);
          return match === 'data' ? 'flushed' : 'replayed';
        }
      } else {
        // the draws changed: so are the recordings for other targets
        list.reset();
      }
      // as drawn directly, the caps (and their clears) included
      list.record(globals, () => renderPassRows(pass, table, part));
      // drawing may change values (e.g. what culling writes): the snapshot after it is the recorded one
      walkPass(pass, table, snapshot, true, part);
      if (ctx.captureDraws) {
        for (let i = 0, il = snapshot.count; i < il; ++i) snapshot.digests[i] = captureDraw(snapshot, i, pass.variant);
      }
      setChunkInputs(chunk, table);
      snapshot.valid = true;
      // a chunk drawing nothing keeps no recording
      if (isEmptySnapshot(snapshot)) list.reset();
      return 'recorded';
    };

    /**
     * The chunks of `drawList` replayed as recorded (those drawing nothing skipped), the caps chunk last (index
     * `chunks.length`): the index of the first that can't replay (no recording for this target, globals bound
     * otherwise), from which the chunks are drawn one by one; -1 where all replayed
     */
    const replayChunks = (drawList: PassDrawList) => {
      const { chunks, caps } = drawList;
      for (let c = 0, il = chunks.length; c < il; ++c) {
        const { snapshot, list } = chunks[c];
        // drawn directly last (`drawChunk`): no recording, its snapshot not of the draws
        if (!snapshot.valid) return c;
        if (!isEmptySnapshot(snapshot) && !list.replay(globals)) return c;
      }
      if (caps && !isEmptySnapshot(caps.snapshot) && !caps.list.replay(globals)) return chunks.length;
      return -1;
    };

    /**
     * Draws `pass` from a draw list (`drawLists`) chunk by chunk (`PassDrawList`): replays each chunk's recording while
     * it draws the recorded draws with unchanged values, caps and globals, records it otherwise (as `renderPass` draws
     * it); false without draw lists
     */
    const renderPassFromDrawList = (pass: PassDraws, group: Scene.Group) => {
      const { table } = group;
      if (!table) return false;
      const drawList = getDrawList(pass, group);
      if (!drawList) return false;
      const { chunks, inputs } = drawList;
      // volatile cells compared, once for all chunks
      table.refresh();
      // nothing any chunk's walk read changed: each replays its recording, from the first that can't one by one
      let first = 0;
      if (isSnapshotCurrent(inputs, table)) {
        first = replayChunks(drawList);
        if (first < 0) {
          ctx.stats.drawLists.replayed += 1;
          ctx.stats.drawLists.unplanned += 1;
          return true;
        }
      }
      inputs.valid = false;
      const chunkCount = Math.ceil(table.count / ChunkSize);
      while (chunks.length < chunkCount) {
        const begin = chunks.length * ChunkSize;
        chunks.push(createChunk({ begin, end: begin + ChunkSize, objects: true, caps: false }));
      }
      // objects removed: their chunks' recordings go
      if (chunks.length > chunkCount) chunks.length = chunkCount;
      let recorded = 0;
      let direct = 0;
      let unplanned = true;
      let flushed = false;
      const count = (result: ReturnType<typeof drawChunk>) => {
        if (result === 'recorded') recorded += 1;
        if (result === 'direct') direct += 1;
        if (result !== 'unplanned' && result !== 'empty') unplanned = false;
        if (result === 'flushed') flushed = true;
      };
      for (let c = first; c < chunkCount; ++c) count(drawChunk(pass, table, chunks[c]));
      if (pass.caps && !pass.caps.inline) {
        drawList.caps ??= createChunk(CapsPart);
        count(drawChunk(pass, table, drawList.caps));
      }
      setSnapshotInputs(inputs, table);
      inputs.valid = true;
      const { drawLists } = ctx.stats;
      drawLists.chunksDirect += direct;
      if (recorded > 0) {
        drawLists.recorded += 1;
        drawLists.chunksRecorded += recorded;
      } else {
        drawLists.replayed += 1;
        if (unplanned) drawLists.unplanned += 1;
        if (flushed) drawLists.flushed += 1;
      }
      return true;
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
          const values: { [k: string]: ValueCell<unknown> } = {
            alpha: ValueCell.create(alpha),
            uAlpha: ValueCell.create(alpha),
          };
          for (const k of Object.keys(fixed)) if (!schema || schema[k]) values[k] = ValueCell.create(fixed[k]);
          return {
            values,
            renderItem: item,
            state: { disposed: false, colorOnly: false },
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
        const state = getObjectState(getObjectFacts(r), variant, flag, passState);
        if (state) requests.push({ item: r.renderItem, variant, state, target });
      };
      for (const r of renderables) {
        const facts = getObjectFacts(r);
        const state = getDrawState(r.state);
        if (r.state.disposed || (select && !select(facts, state)) || (!all && !check(facts, state))) continue;
        for (const { flag, passState } of draws(facts)) request(r, flag, passState);
        if (caps && capsDrawn && hasSolidInteriorCap(facts) && (all || caps.check(facts, state))) {
          request(r, Flag.SolidInteriorMark, getSolidInteriorMarkState(caps.passState));
          const { resolve, fill } = getSolidInteriorFillStates(caps.passState, caps.mode);
          if (resolve) request(r, Flag.SolidInteriorFill, resolve);
          request(r, Flag.SolidInteriorFill, fill);
        }
      }
      return requests;
    };

    /** `pass` replayed or recorded with the draw lists (`drawLists`) where it can be, drawn directly otherwise */
    const drawPass = (pass: PassDraws, group: Scene.Group) => {
      if (p.drawLists && renderPassFromDrawList(pass, group)) return;
      renderPass(pass, group);
    };

    const renderPick = (group: Scene.Group, camera: ICamera, variant: 'pick' | 'depth', pickType: PickType) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uPickType, pickType);
      drawPass(variant === 'pick' ? PassDraws.pick : PassDraws.pickDepth, group);
    };

    const renderDepth = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      drawPass(PassDraws.depth, group);
    };

    const renderDepthOpaque = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      drawPass(PassDraws.depthOpaque, group);
    };

    const renderDepthOpaqueBack = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Opaque, false);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, true);
      renderPass(PassDraws.depthOpaqueBack, group);
      ValueCell.updateIfChanged(globalUniforms.uDepthBack, false);
    };

    const renderDepthTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.depthTransparent, group);
    };

    const renderMarkingDepth = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.All, false);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Depth);
      drawPass(PassDraws.markingDepth, group);
    };

    const renderMarkingMask = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, null, depthTexture, Mask.All, !!depthTexture);
      ValueCell.updateIfChanged(globalUniforms.uMarkingType, MarkingType.Mask);
      drawPass(PassDraws.markingMask, group);
    };

    const renderEmissiveOpaque = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture,
      occludeWithOpaqueDepth: boolean,
    ) => {
      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Opaque, false);
      drawPass(occludeWithOpaqueDepth ? PassDraws.emissiveOpaqueOccluded : PassDraws.emissiveOpaque, group);
    };

    const renderEmissiveTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      // tDepthAlpha = front transparent depth+alpha; emissive shader dims emitters behind it
      updateInternal(group, camera, null, depthTexture, Mask.Transparent, false);
      drawPass(PassDraws.emissiveTransparent, group);
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
      drawPass(PassDraws.opaque, group);
    };

    const renderBlendedTransparent = (group: Scene.Group, camera: ICamera) => {
      updateInternal(group, camera, null, null, Mask.Transparent, false);
      drawPass(transparentBackground ? PassDraws.blendedTransparentBackground : PassDraws.blendedTransparent, group);
    };

    const renderVolume = (group: Scene.Group, camera: ICamera, depthTexture: Texture) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.volume, group);
    };

    const renderWboitTransparent = (group: Scene.Group, camera: ICamera, depthTexture: Texture | null) => {
      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.wboitTransparent, group);
    };

    const renderDpoitTransparent = (
      group: Scene.Group,
      camera: ICamera,
      depthTexture: Texture,
      dpoitTextures: { depth: Texture; frontColor: Texture },
    ) => {
      arrayMapUpsert(sharedTexturesList, 'tDpoitDepth', dpoitTextures.depth);
      arrayMapUpsert(sharedTexturesList, 'tDpoitFrontColor', dpoitTextures.frontColor);

      updateInternal(group, camera, depthTexture, null, Mask.Transparent, false);
      drawPass(PassDraws.dpoitTransparent, group);
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
      prepareCull,
      get gpuCulledCount() {
        return gpuCulledCount;
      },

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

        if (props.drawLists !== undefined && props.drawLists !== p.drawLists) {
          p.drawLists = props.drawLists;
          resetDrawLists();
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
      setOcclusionPyramid: (get: (() => OcclusionPyramid | undefined) | null) => {
        getOcclusionPyramid = get;
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
