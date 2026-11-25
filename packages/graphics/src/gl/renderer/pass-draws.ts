/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

// The renderer's passes as tables: per pass what it draws (the checks on an object's facts and state), how (its draws
// and their render states, per object its cull and front face) and its solid interior caps

import { DrawFact, DrawState, getCheckFacts } from '../draw-table.js';
import type { GraphicsRenderable } from '../renderable.js';
import type { GraphicsRenderVariant } from '../render-item.js';
import {
  BlendState,
  type CullMode,
  type DepthState,
  type FrontFace,
  RenderState,
  type StencilState,
} from '../render-state.js';

export const enum Flag {
  None = 0,
  BlendedFront = 1,
  BlendedBack = 2,
  DepthBack = 3,
  SolidInteriorMark = 4,
  SolidInteriorFill = 5,
}

export type SolidInteriorMode = 'opaque' | 'back' | 'blended' | 'oit';

// depth is reversed (near 1, far 0, see `DepthRangeTag`): nearer is greater
const DepthNearer: DepthState = { compare: 'greater', write: true };
const DepthNearerReadOnly: DepthState = { compare: 'greater', write: false };

export const OpaqueState: RenderState = { depth: DepthNearer };
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
export function getObjectState(
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
export function getSolidInteriorMarkState(passState: RenderState) {
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
export function getSolidInteriorFillStates(passState: RenderState, mode: SolidInteriorMode): SolidInteriorFillStates {
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
export type FactTest = (facts: number, state: number) => boolean;

const OpaqueBlockers =
  DrawFact.AlphaOne |
  DrawFact.TransparencyFull |
  DrawFact.DirectVolume |
  DrawFact.FuzzyPoints |
  DrawFact.Xray |
  DrawFact.GaussianSplats;
const TransparentReasons =
  DrawFact.AlphaBelowOne |
  DrawFact.TransparencyAny |
  DrawFact.DirectVolume |
  DrawFact.FuzzyPoints |
  DrawFact.TextOrImage |
  DrawFact.Xray |
  DrawFact.Occludable;

/**
 * alpha 1, not fully transparent, not a direct volume, fuzzy points, x-ray shaded or gaussian splats; or back faces
 * opaque (not splats)
 */
const checkOpaque: FactTest = (f) =>
  (f & OpaqueBlockers) === DrawFact.AlphaOne ||
  ((f & DrawFact.BackfacesOpaque) !== 0 && (f & DrawFact.GaussianSplats) === 0);

/**
 * not alpha 0, and alpha below 1, transparency, a direct volume, fuzzy points, text, an image or x-ray shaded; not
 * gaussian splats (their own pass, `checkGaussianSplats`)
 */
const checkTransparent: FactTest = (f) =>
  (f & (DrawFact.AlphaZero | DrawFact.GaussianSplats)) === 0 && (f & TransparentReasons) !== 0;

/** gaussian splats, not alpha 0: always blended, back to front */
const checkGaussianSplats: FactTest = (f) =>
  (f & (DrawFact.AlphaZero | DrawFact.GaussianSplats)) === DrawFact.GaussianSplats;

const checkPickable: FactTest = (_, s) => (s & DrawState.ColorOnly) === 0;

const checkMarkingDepth: FactTest = (f) =>
  (f & (DrawFact.AlphaZero | DrawFact.TransparencyFull | DrawFact.MarkerFull)) === 0;

const checkMarkingMask: FactTest = (f) => (f & DrawFact.MarkerAny) !== 0;

const checkClipGhost: FactTest = (f) => (f & DrawFact.ClipGhost) !== 0;

const checkEmissive: FactTest = (f) => (f & DrawFact.Emissive) !== 0;

const isDirectVolume: FactTest = (f) => (f & DrawFact.DirectVolume) !== 0;

const isGaussianSplats: FactTest = (f) => (f & DrawFact.GaussianSplats) !== 0;

/** caps are drawn from back faces, which transparent passes discard when `transparentBackfaces` is 'off' */
const checkTransparentCap: FactTest = (f, s) => checkTransparent(f, s) && (f & DrawFact.BackfacesOff) === 0;

const checkAll: FactTest = () => true;

/** The facts of an object outside a draw table, for the checks (`FactTest`; its state: `getDrawState`) */
export function getObjectFacts(r: GraphicsRenderable) {
  return getCheckFacts(r.values);
}

/**
 * whether an object with `state` draws in `pass`: the pick passes (the pick and its depth) skip objects that are not
 * pickable, the others pick-only objects (hit areas)
 */
export function isPassDrawn(pass: PassDraws, state: number) {
  return pass.pick ? (state & DrawState.Pickable) !== 0 : (state & DrawState.PickOnly) === 0;
}

/** A draw of an object in a pass */
type ObjectDraw = { readonly flag: Flag; readonly passState: RenderState };

/** The draws of a renderer pass, shared by drawing (`renderPass`) and collecting its pipeline requests (`collect`) */
export interface PassDraws {
  readonly variant: GraphicsRenderVariant;
  /** a pick pass (the pick or its depth): pickable objects, pick-only ones included (`isPassDrawn`) */
  readonly pick?: boolean;
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

export const PassDraws = {
  pick: { ...simplePassDraws('pick', checkPickable, OpaqueState, { mode: 'opaque' }), pick: true },
  pickDepth: { ...simplePassDraws('depth', checkPickable, OpaqueState, { mode: 'opaque' }), pick: true },
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
  // no caps: the ghost is the removed parts' surfaces
  clipGhost: { variant: 'depth', check: checkClipGhost, draws: singleDraw(OpaqueState) } as PassDraws,
  markingDepth: simplePassDraws('marking', checkMarkingDepth, OpaqueState, { mode: 'opaque' }),
  markingMask: simplePassDraws('marking', checkMarkingMask, OpaqueState, { mode: 'opaque' }),
  emissiveOpaque: createEmissiveOpaquePassDraws(false),
  emissiveOpaqueOccluded: createEmissiveOpaquePassDraws(true),
  emissiveTransparent: {
    ...simplePassDraws('emissive', (f, s) => checkTransparent(f, s) && checkEmissive(f, s), EmissiveTransparentState, {
      check: (f, s) => checkTransparentCap(f, s) && checkEmissive(f, s),
      mode: 'blended',
    }),
    select: (f, s) => !isDirectVolume(f, s),
  },
  tracing: simplePassDraws('tracing', checkOpaque, OpaqueState, { mode: 'opaque' }),
  tracingSecond: { variant: 'tracing', check: checkOpaque, draws: singleDraw(OpaqueState) } as PassDraws,
  opaque: simplePassDraws('color', checkOpaque, OpaqueState, { mode: 'opaque' }),
  opaqueSecond: { variant: 'color', check: checkOpaque, draws: singleDraw(OpaqueState) } as PassDraws,
  blendedTransparent: createBlendedTransparentPassDraws(false),
  blendedTransparentBackground: createBlendedTransparentPassDraws(true),
  volume: { variant: 'color', select: isDirectVolume, check: checkAll, draws: singleDraw(VolumeState) },
  // after the other transparency (its composite), on the opaque depth: as blended transparency, no faces nor caps
  gaussianSplats: {
    variant: 'color',
    select: isGaussianSplats,
    check: checkGaussianSplats,
    draws: singleDraw(BlendedTransparentState),
  } as PassDraws,
  gaussianSplatsBackground: {
    variant: 'color',
    select: isGaussianSplats,
    check: checkGaussianSplats,
    draws: singleDraw(BlendedTransparentBackgroundState),
  } as PassDraws,
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
        | 'clipGhost'
        | 'markingDepth'
        | 'markingMask'
        | 'emissiveTransparent'
        | 'tracing'
        | 'opaque'
        | 'blendedTransparent'
        | 'volume'
        | 'gaussianSplats'
        | 'wboitTransparent'
        | 'dpoitTransparent';
    };

export function getPassDraws(pass: Pass, transparentBackground: boolean): PassDraws {
  switch (pass.name) {
    case 'pick':
      return pass.variant === 'pick' ? PassDraws.pick : PassDraws.pickDepth;
    case 'emissiveOpaque':
      return pass.occludeWithOpaqueDepth ? PassDraws.emissiveOpaqueOccluded : PassDraws.emissiveOpaque;
    case 'blendedTransparent':
      return transparentBackground ? PassDraws.blendedTransparentBackground : PassDraws.blendedTransparent;
    case 'gaussianSplats':
      return transparentBackground ? PassDraws.gaussianSplatsBackground : PassDraws.gaussianSplats;
    default:
      return PassDraws[pass.name];
  }
}
