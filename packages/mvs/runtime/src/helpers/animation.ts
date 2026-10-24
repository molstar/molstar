/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { murmurHash3_128 } from '@molstar/core/data/util/hash-functions';
import { EasingFunctions } from '@molstar/core/math/easing';
import { clamp, lerp } from '@molstar/core/math/interpolate';
import { EPSILON, Mat3, Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { Color } from '@molstar/core/util/color';
import { deepEqual } from '@molstar/core/util';
import { ObjectMutator } from '@molstar/core/util/object-mutator';
import { areDistinctCanonical, CanonicalValues } from './canonical-values.js';
import { makeValuePaletteFunction, type ValuePaletteProps } from '@molstar/core/util/color/value-palette';
import { SplitColorProp } from '@molstar/mvs/components/split-uniform-color-theme';
import { SplitColor } from './utils.js';
import { palettePropsFromMVSPalette } from '@molstar/mvs/load-helpers';
import type { Snapshot } from '@molstar/mvs-builder/mvs-data';
import { type MVSAnimationNode, MVSAnimationSchema } from '@molstar/mvs-builder/tree/animation/animation-tree';
import type { Tree } from '@molstar/mvs-builder/tree/generic/tree-schema';
import { addDefaults } from '@molstar/mvs-builder/tree/generic/tree-utils';
import type { MVSTree } from '@molstar/mvs-builder/tree/mvs/mvs-tree';
import { ColorT, EasingT } from '@molstar/mvs-builder/tree/mvs/param-types';

/** Expand the animation of `snapshot` into MVS trees, one per distinct frame, each with its duration (`animation.md §1`). */
export async function generateStateTransition(
  ctx: RuntimeContext,
  snapshot: Snapshot,
  snapshotIndex: number,
  snapshotCount: number,
) {
  if (!snapshot.animation) return undefined;

  const tree = addDefaults(snapshot.animation, MVSAnimationSchema);
  const transitions = tree.children?.filter((child) => child.kind === 'interpolate');
  if (!transitions?.length) return undefined;

  const duration = Math.max(
    snapshot.animation.params?.duration_ms ?? 0,
    ...transitions.map((t) => (t.params.start_ms ?? 0) + t.params.duration_ms),
  );

  const frames: [tree: MVSTree, time: number][] = [];
  const dt = tree.params?.frame_time_ms ?? 1000 / 60;
  const N = Math.ceil(duration / dt);

  const { mutator, groups } = compileTransitions(snapshot.root, transitions);

  let prevRoot: MVSTree | undefined;
  for (let i = 0; i <= N; i++) {
    const t = i * dt;
    for (const group of groups) evaluateGroup(group, t, mutator);
    // (the mutator returns the previous frame's root when no value changed)
    const root = mutator.commit();

    if (root === prevRoot) {
      frames[frames.length - 1][1] += dt;
    } else {
      frames.push([root, dt]);
    }

    prevRoot = root;

    if (ctx.shouldUpdate) {
      await ctx.update({
        message: `Generating transition for snapshot ${snapshotIndex + 1}/${snapshotCount}`,
        current: i + 1,
        max: N,
      });
    }
  }

  return { tree, frametimeMs: dt, frames };
}

const EasingFnMap: Record<EasingT, (t: number) => number> = EasingFunctions;

interface PaletteFunction {
  hasSecondary: boolean;
  apply: (value: number, isSecondary: boolean) => Color;
}
const PaletteFunction = {
  /** Apply palette to value and return the resulting color as hex string. If the palette uses secondary colors, return primary and secondary color as "#ffffff/#eeeeee" */
  getAsHexStyle(palette: PaletteFunction, value: number): ColorT {
    if (palette.hasSecondary) {
      return `${Color.toHexStyle(palette.apply(value, false))}/${Color.toHexStyle(palette.apply(value, true))}`;
    } else {
      return Color.toHexStyle(palette.apply(value, false)) as ColorT;
    }
  },
};

interface InterpolationCacheEntry {
  paletteFn?: PaletteFunction;
  startColors?: DecodedSplitColors;
  endColors?: DecodedSplitColors;
  colorMap?: ColorMapPlan;
  /** Shared by the transitions of the animation (`animation.md §4`): released with them after the expansion */
  colorCaches: ColorCaches;
  rotation?: { axis: Vec3; angle: number; start: Quat; end: Quat };
}

function getTransitionKey(transition: MVSAnimationNode<'interpolate'>) {
  const prop = transition.params.property;
  if (Array.isArray(prop)) {
    return `${transition.params.target_ref}:${prop.join('.')}`;
  }
  return `${transition.params.target_ref}:${prop}`;
}

function groupTransitions(transitions: MVSAnimationNode<'interpolate'>[]) {
  const map = new Map<string, MVSAnimationNode<'interpolate'>[]>();
  const groups: MVSAnimationNode<'interpolate'>[][] = [];
  for (const t of transitions) {
    const key = getTransitionKey(t);
    if (!map.has(key)) {
      const group: MVSAnimationNode<'interpolate'>[] = [];
      map.set(key, group);
      groups.push(group);
    }
    map.get(key)!.push(t);
  }
  for (const group of groups) {
    group.sort((a, b) => {
      const s = (a.params.start_ms ?? 0) - (b.params.start_ms ?? 0);
      if (s !== 0) return s;
      return a.params.duration_ms - b.params.duration_ms;
    });
  }
  return groups;
}

/** The transitions of one property of one node (same `target_ref` and `property`), sorted by start. */
interface TransitionGroup {
  /** Where the property is in the snapshot tree. */
  slot: number;
  /** The property's value in the frame being built. */
  current: () => any;
  transitions: MVSAnimationNode<'interpolate'>[];
  startMs: number[];
  caches: InterpolationCacheEntry[];
}

/** Resolve each group's target once into a slot of an `ObjectMutator` over the snapshot tree. */
function compileTransitions(root: MVSTree, transitions: MVSAnimationNode<'interpolate'>[]) {
  const colorCaches = createColorCaches();
  const nodeMap = makeNodeMap(root, new Map(), []);
  const plan = ObjectMutator.plan(root);
  const groups: TransitionGroup[] = [];
  for (const transitionGroup of groupTransitions(transitions)) {
    const pivot = transitionGroup[0];
    const nodePath = nodeMap.get(pivot.params.target_ref);
    if (!nodePath) continue;

    const property = pivot.params.property;
    const isCustom = property[0] === 'custom';
    const node = select(root, nodePath, 0);
    const target = isCustom ? node?.custom : node?.params;
    // (a target missing in the snapshot tree: the group is skipped, as when it was resolved on every frame)
    if (!target) continue;

    plan.pushPrefix(...nodePath, isCustom ? 'custom' : 'params');
    const slot = plan.slot(typeof property === 'string' ? property : property.slice(isCustom ? 1 : 0));
    plan.popPrefix();
    groups.push({
      slot,
      current: undefined as any,
      transitions: transitionGroup,
      startMs: transitionGroup.map((t) => t.params.start_ms ?? 0),
      caches: transitionGroup.map(() => ({ colorCaches })),
    });
  }
  // (color maps are made once per content: two different ones are different, `animation.md §5`)
  const mutator = ObjectMutator.create(plan, { equals: (a, b) => !areDistinctCanonical(a, b) && deepEqual(a, b) });
  for (const group of groups) group.current = () => mutator.get(group.slot);
  return { mutator, groups };
}

/** Set the value of the group's property at `time`. */
function evaluateGroup(group: TransitionGroup, time: number, mutator: ObjectMutator<MVSTree>) {
  let index = 0;
  for (let i = group.transitions.length - 1; i > 0; i--) {
    if (time >= group.startMs[i]) {
      index = i;
      break;
    }
  }
  const transition = group.transitions[index];
  const previous = index > 0 ? group.transitions[index - 1] : undefined;
  const cacheEntry = group.caches[index];

  const startTime = group.startMs[index];
  const durationMs: number = transition.params.duration_ms ?? 0;
  const t = (time - startTime) / durationMs;

  const current = group.current;
  const next =
    transition.params.kind === 'transform_matrix'
      ? processTransformMatrix(transition, current, clamp(t, 0, 1), cacheEntry, previous)
      : processScalarLike(transition, current, t, cacheEntry, previous);
  if (next !== undefined) mutator.set(group.slot, next);
}

function applyFrequency(t: number, frequency: number, alternate: boolean) {
  let v = t * (frequency || 1);
  if (v < 1) return v;

  if (!alternate) {
    v = v % 1;
    if (v === 0) return 1;
    return v;
  }

  if (Math.abs(v - 1) < EPSILON) return 1;
  v = v % 2;
  if (v > 1) return 2 - v;
  return v;
}

function getPreviousScalarEnd(previous: MVSAnimationNode<'interpolate'> | undefined) {
  if (!previous || previous.params.kind === 'transform_matrix') return undefined;
  return previous.params.end;
}

function processScalarLike(
  transition: MVSAnimationNode<'interpolate'>,
  current: () => any,
  time: number,
  cacheEntry: InterpolationCacheEntry,
  previous: MVSAnimationNode<'interpolate'> | undefined,
) {
  if (transition.params.kind === 'transform_matrix') return;
  if (previous && previous.params.kind === 'transform_matrix') return;

  const startValue = transition.params.start ?? getPreviousScalarEnd(previous) ?? current();
  if (transition.params.kind === 'color' && !cacheEntry.paletteFn) {
    cacheEntry.paletteFn = makePaletteFunction(transition);
  }

  const endValue: any = transition.params.end;

  if (time <= 0 && transition.params.kind !== 'color') return startValue;
  if (time >= 1 - EPSILON && !transition.params.alternate_direction && transition.params.kind !== 'color')
    return endValue;

  let t = clamp(time, 0, 1);
  t = applyFrequency(t, transition.params.frequency ?? 1, !!transition.params.alternate_direction);

  const easing = EasingFnMap[transition.params.easing ?? 'linear'] ?? EasingFnMap['linear'];
  t = easing(t);

  if (transition.params.kind === 'scalar') {
    return interpolateScalars(
      startValue,
      endValue,
      t,
      transition.params.noise_magnitude ?? 0,
      !!transition.params.discrete,
    );
  } else if (transition.params.kind === 'vec3') {
    return interpolateVectors(
      startValue,
      endValue,
      t,
      transition.params.noise_magnitude ?? 0,
      !!transition.params.spherical,
    );
  } else if (transition.params.kind === 'rotation_matrix') {
    return interpolateRotation(startValue, endValue, t, transition.params.noise_magnitude ?? 0, cacheEntry);
  } else if (transition.params.kind === 'color') {
    if (cacheEntry.paletteFn) {
      return PaletteFunction.getAsHexStyle(cacheEntry.paletteFn, t);
    }

    const baseColors = typeof startValue === 'object' ? current() : undefined;
    return interpolateColors(startValue, endValue, t, cacheEntry, baseColors);
  }
}

function getPreviousMatrixEnd(
  previous: MVSAnimationNode<'interpolate'> | undefined,
  prop: 'rotation_start' | 'translation_start' | 'scale_start',
) {
  if (!previous || previous.params.kind !== 'transform_matrix') return undefined;
  return previous.params[prop];
}

const TransformState = {
  pivotTranslation: Mat4(),
  pivotTranslationInv: Mat4(),
  rotation: Mat4(),
  scale: Mat4(),
  translation: Mat4(),
  pivotNeg: Vec3(),
  temp: Mat4(),
};
function processTransformMatrix(
  transition: MVSAnimationNode<'interpolate'>,
  current: () => any,
  time: number,
  cache: InterpolationCacheEntry,
  previous: MVSAnimationNode<'interpolate'> | undefined,
) {
  if (transition.params.kind !== 'transform_matrix') return;
  if (previous && previous.params.kind !== 'transform_matrix') return;

  const transform = current() ?? Mat4.identity();

  const startRotation =
    transition.params.rotation_start ??
    getPreviousMatrixEnd(previous, 'rotation_start') ??
    Mat3.fromMat4(Mat3(), transform);
  const startTranslation =
    transition.params.translation_start ??
    getPreviousMatrixEnd(previous, 'translation_start') ??
    Mat4.getTranslation(Vec3(), transform);
  const startScale =
    transition.params.scale_start ??
    getPreviousMatrixEnd(previous, 'scale_start') ??
    Mat4.getScaling(Vec3(), transform);

  const endRotation = transition.params.rotation_end;
  const endTranslation = transition.params.translation_end;
  const endScale = transition.params.scale_end;

  let rotation, translation, scale;

  if (time <= 0) {
    rotation = startRotation as Mat3;
    translation = startTranslation as Vec3;
    scale = startScale as Vec3;
  } else {
    const clampedTime = clamp(time, 0, 1);

    let t = applyFrequency(
      clampedTime,
      transition.params.rotation_frequency ?? 1,
      !!transition.params.rotation_alternate_direction,
    );
    let easing = EasingFnMap[transition.params.rotation_easing ?? 'linear'] ?? EasingFnMap['linear'];
    rotation = interpolateRotation(
      startRotation as Mat3,
      endRotation as Mat3,
      easing(t),
      transition.params.rotation_noise_magnitude ?? 0,
      cache,
    );

    t = applyFrequency(
      clampedTime,
      transition.params.translation_frequency ?? 1,
      !!transition.params.translation_alternate_direction,
    );
    easing = EasingFnMap[transition.params.translation_easing ?? 'linear'] ?? EasingFnMap['linear'];
    translation = interpolateVec3(
      startTranslation as Vec3,
      endTranslation as Vec3 | undefined,
      easing(t),
      transition.params.translation_noise_magnitude ?? 0,
      false,
    );

    t = applyFrequency(
      clampedTime,
      transition.params.scale_frequency ?? 1,
      !!transition.params.scale_alternate_direction,
    );
    easing = EasingFnMap[transition.params.scale_easing ?? 'linear'] ?? EasingFnMap['linear'];
    scale = interpolateVec3(
      startScale as Vec3,
      endScale as Vec3 | undefined,
      easing(t),
      transition.params.scale_noise_magnitude ?? 0,
      false,
    );
  }

  const pivot = transition.params.pivot ?? Vec3.zero();

  Mat4.fromTranslation(TransformState.translation, translation);
  Mat4.fromScaling(TransformState.scale, scale);
  Mat4.setIdentity(TransformState.rotation);
  Mat4.fromMat3(TransformState.rotation, rotation);
  Mat4.fromTranslation(TransformState.pivotTranslation, pivot as Vec3);
  Mat4.fromTranslation(TransformState.pivotTranslationInv, Vec3.negate(TransformState.pivotNeg, pivot as Vec3));

  // translation . pivot . rotation . scale . pivotInv
  const result = Mat4();
  Mat4.mul(result, TransformState.scale, TransformState.pivotTranslationInv);
  Mat4.mul(result, TransformState.rotation, result);
  Mat4.mul(result, TransformState.translation, result);

  return result;
}

function interpolateScalars(
  start: number | number[],
  end: number | number[] | undefined,
  t: number,
  noise: number,
  discrete: boolean,
) {
  if (Array.isArray(start)) {
    const ret = Array.from<number>({ length: start.length }).fill(0.1);
    if (!end || !Array.isArray(end)) {
      for (let i = 0; i < start.length; i++) {
        ret[i] = interpolateScalar(start[i], end, t, noise, discrete);
      }
      return ret;
    }

    for (let i = 0; i < start.length; i++) {
      ret[i] = interpolateScalar(start[i], end[i], t, noise, discrete);
    }
    return ret;
  }

  if (Array.isArray(end)) {
    const ret = Array.from<number>({ length: end.length }).fill(0.1);
    for (let i = 0; i < end.length; i++) {
      ret[i] = interpolateScalar(start, end[i], t, noise, discrete);
    }
    return ret;
  }

  return interpolateScalar(start, end, t, noise, discrete);
}

function interpolateScalar(start: number, end: number | undefined, t: number, noise: number, discrete: boolean) {
  let v = typeof end === 'number' ? lerp(start, end, t) : start;
  if (noise) {
    v += (Math.random() - 0.5) * noise;
  }
  if (discrete) {
    v = Math.round(v);
  }
  return v;
}

const InterpolateVectorsState = {
  start: Vec3(),
  end: Vec3(),
  v: Vec3(),
};
function interpolateVectors(
  start: number[],
  end: number[] | undefined,
  t: number,
  noise: number,
  isSpherical: boolean,
) {
  if ((!end || start === end) && !noise) return start;

  const ret: number[] = Array.from<number>({ length: start.length }).fill(0.1);

  for (let i = 0; i < start.length; i += 3) {
    const s = Vec3.fromArray(InterpolateVectorsState.start, start, i);

    let v: Vec3;
    if (end) {
      const e = Vec3.fromArray(InterpolateVectorsState.end, end, i);
      v = isSpherical ? Vec3.slerp(InterpolateVectorsState.v, s, e, t) : Vec3.lerp(InterpolateVectorsState.v, s, e, t);
    } else {
      v = Vec3.clone(s);
    }

    if (noise && t <= 1 - EPSILON) {
      Vec3.random(Vec3Noise, noise);
      Vec3.add(v, v, Vec3Noise);
    }

    Vec3.toArray(v, ret, i);
  }

  return ret;
}

const Vec3Noise = Vec3();
function interpolateVec3(start: Vec3, end: Vec3 | undefined, t: number, noise: number, isSpherical: boolean) {
  if ((!end || start === end) && !noise) return start;

  let v: Vec3;

  if (end) {
    v = isSpherical ? Vec3.slerp(Vec3(), start, end, t) : Vec3.lerp(Vec3(), start, end, t);
  } else {
    v = Vec3.clone(start);
  }

  if (noise && t <= 1 - EPSILON) {
    Vec3.random(Vec3Noise, noise);
    Vec3.add(v, v, Vec3Noise);
  }
  return v;
}

const RotationState = {
  start: Quat(),
  end: Quat(),
  v: Quat(),
  noise: Quat(),
  axis: Vec3(),
  temp: Mat4(),
};
function interpolateRotation(
  start: Mat3,
  end: Mat3 | undefined,
  t: number,
  noise: number,
  cache: InterpolationCacheEntry,
) {
  if ((!end || start === end) && !noise) return start;

  if (end) {
    if (!cache.rotation) {
      cache.rotation = {
        ...relativeAxisAngle(start, end),
        start: Quat.fromMat3(Quat(), start),
        end: Quat.fromMat3(Quat(), end),
      };
    }

    const { axis, angle, start: startQ, end: endQ } = cache.rotation;

    if (angle < 1e-6) {
      // start ≈ end: make a clean spin about the detected (or default) axis
      Quat.setAxisAngle(RotationState.v, axis, t * 2 * Math.PI); // Make a full turn
    } else {
      // Normal case: stick with your existing slerp between start/end
      Quat.slerp(RotationState.v, startQ, endQ, t);
    }
  } else {
    Quat.fromMat3(RotationState.v, start);
  }

  if (noise && t <= 1 - EPSILON) {
    Vec3.random(RotationState.axis, 1);
    Quat.setAxisAngle(RotationState.noise, RotationState.axis, 2 * Math.PI * noise * (Math.random() - 0.5));
    Quat.multiply(RotationState.v, RotationState.noise, RotationState.v);
  }
  Mat4.fromQuat(RotationState.temp, RotationState.v);
  return Mat3.fromMat4(Mat3(), RotationState.temp);
}

const _fallbackSplitColor: [Color, Color] = [Color(0), Color(0)]; // invalid color strings fall back to black rather then undefined to make debugging easier
type DecodedSplitColors =
  | { kind: 'dict'; dict: Record<string, [Color, Color]> }
  | { kind: 'value'; value: [Color, Color] };

/** Decoded colors by string (color maps repeat a few colors many times); the same tuple for the same string. */
const DecodedColors = new Map<unknown, [Color, Color]>();
const DecodedColorsStrict = new Map<unknown, [Color, Color]>();
function decodeCached(cache: Map<unknown, [Color, Color]>, decode: (c: any) => [Color, Color] | undefined, c: unknown) {
  let decoded = cache.get(c);
  if (!decoded) {
    if (cache.size > 10_000) cache.clear();
    decoded = decode(c) ?? _fallbackSplitColor;
    cache.set(c, decoded);
  }
  return decoded;
}

function decodeSplitColors(
  color: ColorT | Record<number | string, ColorT> | undefined,
  baseColors: Record<number | string, ColorT> | undefined,
): DecodedSplitColors | undefined {
  if (color === undefined || color === null) return undefined;

  if (typeof color === 'object') {
    const ret: Record<number | string, [Color, Color]> = {};
    if (baseColors) {
      for (const key of Object.keys(baseColors)) {
        ret[key] = decodeCached(DecodedColors, (c) => SplitColor.decode(c), baseColors[key]);
      }
    }
    for (const key of Object.keys(color)) {
      ret[key] = decodeCached(DecodedColorsStrict, (c) => SplitColor.decodeStrict(c), color[key]);
    }
    return { kind: 'dict', dict: ret };
  }

  return { kind: 'value', value: SplitColor.decode(color) ?? _fallbackSplitColor };
}

function interpolateColors(
  start: ColorT | Record<number | string, ColorT> | undefined,
  end: ColorT | Record<number | string, ColorT> | undefined,
  time: number,
  cacheEntry: InterpolationCacheEntry,
  baseColors: Record<number, ColorT> | undefined,
) {
  const t = clamp(time, 0, 1);

  if (cacheEntry.paletteFn) {
    return PaletteFunction.getAsHexStyle(cacheEntry.paletteFn, t);
  }

  const caches = cacheEntry.colorCaches;
  if (cacheEntry.colorMap) return ColorMapPlan.at(caches, cacheEntry.colorMap, t, baseColors!);
  if (isColorMap(start) && isColorMap(end) && isColorMap(baseColors)) {
    const plan = ColorMapPlan.fromRaw(caches, baseColors, start, end);
    if (plan) return ColorMapPlan.at(caches, (cacheEntry.colorMap = plan), t, baseColors);
  }

  const startColors = (cacheEntry.startColors ??= decodeSplitColors(start, baseColors));
  const endColors = (cacheEntry.endColors ??= decodeSplitColors(end, undefined));

  if (startColors?.kind === 'dict') {
    if (baseColors === undefined || baseColors === null) throw new Error('Missing color value in the target node');
    if (typeof baseColors !== 'object') throw new Error('Cannot interpolate from scalar color to color mapping');
    const plan = (cacheEntry.colorMap ??= ColorMapPlan.get(caches, baseColors, startColors.dict, endColors));
    return ColorMapPlan.at(caches, plan, t, baseColors);
  } else if (startColors?.kind === 'value') {
    if (endColors?.kind === 'dict') {
      throw new Error('Cannot interpolate from scalar color to color mapping');
    } else if (endColors?.kind === 'value') {
      return interpolateSplitColor(startColors.value, endColors.value, t);
    } else {
      return start;
    }
  } else {
    // no start value (a key the target's color map does not have, no `start`): nothing to interpolate, the property is
    // left as it is (as before split coloring, which made this an error and broke documents that relied on it)
    return undefined;
  }
}

/**
 * A color map transition (`animation.md §4`): the keys in the order they are set, each either a constant color or one of
 * a few groups of keys that go from the same start to the same end color. Prepared once per content (base map object,
 * start, end), so transitions with the same colors share it; the maps it returns are kept by their group colors, so
 * equal maps are the same object (frame dedup, hashing, comparisons and memory by identity).
 */
interface ColorMapPlan {
  base: Record<number | string, ColorT>;
  keys: string[];
  /** Hash of `keys` (their order) */
  keysId: string;
  /** Per key: its constant color, or undefined when it is in a group */
  constants: (ColorT | undefined)[];
  /** Per key: its group, or -1 */
  group: Int32Array;
  groups: { start: [Color, Color]; end: [Color, Color] }[];
  maps: Map<string, Record<number | string, ColorT>>;
}

/** The color map caches of one animation, by base map (`animation.md §4`). */
interface ColorCaches {
  /** The plans made over a base map and its key orders */
  plans: WeakMap<object, { plans: ColorMapPlan[]; keyOrders: string[][] }>;
  /** The start and end maps a plan was made from */
  raw: WeakMap<object, RawColorMapEntry[]>;
  /** The maps made, by content: equal maps from different transitions are the same object */
  byContent: WeakMap<object, Map<string, Record<number | string, ColorT>>>;
}

function createColorCaches(): ColorCaches {
  return { plans: new WeakMap(), raw: new WeakMap(), byContent: new WeakMap() };
}
interface RawColorMapEntry {
  start: Record<string, unknown>;
  end: Record<string, unknown>;
  startKeys: string[];
  endKeys: string[];
  plan: ColorMapPlan;
}

/** The same keys in the same order with the same values. */
function rawMapsEqual(ka: string[], a: Record<string, unknown>, kb: string[], b: Record<string, unknown>): boolean {
  if (a === b) return true;
  if (ka.length !== kb.length) return false;
  for (let i = 0; i < ka.length; i++) {
    if (ka[i] !== kb[i] || a[ka[i]] !== b[kb[i]]) return false;
  }
  return true;
}

function arraysEqual<T>(a: ArrayLike<T>, b: ArrayLike<T>): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

function tupleEqual(a: [Color, Color], b: [Color, Color]) {
  return a === b || (a[0] === b[0] && a[1] === b[1]);
}

function isColorMap(x: unknown): x is Record<number | string, ColorT> {
  return typeof x === 'object' && x !== null;
}

type AddColorMapKey = (key: string, constant: ColorT | undefined, from?: [Color, Color], to?: [Color, Color]) => void;

const ColorMapPlan = {
  get(
    caches: ColorCaches,
    base: Record<number | string, ColorT>,
    start: Record<string, [Color, Color]>,
    end: DecodedSplitColors | undefined,
  ): ColorMapPlan {
    // the keys in the order the earlier loops set them: start keys not in an end map, then the end map's keys
    return ColorMapPlan.make(caches, base, (add) => {
      if (end?.kind === 'dict') {
        for (const key of Object.keys(start)) {
          if (!(key in end.dict)) add(key, SplitColor.toHexStyle(...start[key]));
        }
        for (const key of Object.keys(end.dict)) {
          if (key in start) add(key, undefined, start[key], end.dict[key]);
          else add(key, SplitColor.toHexStyle(...end.dict[key]));
        }
      } else if (end?.kind === 'value') {
        for (const key of Object.keys(start)) add(key, undefined, start[key], end.value);
      } else {
        for (const key of Object.keys(start)) add(key, SplitColor.toHexStyle(...start[key]));
      }
    });
  },
  /**
   * The plan for these start and end maps as written: the plan of an earlier transition with equal maps (same keys in
   * the same order, equal values; documents repeat the same fades many times), else `fromSameKeys`.
   */
  fromRaw(
    caches: ColorCaches,
    base: Record<number | string, ColorT>,
    start: Record<number | string, ColorT>,
    end: Record<number | string, ColorT>,
  ): ColorMapPlan | undefined {
    let seen = caches.raw.get(base);
    if (!seen) caches.raw.set(base, (seen = []));
    const startKeys = Object.keys(start);
    const endKeys = Object.keys(end);
    for (let k = 0; k < seen.length; k++) {
      const entry = seen[k];
      if (
        rawMapsEqual(entry.endKeys, entry.end, endKeys, end) &&
        rawMapsEqual(entry.startKeys, entry.start, startKeys, start)
      ) {
        // (the most recent match first: fades repeat in runs)
        if (k > 0) seen.unshift(...seen.splice(k, 1));
        return entry.plan;
      }
    }
    const plan = ColorMapPlan.fromSameKeys(caches, base, start, end);
    if (plan) seen.unshift({ start, end, startKeys, endKeys, plan });
    return plan;
  },
  /**
   * `get` for start and end maps over the same keys, which the base's keys are among (a map fading as a whole): every key
   * goes from its start to its end color, in the end map's key order, read in one pass from the maps as they are.
   */
  fromSameKeys(
    caches: ColorCaches,
    base: Record<number | string, ColorT>,
    start: Record<number | string, ColorT>,
    end: Record<number | string, ColorT>,
  ): ColorMapPlan | undefined {
    const endKeys = Object.keys(end);
    let startCount = 0;
    for (const key in start) if (Object.hasOwn(start, key)) startCount++;
    if (startCount !== endKeys.length) return undefined;
    const n = endKeys.length;
    const start0 = n > 0 ? start[endKeys[0]] : undefined;
    const end0 = n > 0 ? end[endKeys[0]] : undefined;
    let uniform = n > 0;
    for (const key of endKeys) {
      if (!Object.hasOwn(start, key)) return undefined;
      if (uniform && (start[key] !== start0 || end[key] !== end0)) uniform = false;
    }
    for (const key in base) if (Object.hasOwn(base, key) && !Object.hasOwn(end, key)) return undefined;
    if (uniform) {
      // (one start and one end color for every key: one group)
      const from = decodeCached(DecodedColorsStrict, (c) => SplitColor.decodeStrict(c), start0);
      const to = decodeCached(DecodedColorsStrict, (c) => SplitColor.decodeStrict(c), end0);
      return ColorMapPlan.intern(caches, base, endKeys, new Array(n).fill(undefined), new Int32Array(n), [
        { start: from, end: to },
      ]);
    }
    return ColorMapPlan.make(caches, base, (add) => {
      for (const key of endKeys) {
        const from = decodeCached(DecodedColorsStrict, (c) => SplitColor.decodeStrict(c), start[key]);
        const to = decodeCached(DecodedColorsStrict, (c) => SplitColor.decodeStrict(c), end[key]);
        add(key, undefined, from, to);
      }
    });
  },
  make(
    caches: ColorCaches,
    base: Record<number | string, ColorT>,
    entries: (add: AddColorMapKey) => void,
  ): ColorMapPlan {
    const keys: string[] = [];
    const constants: (ColorT | undefined)[] = [];
    const groupOf: number[] = [];
    const groups: ColorMapPlan['groups'] = [];
    // (by the decoded tuples, the same object per color string: `decodeCached`)
    const groupIndex = new Map<[Color, Color], Map<[Color, Color], number>>();
    const add = (key: string, constant: ColorT | undefined, from?: [Color, Color], to?: [Color, Color]) => {
      keys.push(key);
      constants.push(constant);
      if (constant !== undefined) {
        groupOf.push(-1);
        return;
      }
      let byEnd = groupIndex.get(from!);
      if (!byEnd) groupIndex.set(from!, (byEnd = new Map()));
      let g = byEnd.get(to!);
      if (g === undefined) {
        g = groups.length;
        groups.push({ start: from!, end: to! });
        byEnd.set(to!, g);
      }
      groupOf.push(g);
    };
    entries(add);
    return ColorMapPlan.intern(caches, base, keys, constants, Int32Array.from(groupOf), groups);
  },
  /** The plan of this content over `base`: an existing one when equal (compared as arrays, no strings built). */
  intern(
    caches: ColorCaches,
    base: Record<number | string, ColorT>,
    keys: string[],
    constants: (ColorT | undefined)[],
    group: Int32Array,
    groups: ColorMapPlan['groups'],
  ): ColorMapPlan {
    let registry = caches.plans.get(base);
    if (!registry) caches.plans.set(base, (registry = { plans: [], keyOrders: [] }));
    let keysIndex = registry.keyOrders.findIndex((k) => arraysEqual(k, keys));
    if (keysIndex < 0) keysIndex = registry.keyOrders.push(keys) - 1;
    const keysId = `${keysIndex}`;
    for (const plan of registry.plans) {
      if (
        plan.keysId === keysId &&
        plan.groups.length === groups.length &&
        plan.groups.every((g, i) => tupleEqual(g.start, groups[i].start) && tupleEqual(g.end, groups[i].end)) &&
        arraysEqual(plan.group, group) &&
        arraysEqual(plan.constants, constants)
      ) {
        return plan;
      }
    }
    const plan: ColorMapPlan = { base, keys, keysId, constants, group, groups, maps: new Map() };
    registry.plans.push(plan);
    return plan;
  },
  /**
   * The map at eased time `t`, over `base`: the plan's own base (a map kept by content), or another one when the base
   * changes from step to step (a map slot inside another animated slot): a new map over it, as each step made before.
   */
  at(
    caches: ColorCaches,
    plan: ColorMapPlan,
    t: number,
    base: Record<number | string, ColorT>,
  ): Record<number | string, ColorT> {
    const colors = plan.groups.map((g) => interpolateSplitColor(g.start, g.end, t));
    if (base !== plan.base) {
      const map: Record<number | string, ColorT> = Array.isArray(base)
        ? (base.slice() as Record<number, ColorT>)
        : { ...base };
      for (let i = 0; i < plan.keys.length; i++) {
        map[plan.keys[i]] = plan.group[i] >= 0 ? colors[plan.group[i]] : plan.constants[i]!;
      }
      return map;
    }
    const id = colors.join('|');
    let map = plan.maps.get(id);
    if (map) return map;
    // the content: the keys' order, and each key's color as an index into the distinct colors
    const { keys, constants, group } = plan;
    const values: ColorT[] = new Array(keys.length);
    const distinct = new Map<ColorT, number>();
    let content = '';
    for (let i = 0; i < keys.length; i++) {
      const value = group[i] >= 0 ? colors[group[i]] : constants[i]!;
      values[i] = value;
      let d = distinct.get(value);
      if (d === undefined) distinct.set(value, (d = distinct.size));
      content += d + ',';
    }
    const contentId = `${plan.keysId}|${Array.from(distinct.keys()).join(',')}|${murmurHash3_128(content, 42)}`;
    let byContent = caches.byContent.get(plan.base);
    if (!byContent) caches.byContent.set(plan.base, (byContent = new Map()));
    map = byContent.get(contentId);
    if (!map) {
      map = Array.isArray(plan.base) ? (plan.base.slice() as Record<number, ColorT>) : { ...plan.base };
      for (let i = 0; i < keys.length; i++) map[keys[i]] = values[i];
      byContent.set(contentId, map);
      CanonicalValues.add(map);
    }
    plan.maps.set(id, map);
    return map;
  },
};

function interpolateSplitColor(start: [Color, Color], end: [Color, Color], t: number): ColorT {
  const out0 = Color.interpolate(start[0], end[0], t);
  if (start[0] === start[1] && end[0] === end[1]) {
    return SplitColor.toHexStyle(out0);
  } else {
    const out1 = Color.interpolate(start[1], end[1], t);
    return SplitColor.toHexStyle(out0, out1);
  }
}

function select(params: any, path: string | (string | number)[], offset: number) {
  if (typeof path === 'string') {
    return params?.[path];
  }

  let f = params;
  for (let i = offset; i < path.length; i++) {
    if (!f) break;
    f = f[path[i]];
  }

  return f;
}

function makeNodeMap(tree: Tree, map: Map<string, (string | number)[]>, currentPath: (string | number)[]) {
  if (tree.ref) {
    map.set(tree.ref, [...currentPath]);
  }

  if (!tree.children) return map;

  currentPath.push('children');
  for (let i = 0; i < tree.children.length; i++) {
    const child = tree.children[i];
    currentPath.push(i);
    makeNodeMap(child, map, currentPath);
    currentPath.pop();
  }
  currentPath.pop();

  return map;
}

function makePaletteFunction(props: MVSAnimationNode<'interpolate'>): PaletteFunction | undefined {
  if (props.params.kind !== 'color' || !props.params.palette) return undefined;

  const params = palettePropsFromMVSPalette(props.params.palette);
  if (params.name === 'discrete') return makePaletteFunctionNumeric(params);
  if (params.name === 'continuous') return makePaletteFunctionNumeric(params);
  throw new Error(`NotImplementedError: makePaletteFunction for ${(props as any).name}`);
}

/** The palette over the value clamped to [0, 1] (`absolute`: the value is already normalized), black where it gives no color. */
function makePaletteFunctionNumeric(
  params: Extract<ValuePaletteProps, { name: 'discrete' | 'continuous' }>,
): PaletteFunction {
  const defaultColor = Color(0x0);
  const palette = makeValuePaletteFunction(
    { ...params, params: { ...params.params, mode: 'absolute' } } as ValuePaletteProps,
    [],
  );
  const colors = params.params.colors.map((item) => SplitColorProp.toTuple(item.color));
  const flows =
    params.name === 'continuous'
      ? [SplitColorProp.toTuple(params.params.underflowColor), SplitColorProp.toTuple(params.params.overflowColor)]
      : [];
  const hasSecondary = colors.some((col) => col[0] !== col[1]) || flows.some((col) => col[0] !== col[1]);

  return {
    hasSecondary,
    apply: (value: number, isSecondary: boolean) => palette(clamp(value, 0, 1), isSecondary) ?? defaultColor,
  };
}

const RelativeAxisAngleState = {
  Rt: Mat3(),
  R: Mat3(),
};
function relativeAxisAngle(start: Mat3, end: Mat3): { axis: Vec3; angle: number } {
  // R_rel = end * start^T
  const R0 = start,
    R1 = end;
  const Rt = Mat3.transpose(RelativeAxisAngleState.Rt, R0);
  const R = Mat3.mul(RelativeAxisAngleState.R, R1, Rt);

  const tr = R[0] + R[4] + R[8]; // trace
  let angle = Math.acos(clamp((tr - 1) * 0.5, -1, 1)); // in [0, π]
  const axis = Vec3();

  const eps = 1e-6;
  const sinA = Math.sin(angle);

  if (angle < eps) {
    // Near identity: axis undefined; return any unit axis (choose something stable)
    Vec3.set(axis, 0, 0, 1);
    angle = 0.0;
    return { axis, angle };
  }

  if (Math.PI - angle > 1e-4) {
    // General case
    axis[0] = (R[5] - R[7]) / (2 * sinA); // (r32 - r23)
    axis[1] = (R[6] - R[2]) / (2 * sinA); // (r13 - r31)
    axis[2] = (R[1] - R[3]) / (2 * sinA); // (r21 - r12)
    Vec3.normalize(axis, axis);
    return { axis, angle };
  }

  // angle ~ π: use diagonal-based extraction for stability
  // Compute squared components then pick the largest to avoid precision loss
  const xx = Math.max(0, (R[0] + 1) * 0.5);
  const yy = Math.max(0, (R[4] + 1) * 0.5);
  const zz = Math.max(0, (R[8] + 1) * 0.5);

  let x = Math.sqrt(xx),
    y = Math.sqrt(yy),
    z = Math.sqrt(zz);

  if (x >= y && x >= z) {
    x = Math.max(x, 1e-8);
    y = (R[1] + R[3]) / (4 * x);
    z = (R[2] + R[6]) / (4 * x);
    Vec3.set(axis, x, y, z);
  } else if (y >= x && y >= z) {
    y = Math.max(y, 1e-8);
    x = (R[1] + R[3]) / (4 * y);
    z = (R[5] + R[7]) / (4 * y);
    Vec3.set(axis, x, y, z);
  } else {
    z = Math.max(z, 1e-8);
    x = (R[2] + R[6]) / (4 * z);
    y = (R[5] + R[7]) / (4 * z);
    Vec3.set(axis, x, y, z);
  }

  Vec3.normalize(axis, axis);
  return { axis, angle: Math.PI };
}
