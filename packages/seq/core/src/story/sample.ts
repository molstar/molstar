/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Sampling a plan (_spec/README.md "Stories", Sampling): pure, deterministic, allocation-free once warm (a plan's frame
// is kept and written again), the same function for the viewer's port, scrubbing and the headless export. What a frame
// shows over the shown state (`frameRoot`, `frameView`) is made from it, never committed.

import type { Op } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { colorFrom, type StoryState, type TransitionPlan, type ViewState } from './plan.js';
import { applyLoose, type StepViewport } from './state.js';

/** Where a channel is at a frame's time. */
export const ChannelPhase = { before: 0, moving: 1, done: 2 } as const;

/** A plan at a time: each channel's values (by its `slot`) and phase; reused by `sample`. */
export interface Frame {
  /** The plan's time, [0, 1]. */
  t: number;
  readonly values: Float64Array;
  /** By channel: `ChannelPhase`. A cut's is `done` once switched. */
  readonly phase: Uint8Array;
  /** Channels holding their start until their data arrives. */
  waiting: number;
}

export function createFrame(plan: TransitionPlan): Frame {
  return { t: 0, values: new Float64Array(plan.slots), phase: new Uint8Array(plan.channels.length), waiting: 0 };
}

const Frames = new WeakMap<TransitionPlan, Frame>();

/**
 * The plan at `t` (clamped to [0, 1]), each channel eased: written into `out`, else the plan's own frame (kept, so a
 * warm call allocates nothing). `sample(plan, 0)` is the start, `sample(plan, 1)` the end, exactly.
 */
export function sample(plan: TransitionPlan, t: number, out?: Frame): Frame {
  let f = out;
  if (!f) {
    f = Frames.get(plan);
    if (!f) Frames.set(plan, (f = createFrame(plan)));
  }
  const tt = t > 1 ? 1 : t >= 0 ? t : 0;
  f.t = tt;
  f.waiting = 0;
  const v = f.values;
  const cs = plan.channels;
  for (let i = 0; i < cs.length; i++) {
    const c = cs[i];
    // the channel's own time and its easing, inlined on purpose: a double a call returns from a site of many channel
    // shapes is boxed, and this loop must allocate nothing
    let u: number;
    if (c.pending) {
      u = 0;
      f.waiting++;
    } else if (c.t1 <= c.t0) u = tt >= c.t1 ? 1 : 0;
    else {
      u = (tt - c.t0) / (c.t1 - c.t0);
      u = u <= 0 ? 0 : u >= 1 ? 1 : u;
    }
    let e: number;
    if (u <= 0) e = 0;
    else if (u >= 1) e = 1;
    else if (c.easing === 'linear') e = u;
    else if (c.easing === 'ease-in') e = u * u * u;
    else if (c.easing === 'ease-out') e = 1 - (1 - u) * (1 - u) * (1 - u);
    else e = u < 0.5 ? 4 * u * u * u : 1 - ((2 - 2 * u) * (2 - 2 * u) * (2 - 2 * u)) / 2;
    f.phase[i] = u <= 0 ? ChannelPhase.before : u >= 1 ? ChannelPhase.done : ChannelPhase.moving;
    const s = c.slot;
    switch (c.kind) {
      case 'view': {
        const n = c.from.length;
        if (e >= 1) {
          for (let k = 0; k < n; k++) v[s + k] = c.to[k];
          break;
        }
        // the window: its width geometric (zooming reads evenly), its centre along it; out-across-in passes through
        // the width holding both
        const s0 = c.from[0];
        const e0 = c.from[1];
        const s1 = c.to[0] + c.shift;
        const e1 = c.to[1] + c.shift;
        const w0 = e0 - s0;
        const w1 = e1 - s1;
        let w: number;
        if (c.path === 'out-across-in') {
          // at the middle the window is centred between the two and holds both
          const mid = (s0 + e0 + s1 + e1) / 4;
          const span = 2 * Math.max(mid - Math.min(s0, s1), Math.max(e0, e1) - mid);
          w =
            w0 > 0 && w1 > 0 && span > 0
              ? e < 0.5
                ? w0 * (span / w0) ** (2 * e)
                : span * (w1 / span) ** (2 * e - 1)
              : w0 + (w1 - w0) * e;
        } else w = w0 > 0 && w1 > 0 ? w0 * (w1 / w0) ** e : w0 + (w1 - w0) * e;
        const c0 = (s0 + e0) / 2;
        const cm = c0 + ((s1 + e1) / 2 - c0) * e;
        v[s] = cm - w / 2;
        v[s + 1] = cm + w / 2;
        for (let k = 2; k < n; k++) {
          const a = c.from[k];
          const b = c.to[k];
          v[s + k] = c.geometric[k] ? a * (b / a) ** e : a + (b - a) * e;
        }
        break;
      }
      case 'number':
        // an angle ends where its short way round does (the same angle as `to`, so it never jumps back a turn)
        // two stores, not one of a conditional: a value merged from a field and a sum is boxed
        if (e >= 1 && c.hold === undefined) v[s] = c.to;
        else v[s] = c.from + c.delta * e;
        break;
      case 'color':
        if (e >= 1) {
          v[s] = c.b[0];
          v[s + 1] = c.b[1];
          v[s + 2] = c.b[2];
        } else {
          v[s] = c.a[0] + (c.b[0] - c.a[0]) * e;
          v[s + 1] = c.a[1] + (c.b[1] - c.a[1]) * e;
          v[s + 2] = c.a[2] + (c.b[2] - c.a[2]) * e;
        }
        break;
      case 'enter':
      case 'leave':
        v[s] = e;
        // with no height to grow or shrink, it comes at its window's start and goes at its end: what names it (a
        // track naming entering data) never dangles
        if (c.height === undefined)
          f.phase[i] = (c.kind === 'enter' ? u > 0 : u >= 1) ? ChannelPhase.done : ChannelPhase.before;
        break;
      default:
        // a cut switches at its time within its window
        v[s] = e;
        f.phase[i] = u > 0 && u >= c.cut ? ChannelPhase.done : ChannelPhase.before;
    }
  }
  return f;
}

/** Whether every channel of the frame is at its end. */
export function frameDone(plan: TransitionPlan, frame: Frame): boolean {
  if (frame.waiting) return false;
  for (let i = 0; i < plan.channels.length; i++) if (frame.phase[i] !== ChannelPhase.done) return false;
  return true;
}

/** Whether no channel of the frame has started. */
function frameAtStart(plan: TransitionPlan, frame: Frame): boolean {
  for (let i = 0; i < plan.channels.length; i++) if (frame.phase[i] !== ChannelPhase.before) return false;
  return true;
}

/**
 * The spec a frame shows: the start's root with each channel's value (a tween's, a node entering or going, a cut once
 * switched); the end's root exactly once every channel is done, the start's before any began.
 */
export function frameRoot(plan: TransitionPlan, frame: Frame): Node {
  if (frameDone(plan, frame)) return plan.to.root;
  if (frameAtStart(plan, frame)) return plan.from.root;
  const ops: Op[] = [];
  const v = frame.values;
  plan.channels.forEach((c, i) => {
    const phase = frame.phase[i];
    switch (c.kind) {
      case 'view':
        return;
      case 'number':
        ops.push({ op: 'set', at: c.at, key: c.key, value: v[c.slot], ...(c.where ? { where: c.where } : {}) });
        if (c.hold && phase === ChannelPhase.moving)
          ops.push({ op: 'set', at: c.at, key: c.hold.key, value: c.hold.value });
        return;
      case 'color':
        ops.push({
          op: 'set',
          at: c.at,
          key: c.key,
          value: colorFrom(c.space, v[c.slot], v[c.slot + 1], v[c.slot + 2]),
          ...(c.where ? { where: c.where } : {}),
        });
        return;
      case 'enter':
        if (phase === ChannelPhase.done) ops.push(c.op);
        else if (c.height !== undefined && phase === ChannelPhase.moving) {
          const node = c.op.node;
          ops.push({ ...c.op, node: { ...node, params: { ...node.params, height: c.height * v[c.slot] } } });
        }
        return;
      case 'leave':
        if (phase === ChannelPhase.done) ops.push({ op: 'remove', at: c.at });
        else if (c.height !== undefined && phase === ChannelPhase.moving)
          ops.push({ op: 'set', at: c.at, key: 'height', value: c.height * (1 - v[c.slot]) });
        return;
      default:
        if (phase === ChannelPhase.done && c.op) ops.push(c.op);
    }
  });
  return ops.length ? applyLoose(plan.from.root, { ops }, { maxOps: Number.POSITIVE_INFINITY }).root : plan.from.root;
}

/** The view a frame shows: the view (switched at its cut) and each section's viewport at the frame's values. */
export function frameView(plan: TransitionPlan, frame: Frame): ViewState | undefined {
  const base: StoryState['view'] = frameDone(plan, frame) ? plan.to.view : plan.from.view;
  if (!base && !plan.to.view) return undefined;
  let view = plan.from.view?.view;
  const sections: Record<string, StepViewport> = { ...(base?.sections ?? {}) };
  plan.channels.forEach((c, i) => {
    if (c.kind === 'cut' && c.view) {
      if (frame.phase[i] === ChannelPhase.done) view = c.view.to;
      return;
    }
    if (c.kind !== 'view') return;
    const v = frame.values;
    // the base's fields first (one only the start holds is kept until the end), then the channel's
    const out: { -readonly [K in keyof StepViewport]: StepViewport[K] } = {
      ...base?.sections[c.section],
      start: v[c.slot],
      end: v[c.slot + 1],
    };
    // copies of the nested records: the channel's values are written into them, never into the state's
    if (out.bands) out.bands = Object.fromEntries(Object.entries(out.bands).map(([k, b]) => [k, { ...b }]));
    if (out.y) out.y = { ...out.y };
    for (let k = 2; k < c.fields.length; k++) {
      const name = c.fields[k];
      const x = v[c.slot + k];
      if (name === 'top' || name === 'scale') out[name] = x;
      else if (name === 'y.start' || name === 'y.end') {
        const y = (out.y ??= { start: 0, end: 0 }) as { start: number; end: number };
        y[name === 'y.start' ? 'start' : 'end'] = x;
      } else {
        const [, track, field] = /^bands\.(.*)\.(bandTop|rowScale)$/.exec(name) ?? [];
        if (!track) continue;
        const bands = (out.bands ??= {}) as Record<string, { bandTop?: number; rowScale?: number }>;
        (bands[track] ??= {})[field as 'bandTop' | 'rowScale'] = x;
      }
    }
    sections[c.section] = out;
  });
  if (frameDone(plan, frame)) view = plan.to.view?.view;
  return { ...(view !== undefined ? { view } : {}), sections };
}
