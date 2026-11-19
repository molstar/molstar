/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A transition planned from the diff of two states (_spec/README.md "Stories", Plans): plain data, a list of channels,
// each one property that changes with its kind, its ends and its time window. An application may inspect a plan,
// retime it, drop or add a channel; `sample` reads it at any time.

import { diff, type InsertOp, type MoveOp, type Op } from '@molstar/core/repr/patch';
import { type Address, type Node, resolveAddress, type Value } from '@molstar/core/repr/tree';
import { Color } from '@molstar/core/util/color/color';
import { Oklab } from '@molstar/core/util/color/spaces/oklab';
import { DefaultStoryConfig, type Easing, type StoryConfig } from './config.js';
import type { StepViewport } from './state.js';

/** The view a state shows: the view's address and each section's viewport, in the units the caller reads them. */
export interface ViewState {
  readonly view?: string;
  readonly sections: Readonly<Record<string, StepViewport>>;
}

/** A state a transition goes between: a spec tree and what it shows. */
export interface StoryState {
  readonly root: Node;
  readonly view?: ViewState;
}

interface ChannelTime {
  /** The channel's window within the transition, fractions of it. */
  readonly t0: number;
  readonly t1: number;
  readonly easing: Easing;
  /** Where its values start in a frame's `values`. */
  readonly slot: number;
  /** It needs data not held yet: it holds its start and the frame says it waits. */
  readonly pending?: boolean;
}

/** A section's viewport, field by field (`start`, `end` first; `top`, `scale`, `bands.<track>.bandTop`, `y.start`...). */
export interface ViewChannel extends ChannelTime {
  readonly kind: 'view';
  readonly section: string;
  readonly fields: readonly string[];
  readonly from: readonly number[];
  readonly to: readonly number[];
  /** Which fields move geometrically (a zoom: `scale`, `rowScale`). */
  readonly geometric: readonly boolean[];
  /** The window's path: `pan` (its width geometric, its centre along), or `out-across-in`. */
  readonly path: 'pan' | 'out-across-in';
  /** On a cyclic axis, what the centre's end is moved by so the window goes round the shortest arc. */
  readonly shift: number;
}

/** A numeric param, tweened; an angle (a ring's `rotation`) the short way round, with a param held while it moves. */
export interface NumberChannel extends ChannelTime {
  readonly kind: 'number';
  readonly at: Address;
  readonly key: string;
  readonly where?: 'params' | 'custom';
  readonly from: number;
  readonly to: number;
  /** What it moves by: `to − from`, or the short way round for an angle. */
  readonly delta: number;
  /** A param set while it moves, back as the end state has it at the end (a ring's labels hidden: placed once). */
  readonly hold?: { readonly key: string; readonly value: Value };
}

/** A color param (hex), blended in `space`: three slots. */
export interface ColorChannel extends ChannelTime {
  readonly kind: 'color';
  readonly at: Address;
  readonly key: string;
  readonly where?: 'params' | 'custom';
  readonly from: string;
  readonly to: string;
  readonly space: 'oklab' | 'srgb';
  readonly a: readonly [number, number, number];
  readonly b: readonly [number, number, number];
}

/** A node that appears: inserted at its window's start, its `height` grown from 0 when it has one, else at the cut. */
export interface EnterChannel extends ChannelTime {
  readonly kind: 'enter';
  readonly op: InsertOp;
  readonly height?: number;
  readonly cut: number;
}

/** A node that goes: its `height` shrunk to 0 and removed at its window's end when it has one, else at the cut. */
export interface LeaveChannel extends ChannelTime {
  readonly kind: 'leave';
  readonly at: Address;
  readonly height?: number;
  readonly cut: number;
}

/** Children reordered: each moves to its new place at the cut. */
export interface OrderChannel extends ChannelTime {
  readonly kind: 'order';
  readonly op: MoveOp;
  readonly cut: number;
}

/** Anything else (a layer's type, a data source, a text, the view shown), switched at the cut. */
export interface CutChannel extends ChannelTime {
  readonly kind: 'cut';
  readonly op?: Op;
  /** The view shown switches (a crossfade is drawn as a cut). */
  readonly view?: { readonly from?: string; readonly to?: string };
  readonly cut: number;
}

export type Channel =
  | ViewChannel
  | NumberChannel
  | ColorChannel
  | EnterChannel
  | LeaveChannel
  | OrderChannel
  | CutChannel;
export type ChannelKind = Channel['kind'];

export interface TransitionPlan {
  readonly from: StoryState;
  readonly to: StoryState;
  /** Its natural duration, ms (the step's timing, else the configuration's). */
  readonly duration: number;
  readonly channels: readonly Channel[];
  /** How many values a frame holds. */
  readonly slots: number;
}

export interface PlanOptions {
  readonly duration?: number;
  readonly easing?: Easing;
  readonly config?: StoryConfig;
  /** Every channel a cut at the end. */
  readonly reducedMotion?: boolean;
  /** A cyclic section's axis total, by address: its window goes round the shortest arc. */
  readonly cyclic?: Readonly<Record<string, number>>;
  /** Sections whose viewport needs data not held yet (an axis whose length is pending). */
  readonly pending?: (section: string) => boolean;
  /**
   * A color as drawn, hex, for a value a plan cannot read itself (a token: `token:cat-2`); the viewer's token table.
   * Without it a token is a cut.
   */
  readonly resolveColor?: (color: string) => string | undefined;
}

const isNum = (x: unknown): x is number => typeof x === 'number' && Number.isFinite(x);
const HexColor = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i;

/** A color as hex: itself, or what `resolve` makes of it (a token); undefined when neither reads. */
function hexOf(c: string, resolve: ((c: string) => string | undefined) | undefined): string | undefined {
  if (HexColor.test(c)) return c;
  const r = resolve?.(c);
  return r !== undefined && HexColor.test(r) ? r : undefined;
}

function colorOf(hex: string): Color {
  const h = hex.length === 4 ? `#${hex[1]}${hex[1]}${hex[2]}${hex[2]}${hex[3]}${hex[3]}` : hex;
  return Color.fromHexStyle(h);
}

// pure: a composition with no story part bundles none of this module
const _lab = /* @__PURE__ */ Oklab.zero();
function components(hex: string, space: 'oklab' | 'srgb'): [number, number, number] {
  const c = colorOf(hex);
  if (space === 'srgb') return [(c >> 16) & 255, (c >> 8) & 255, c & 255];
  Oklab.fromColor(_lab, c);
  return [_lab[0], _lab[1], _lab[2]];
}

/** A color from a channel's three values, as hex. */
export function colorFrom(space: 'oklab' | 'srgb', x: number, y: number, z: number): string {
  if (space === 'srgb') {
    const cl = (v: number) => Math.max(0, Math.min(255, Math.round(v)));
    return Color.toHexStyle(Color.fromRgb(cl(x), cl(y), cl(z)));
  }
  _lab[0] = x;
  _lab[1] = y;
  _lab[2] = z;
  return Color.toHexStyle(Oklab.toColor(_lab));
}

const valueAt = (root: Node, at: Address, key: string, where?: 'params' | 'custom'): Value | undefined => {
  const n = resolveAddress(root, at);
  return (where === 'custom' ? n?.custom : n?.params)?.[key];
};

const heightOf = (n: Node | undefined) =>
  isNum(n?.params?.height) && n.params.height > 0 ? n.params.height : undefined;

/** The fields of two viewports both hold, as numbers, with whether each moves geometrically. */
function viewFields(a: StepViewport, b: StepViewport) {
  const fields = ['start', 'end'];
  const from = [a.start, a.end];
  const to = [b.start, b.end];
  const geometric = [false, false];
  const add = (name: string, x: number | undefined, y: number | undefined, geo: boolean) => {
    if (!isNum(x) || !isNum(y)) return;
    fields.push(name);
    from.push(x);
    to.push(y);
    geometric.push(geo && x > 0 && y > 0);
  };
  if (a.top !== undefined || b.top !== undefined) add('top', a.top ?? 0, b.top ?? 0, false);
  add('scale', a.scale, b.scale, true);
  for (const t of Object.keys(b.bands ?? {}).sort()) {
    const x = a.bands?.[t];
    const y = b.bands![t];
    add(`bands.${t}.bandTop`, x?.bandTop, y.bandTop, false);
    add(`bands.${t}.rowScale`, x?.rowScale, y.rowScale, true);
  }
  if (a.y && b.y) {
    add('y.start', a.y.start, b.y.start, false);
    add('y.end', a.y.end, b.y.end, false);
  }
  return { fields, from, to, geometric };
}

/**
 * The transition from one state to another: a channel per property that changes, from their diff and their views.
 * Plain data; `sample` reads it.
 */
export function plan(from: StoryState, to: StoryState, options: PlanOptions = {}): TransitionPlan {
  const config = options.config ?? DefaultStoryConfig;
  const easing = options.easing ?? config.transition.easing;
  const cutAt = config.cut.at;
  // reduced motion: every channel a cut at the end
  const t0 = options.reducedMotion ? 1 : 0;
  const time = { t0, t1: 1, easing };
  const channels: Channel[] = [];
  let slot = 0;
  const take = (n: number) => {
    const s = slot;
    slot += n;
    return s;
  };

  // the view shown, and each section's viewport both states name
  const va = from.view;
  const vb = to.view;
  if (va?.view !== vb?.view)
    channels.push({ kind: 'cut', ...time, slot: take(1), view: { from: va?.view, to: vb?.view }, cut: cutAt });
  if (va && vb && va.view === vb.view)
    for (const [section, b] of Object.entries(vb.sections)) {
      const a = va.sections[section];
      if (!a) continue;
      const f = viewFields(a, b);
      let same = true;
      for (let i = 0; i < f.from.length; i++) if (f.from[i] !== f.to[i]) same = false;
      if (same) continue;
      const total = options.cyclic?.[section];
      let shift = 0;
      if (total && total > 0) {
        const d = (b.start + b.end) / 2 - (a.start + a.end) / 2;
        if (d > total / 2) shift = -total;
        else if (d < -total / 2) shift = total;
      }
      const s0 = a.start;
      const e0 = a.end;
      const s1 = b.start + shift;
      const e1 = b.end + shift;
      // windows that touch pan too: only a gap between them goes out and back in
      const overlap = Math.max(s0, s1) <= Math.min(e0, e1);
      const path =
        config.view.path === 'auto'
          ? overlap
            ? 'pan'
            : 'out-across-in'
          : (config.view.path as 'pan' | 'out-across-in');
      channels.push({
        kind: 'view',
        ...time,
        slot: take(f.fields.length),
        section,
        ...f,
        path,
        shift,
        ...(options.pending?.(section) ? { pending: true } : {}),
      });
    }

  // the spec: a channel per op of the diff, in its order (so a frame's ops apply in it)
  const ops = diff(stripStory(from.root), stripStory(to.root)).ops;
  for (const op of ops) {
    switch (op.op) {
      case 'set': {
        const where = op.where === 'custom' ? 'custom' : undefined;
        const a = valueAt(from.root, op.at, op.key, where);
        const b = op.value;
        if (isNum(a) && isNum(b)) {
          const n = resolveAddress(to.root, op.at);
          // a ring's rotation: the short way round, its labels hidden while it turns (placed once, at its end)
          const angle = op.key === 'rotation' && n?.kind === 'section';
          let delta = b - a;
          if (angle) delta = ((((delta + 180) % 360) + 360) % 360) - 180;
          channels.push({
            kind: 'number',
            ...time,
            slot: take(1),
            at: op.at,
            key: op.key,
            ...(where ? { where } : {}),
            from: a,
            to: b,
            delta,
            ...(angle ? { hold: { key: 'labels', value: 'none' } } : {}),
          });
          break;
        }
        // a theme's resolved colors blend: a token through the color it resolves to
        const ha = typeof a === 'string' ? hexOf(a, options.resolveColor) : undefined;
        const hb = typeof b === 'string' ? hexOf(b, options.resolveColor) : undefined;
        if (typeof a === 'string' && typeof b === 'string' && ha && hb) {
          const space = config.color.space;
          channels.push({
            kind: 'color',
            ...time,
            slot: take(3),
            at: op.at,
            key: op.key,
            ...(where ? { where } : {}),
            from: a,
            to: b,
            space,
            a: components(ha, space),
            b: components(hb, space),
          });
          break;
        }
        channels.push({ kind: 'cut', ...time, slot: take(1), op, cut: cutAt });
        break;
      }
      case 'insert': {
        const height = heightOf(op.node);
        channels.push({ kind: 'enter', ...time, slot: take(1), op, ...(height ? { height } : {}), cut: cutAt });
        break;
      }
      case 'remove': {
        const height = heightOf(resolveAddress(from.root, op.at));
        channels.push({ kind: 'leave', ...time, slot: take(1), at: op.at, ...(height ? { height } : {}), cut: cutAt });
        break;
      }
      case 'move':
        if (!op.rename && sameAddress(op.parent, op.at.slice(0, -1))) {
          channels.push({ kind: 'order', ...time, slot: take(1), op, cut: cutAt });
          break;
        }
        channels.push({ kind: 'cut', ...time, slot: take(1), op, cut: cutAt });
        break;
      default:
        channels.push({ kind: 'cut', ...time, slot: take(1), op, cut: cutAt });
    }
  }
  // a cut takes no time: every channel switches at once
  const duration = options.reducedMotion ? 0 : (options.duration ?? config.transition.duration);
  return { from, to, duration, channels, slots: slot };
}

function sameAddress(a: Address, b: Address): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** A spec without its story (the steps are not what a transition shows). */
export function stripStory(root: Node): Node {
  const c = root.children;
  if (!c?.some((n) => n.kind === 'snapshots')) return root;
  return { ...root, children: c.filter((n) => n.kind !== 'snapshots') };
}
