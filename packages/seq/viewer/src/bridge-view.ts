/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's bridge chunk (_spec/README.md "Bridges"; across-sections, Moving two sections joined by a bridge): a
// gesture carried through locks, each section moved once (follow, group, lock); a drag or a wheel over a band moving
// both its sides; `bridge.align`, `bridge.selectAcross` and the lock toggle; the keyboard cursor across a band; a
// ribbon's card; the bands' ARIA words. Lazy: `models/bridge-hooks.ts` loads it when a spec holds a bridge.

import { batch } from '@molstar/core/reactive';
import { type Location, location } from '@molstar/seqx-core/query';
import type { BridgeRibbon, ViewLayout } from '@molstar/seqx-core/render';
import { placesOf, type ResolvedBridge, type ResolvedSection, type ResolvedView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { HoverCard } from './hover-card.js';
import type { Place } from './describe.js';
import type { Cursor } from './models/selection.js';
import type { ViewModel, Viewport } from './models/view.js';
import type { SeqViewer } from './viewer.js';

/** A gesture's carrying stops after this many sections (a figure has a handful; a cycle ends earlier). */
const MaxCarried = 64;

/** The bridge at an address, drawn (placed and shown). */
function bridgeAt(view: ResolvedView | undefined, address: string): ResolvedBridge | undefined {
  const all = view?.bridges;
  if (!all) return undefined;
  for (const b of all) if (b.address === address && b.placed && !b.hidden) return b;
  return undefined;
}

/** A bridge's upper side (the section above its band) and its lower. */
const upperOf = (b: ResolvedBridge) => (b.above === 'a' ? b.a : b.b);
const lowerOf = (b: ResolvedBridge) => (b.above === 'a' ? b.b : b.a);

/**
 * What moves a section (_spec/README.md "Bridges", One move a section a gesture): up its attached follows to the
 * first leader (a follower moves with it), then that section's scroll group (its first member).
 */
function moverOf(view: ViewModel, section: string): string {
  const sections = view.resolved.value?.sections;
  let s = section;
  for (let d = 0; d < 16; d++) {
    const r = sections?.find((x) => x.address === s);
    if (!r?.follow || !view.isFollowing(s)) break;
    s = r.follow.leader;
  }
  return view.groupOf(s)[0] ?? s;
}

/**
 * The sections a gesture asked of `from` moves (across-sections, One move a section a gesture): `from` first, then
 * each section a `lock: px` bridge reaches whose mover the closure has not moved yet; a section a follow or a group
 * already moves (a lock then redundant), or one reached again round a cycle, is skipped.
 */
export function carriedSections(view: ViewModel, from: readonly string[]): string[] {
  const moved = new Set<string>();
  const out: string[] = [];
  const take = (s: string) => {
    const m = moverOf(view, s);
    if (moved.has(m)) return false;
    moved.add(m);
    out.push(s);
    return true;
  };
  for (const s of from) take(s);
  const locks = (view.resolved.value?.bridges ?? []).filter((b) => b.lock === 'px' && b.placed);
  for (let grew = locks.length > 0; grew && out.length < MaxCarried; ) {
    grew = false;
    for (const b of locks) {
      const ma = moved.has(moverOf(view, b.a));
      const mb = moved.has(moverOf(view, b.b));
      if (ma !== mb) grew = take(ma ? b.b : b.a) || grew;
    }
  }
  return out;
}

/** A gesture's carrying, kept a resolved view and section: its sections, and the follows it was found under. */
interface Carry {
  readonly to: readonly string[];
  /** Each follower's `isFollowing` when found (a detach or an attach since finds it again). */
  readonly following: readonly boolean[];
}

const Carries = new WeakMap<ResolvedView, Map<string, Carry>>();
/** Each resolved view's followers (their follow state is what a kept carrying is checked against). */
const Followers = new WeakMap<ResolvedView, readonly string[]>();

function followersOf(r: ResolvedView): readonly string[] {
  let f = Followers.get(r);
  if (!f) Followers.set(r, (f = r.sections.filter((s) => s.follow).map((s) => s.address)));
  return f;
}

/**
 * The sections a gesture on `section` moves (a band's: both its sides), kept per resolved view and found again only
 * when a follower was detached or attached since: a drag's moves allocate nothing (_spec/README.md "Bridges").
 */
function carriedOf(
  view: ViewModel,
  r: ResolvedView,
  section: string,
  band: ResolvedBridge | undefined,
): readonly string[] {
  let m = Carries.get(r);
  if (!m) Carries.set(r, (m = new Map()));
  const fs = followersOf(r);
  const c = m.get(section);
  if (c) {
    let same = true;
    for (let i = 0; i < fs.length && same; i++) same = c.following[i] === view.isFollowing(fs[i]);
    if (same) return c.to;
  }
  const to = carriedSections(view, band ? [upperOf(band), lowerOf(band)] : [section]);
  m.set(section, { to, following: fs.map((f) => view.isFollowing(f)) });
  return to;
}

/** A section's (or a band's) viewport: a band's is its upper side's. */
function get(view: ViewModel, section: string): Viewport | undefined {
  const b = bridgeAt(view.resolved.value, section);
  return b ? view.get(upperOf(b)) : undefined;
}

/** The move in progress, read by `runMove` inside one batch (no closure a move). */
const Move = {
  view: undefined as ViewModel | undefined,
  to: [] as readonly string[],
  section: '',
  zoom: false,
  dx: 0,
  dy: 0,
  x: undefined as number | undefined,
  factor: 1,
  anchor: undefined as number | undefined,
  anchorY: undefined as number | undefined,
  origin: undefined as Viewport | undefined,
  moved: false,
};

function runMove(): void {
  const m = Move;
  const view = m.view!;
  for (let i = 0; i < m.to.length; i++) {
    const t = m.to[i];
    const own = t === m.section;
    // a band has no y: its sides pan across only; a carried section takes the gesture's dx alone, and zooms about
    // the same page x
    const moved = m.zoom
      ? view.ownZoom(t, m.factor, own ? m.anchor : sameX(m.origin, view.get(t), m.anchor), own ? m.anchorY : undefined)
      : view.ownPanPx(t, m.dx, own ? m.dy : 0, m.x);
    m.moved = moved || m.moved;
  }
}

function carried(view: ViewModel, section: string): readonly string[] | undefined {
  const r = view.resolved.value;
  if (!r?.bridges) return undefined;
  const band = bridgeAt(r, section);
  const to = carriedOf(view, r, section, band);
  return !band && to.length === 1 ? undefined : to;
}

/** A pan by px carried through the locks (a band's: both its sides, whatever their lock). */
function panPx(view: ViewModel, section: string, dx: number, dy: number, x?: number): boolean {
  const to = carried(view, section);
  if (!to) return view.ownPanPx(section, dx, dy, x);
  const m = Move;
  m.view = view;
  m.to = to;
  m.section = section;
  m.zoom = false;
  m.dx = dx;
  m.dy = dy;
  m.x = x;
  m.moved = false;
  batch(runMove);
  m.view = undefined;
  return m.moved;
}

/** An offset of `v` at the page x of offset `anchor` of `from` (same data x: the sides share the label column). */
function sameX(from: Viewport | undefined, v: Viewport | undefined, anchor: number | undefined): number | undefined {
  if (anchor === undefined || !from || !v || !(from.width > 0) || !(v.width > 0) || from.end <= from.start)
    return undefined;
  const px = ((anchor - from.start) * from.width) / (from.end - from.start);
  return v.start + (px * (v.end - v.start)) / v.width;
}

/** A zoom carried through the locks, about the same page x on each section (a band's: both its sides). */
function zoom(view: ViewModel, section: string, factor: number, anchor?: number, anchorY?: number): boolean {
  const to = carried(view, section);
  if (!to) return view.ownZoom(section, factor, anchor, anchorY);
  const m = Move;
  m.view = view;
  m.to = to;
  m.section = section;
  m.zoom = true;
  m.factor = factor;
  m.anchor = anchor;
  m.anchorY = anchorY;
  m.moved = false;
  // the anchor's page x read before any side moved (a band's: its upper side's viewport)
  m.origin = view.get(section);
  batch(runMove);
  m.view = undefined;
  m.origin = undefined;
  return m.moved;
}

/** The bridge whose band holds a page y (its address): a drag or a wheel there moves both sides. */
function bandAt(layout: ViewLayout, y: number): string | undefined {
  for (const b of layout.bridges?.list ?? []) if (y >= b.box.y && y < b.box.y + b.box.height) return b.address;
  return undefined;
}

/** A ribbon's ends as seen from one side: the near end, the far end. */
interface Ends {
  readonly nearSystem: string;
  readonly near0: number;
  readonly near1: number;
  readonly nearX0: number;
  readonly nearX1: number;
  readonly farSystem: string;
  readonly far0: number;
  readonly far1: number;
}

function endsFrom(r: BridgeRibbon, fromA: boolean): Ends {
  return fromA
    ? {
        nearSystem: r.aSystem,
        near0: r.aStart,
        near1: r.aEnd,
        nearX0: r.ax0,
        nearX1: r.ax1,
        farSystem: r.bSystem,
        far0: r.bStart,
        far1: r.bEnd,
      }
    : {
        nearSystem: r.bSystem,
        near0: r.bStart,
        near1: r.bEnd,
        nearX0: r.bx0,
        nearX1: r.bx1,
        farSystem: r.aSystem,
        far0: r.aStart,
        far1: r.aEnd,
      };
}

/** The drawn ribbons of a bridge in the last frame (topmost last), none before it is laid out. */
function ribbonsOf(layout: ViewLayout | undefined, bridge: string): BridgeRibbon[] {
  return layout?.bridges?.list.find((b) => b.address === bridge)?.ribbons() ?? [];
}

/**
 * The keyboard cursor across a band (across-sections, Hover, selection and the keyboard): onto the section beyond
 * it, at the image of its position through the topmost ribbon covering it (its corresponding position on the far end,
 * a reversed ribbon read backwards); undefined when no ribbon covers it (the view's link, else the window, then).
 */
function cross(
  view: ViewModel,
  layout: ViewLayout | undefined,
  c: Cursor,
  to: string,
): { system: string; position: number } | undefined {
  const b = (view.resolved.value?.bridges ?? []).find(
    (x) => x.placed && !x.hidden && ((x.a === c.section && x.b === to) || (x.b === c.section && x.a === to)),
  );
  if (!b) return undefined;
  const rs = ribbonsOf(layout, b.address);
  for (let i = rs.length - 1; i >= 0; i--) {
    // a stub's far end has no place on the far section: never landed on
    if (rs[i].stub) continue;
    const e = endsFrom(rs[i], b.a === c.section);
    if (e.nearSystem !== c.system || c.position < e.near0 || c.position >= e.near1) continue;
    const t = (c.position - e.near0 + 0.5) / Math.max(1, e.near1 - e.near0);
    const p = rs[i].reversed ? e.far1 - t * (e.far1 - e.far0) : e.far0 + t * (e.far1 - e.far0);
    return { system: e.farSystem, position: Math.min(e.far1 - 1, Math.max(e.far0, Math.floor(p))) };
  }
  return undefined;
}

/** A section's units from position `p0` to `p1` on `system` (its first places), or undefined where it shows neither. */
function unitsOf(section: ResolvedSection, system: string, p0: number, p1: number): [number, number] | undefined {
  const u0 = placesOf(section, system, p0)[0]?.unit;
  const u1 = placesOf(section, system, Math.max(p0, p1 - 1))[0]?.unit;
  if (u0 === undefined || u1 === undefined) return undefined;
  // a position's width in units (a block's scale; a reversed block runs backwards)
  const w = p1 - 1 > p0 ? Math.abs(u1 - u0) / (p1 - 1 - p0) : 1;
  return [Math.min(u0, u1), Math.max(u0, u1) + (w || 1)];
}

/**
 * Where `bridge.align` puts the side across `b` from `near` (across-sections, Moving: align): the hull on the far
 * system holding the most aligned length among the ribbons with their near end in the near side's window; undefined
 * when none is in view.
 */
export function alignRange(
  view: ViewModel,
  layout: ViewLayout | undefined,
  b: ResolvedBridge,
  near: string,
): { far: string; start: number; end: number; location: Location } | undefined {
  const far = b.a === near ? b.b : b.a;
  const box = layout?.bridges?.list.find((x) => x.address === b.address)?.box;
  const section = view.resolved.value?.sections.find((s) => s.address === far);
  if (!box || !section) return undefined;
  const by = new Map<string, { length: number; lo: number; hi: number }>();
  for (const r of ribbonsOf(layout, b.address)) {
    // a stub's far end has no place: it aligns nothing
    if (r.stub) continue;
    const e = endsFrom(r, b.a === near);
    if (e.nearX1 < box.x || e.nearX0 > box.x + box.width) continue;
    const h = by.get(e.farSystem);
    if (h) {
      h.length += r.length;
      h.lo = Math.min(h.lo, e.far0);
      h.hi = Math.max(h.hi, e.far1);
    } else by.set(e.farSystem, { length: r.length, lo: e.far0, hi: e.far1 });
  }
  let best: [string, { length: number; lo: number; hi: number }] | undefined;
  for (const x of by) if (!best || x[1].length > best[1].length) best = x;
  if (!best) return undefined;
  const u = unitsOf(section, best[0], best[1].lo, best[1].hi);
  if (!u) return undefined;
  const loc = location({ parts: [{ system: best[0], ranges: [best[1].lo, best[1].hi] }] });
  return { far, start: u[0], end: u[1], location: loc };
}

/** The bridges touching a section (a middle genome has one above and one below), drawn. */
function bridgesOf(view: ViewModel, section: string | undefined): ResolvedBridge[] {
  return (view.resolved.value?.bridges ?? []).filter(
    (b) => b.placed && !b.hidden && (section === undefined || b.a === section || b.b === section),
  );
}

/** Whether `[s, e)` on `system` meets a Location. */
function meets(loc: Location, system: string, s: number, e: number): boolean {
  for (const p of loc.parts) {
    if (p.system !== system) continue;
    for (let i = 0; i < p.ranges.length; i += 2) if (p.ranges[i] < e && p.ranges[i + 1] > s) return true;
  }
  return false;
}

/** The far ends of the ribbons with an end in `shown` (across `bridges`), as one Location; undefined when none. */
export function acrossOf(
  layout: ViewLayout | undefined,
  bridges: readonly ResolvedBridge[],
  shown: Location,
): Location | undefined {
  const ranges = new Map<string, number[]>();
  const add = (system: string, s: number, e: number) => {
    let r = ranges.get(system);
    if (!r) ranges.set(system, (r = []));
    r.push(s, e);
  };
  for (const b of bridges)
    for (const r of ribbonsOf(layout, b.address)) {
      if (r.stub) continue;
      if (meets(shown, r.aSystem, r.aStart, r.aEnd)) add(r.bSystem, r.bStart, r.bEnd);
      if (meets(shown, r.bSystem, r.bStart, r.bEnd)) add(r.aSystem, r.aStart, r.aEnd);
    }
  if (!ranges.size) return undefined;
  return location({ parts: [...ranges].map(([system, r]) => ({ system, ranges: sortedRuns(r) })) });
}

/** Ranges `[s, e, …]` sorted and merged. */
function sortedRuns(r: number[]): number[] {
  const idx = Array.from({ length: r.length / 2 }, (_, i) => i).sort((x, y) => r[2 * x] - r[2 * y]);
  const out: number[] = [];
  for (const i of idx) {
    const s = r[2 * i];
    const e = r[2 * i + 1];
    if (out.length && s <= out[out.length - 1]) out[out.length - 1] = Math.max(out[out.length - 1], e);
    else out.push(s, e);
  }
  return out;
}

/** A side's name: its section's id. */
const sideName = (address: string) => address.slice(address.lastIndexOf('/') + 1);

/** A system's name (a chromosome's), else its id. */
function systemName(view: ViewModel, system: string): string {
  return view.model.value.systems.get(system)?.name ?? system;
}

const count = (n: number) => n.toLocaleString('en-US');

/** A length of positions as a reader reads it (bp, kb, Mb). */
function lengthText(n: number): string {
  return n >= 1e6 ? `${+(n / 1e6).toFixed(2)} Mb` : n >= 1e3 ? `${+(n / 1e3).toFixed(1)} kb` : `${n} bp`;
}

/**
 * A ribbon's card (across-sections, Hover a ribbon): "peach chr2:1,200,001–1,450,000 ↔ grape chr5:880,001–1,130,000",
 * then reversed and its aligned length; undefined off a ribbon.
 */
function card(view: ViewModel, layout: ViewLayout | undefined, place: Place): HoverCard | undefined {
  const bl = layout?.bridges?.list.find((b) => b.id === place.section);
  const key = place.element?.key;
  // the hit's ribbon alone, never the list (a hover's change)
  const r = key === undefined ? undefined : bl?.ribbon(key);
  if (!bl || !r) return undefined;
  const end = (side: string, system: string, s: number, e: number, unplaced: boolean) =>
    `${sideName(side)} ${systemName(view, system)}:${count(s + 1)}–${count(e)}${unplaced ? ', not shown' : ''}`;
  const title = `${end(bl.a, r.aSystem, r.aStart, r.aEnd, r.stub === 'a')} ↔ ${end(bl.b, r.bSystem, r.bStart, r.bEnd, r.stub === 'b')}`;
  // a bin pair: its count and its bins' width, never an aligned length (it stands for many pairs)
  if (r.bin !== undefined) {
    const n = Math.round(r.count ?? 0);
    const words = `${count(n)} pair${n === 1 ? '' : 's'} · ${lengthText(r.bin)} bins · zoom in to see them`;
    return {
      title,
      rows: [
        { label: 'Pairs', value: count(n) },
        { label: 'Bins', value: lengthText(r.bin) },
      ],
      text: `${title} · ${words}`,
    };
  }
  const rows = [
    { label: 'Strand', value: r.reversed ? 'reversed' : 'same' },
    { label: 'Aligned', value: lengthText(r.length) },
    ...(r.stub ? [{ label: 'Not shown', value: sideName(r.stub === 'a' ? bl.a : bl.b) }] : []),
  ];
  const words = [title, ...(r.reversed ? ['reversed'] : []), lengthText(r.length)];
  return { title, rows, text: words.join(' · ') };
}

/**
 * The bands' ARIA words (across-sections, ARIA): each bridge a group labelled by its `label`, described by the ribbons
 * in view ("1,204 ribbons between peach and grape in view"; none: "no ribbons in view"), read when a frame settles.
 */
export function bridgeAria(view: ViewModel, layout: ViewLayout | undefined): { label: string; description: string }[] {
  return (layout?.bridges?.list ?? []).map((b) => {
    const n = b.count;
    const sides = `between ${sideName(b.a)} and ${sideName(b.b)}`;
    return {
      label: b.label,
      description: n ? `${count(n)} ribbon${n === 1 ? '' : 's'} ${sides} in view` : 'no ribbons in view',
    };
  });
}

/** What the bridge commands read: the section (a side) or the bridge they act on, and an align's way. */
export interface BridgeContext {
  readonly section?: string;
  /** A bridge's address (its menu's). */
  readonly bridge?: string;
  /** An align's near side (the bridge's menu: each way alone); absent: the section's every neighbour. */
  readonly to?: string;
}

/** `Shift+A`, checked against the bindings and every built-in key (`_test/bridge-view.test.ts`). */
export const AlignKey = { key: 'A', shift: true } as const;

/** The bridge commands of a viewer (`bridge.align`, `bridge.selectAcross`, `bridge.lock`). */
export function bridgeCommands(viewer: SeqViewer): CommandDefinition<BridgeContext>[] {
  const view = viewer.view;
  const layout = () => viewer.frame?.layout;
  /**
   * The side a command acts at: the one named (a key's context names the pointer's section), else the keyboard
   * cursor's when the keys moved it, else the leading section: the first in the stack a bridge touches. Never none, so
   * an align is one consistent step (in a chain A–B–C, A's neighbour only: never C to B's old window).
   */
  const sideOf = (ctx: BridgeContext) => {
    const c = viewer.selection.cursor.value;
    const named = ctx.section ?? (c?.keyboard ? c.section : undefined);
    if (named !== undefined) return named;
    const r = view.resolved.value;
    const touched = new Set(bridgesOf(view, undefined).flatMap((b) => [b.a, b.b]));
    return r?.sections.find((s) => touched.has(s.address))?.address;
  };
  const result = (message: string): CommandResult => ({ message });
  return [
    {
      id: 'bridge.align',
      title: 'Align across the bridge',
      category: 'View',
      shortcut: { ...AlignKey },
      when: (ctx) => (bridgesOf(view, ctx.bridge ? undefined : sideOf(ctx)).length ? true : 'No bridge here'),
      run: (_, ctx) => {
        // the bridge's menu: one way; a side: each neighbour across a bridge touching it, one jump
        const pairs: [ResolvedBridge, string][] = [];
        const named = ctx.bridge !== undefined ? bridgeAt(view.resolved.value, ctx.bridge) : undefined;
        if (named) pairs.push([named, ctx.to ?? upperOf(named)]);
        else {
          const side = sideOf(ctx);
          if (side !== undefined) for (const b of bridgesOf(view, side)) pairs.push([b, side]);
        }
        const ranges = pairs.map(([b, near]) => alignRange(view, layout(), b, near)).filter((r) => !!r);
        if (!ranges.length) return result('Nothing aligned in view');
        const moved = view.jump(
          () => ranges.reduce((m, r) => view.setRange(r.far, r.start, r.end) || m, false),
          true,
          'focus',
          ranges[0].location,
        );
        return result(moved === 'none' ? '' : 'Aligned');
      },
    },
    {
      id: 'bridge.selectAcross',
      title: 'Select across the bridge',
      category: 'Selection',
      when: (ctx) => (bridgesOf(view, ctx.bridge ? undefined : sideOf(ctx)).length ? true : 'No bridge here'),
      run: (_, ctx) => {
        const named = ctx.bridge !== undefined ? bridgeAt(view.resolved.value, ctx.bridge) : undefined;
        const bridges = named ? [named] : bridgesOf(view, sideOf(ctx));
        const far = acrossOf(layout(), bridges, viewer.selection.value.value.shown);
        if (!far) return result('Nothing across');
        return result(viewer.selection.add(far) ? 'Selected across' : 'Already selected');
      },
    },
    {
      id: 'bridge.lock',
      title: 'Lock the bridge’s sides',
      category: 'View',
      when: (ctx) => (lockedOf(view, ctx) ? true : 'No bridge here'),
      run: (_, ctx) => {
        const b = lockedOf(view, ctx)!;
        const at = b.address.split('/').filter(Boolean);
        const lock = b.lock !== 'px';
        const o = viewer.spec.patch(lock ? 'Lock a bridge' : 'Unlock a bridge', (_root, revision) => ({
          base: revision,
          ops: [lock ? { op: 'set', at, key: 'lock', value: 'px' } : { op: 'unset', at, key: 'lock' }],
        }));
        const said = (x: { ok: boolean; changed?: boolean }): CommandResult =>
          x.ok ? { message: lock ? 'Locked' : 'Unlocked', undoable: !!x.changed } : { message: 'Refused' };
        return o instanceof Promise ? o.then(said) : said(o);
      },
    },
  ];
}

/** The bridge the lock toggle acts on: the one named, else the first touching the section. */
function lockedOf(view: ViewModel, ctx: BridgeContext): ResolvedBridge | undefined {
  return ctx.bridge !== undefined ? bridgeAt(view.resolved.value, ctx.bridge) : bridgesOf(view, ctx.section)[0];
}

const Installed = new WeakSet<SeqViewer>();

/** Adds the bridge commands to a viewer once (its first view with a bridge). */
function install(viewer: SeqViewer): void {
  if (Installed.has(viewer) || !viewer.view.resolved.value?.bridges) return;
  Installed.add(viewer);
  viewer.commands.push(...bridgeCommands(viewer));
}

const AriaGroups = new WeakMap<HTMLElement, HTMLElement>();

/**
 * The bands' ARIA groups in a mounted viewer (across-sections, ARIA): a visually hidden list beside the canvas, a
 * group a bridge labelled by its `label`, its words the ribbons in view (`bridgeAria`); set when a frame settles,
 * removed when no bridge is drawn.
 */
function aria(viewer: SeqViewer, container: HTMLElement): void {
  const words = bridgeAria(viewer.view, viewer.frame?.layout);
  let host = AriaGroups.get(container);
  if (!words.length) {
    host?.replaceChildren();
    return;
  }
  if (!host) {
    host = container.ownerDocument.createElement('div');
    Object.assign(host.style, {
      position: 'absolute',
      width: '1px',
      height: '1px',
      overflow: 'hidden',
      clip: 'rect(0 0 0 0)',
      whiteSpace: 'nowrap',
      pointerEvents: 'none',
    });
    container.append(host);
    AriaGroups.set(container, host);
  }
  const doc = container.ownerDocument;
  host.replaceChildren(
    ...words.map((w) => {
      const g = doc.createElement('div');
      g.setAttribute('role', 'group');
      g.setAttribute('aria-label', w.label);
      g.textContent = w.description;
      return g;
    }),
  );
}

/** What the eager viewer calls (`models/bridge-hooks.ts`). */
export const BridgeViewParts = { get, panPx, zoom, bandAt, cross, card, install, aria };
