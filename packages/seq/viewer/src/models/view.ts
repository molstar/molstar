/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The view's model (_spec/README.md "ViewModel"): the view shown, a viewport per section, scroll groups on equal axes
// (EC-113), follows through the resolved view's follow path (EC-113, EC-114), a move asked of a follower carried back
// to its leader (`follow.followerMoves`: `lead`), viewports kept across updates (EC-111) and moves before a width kept
// as intents (EC-112). Session state: never in the spec or its history.

import { batch, derive, KeyedState, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import { type Node, parseAddress } from '@molstar/core/repr/tree';
import { inversePath, type Path } from '@molstar/seqx-core/resolve';
import {
  type Axis,
  type Follow,
  followRange,
  type ResolvedSection,
  type ResolvedView,
  resolveView,
  resolveViews,
} from '@molstar/seqx-core/view';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { axesEqual, axisKey, axisTotal, nearestOffsetOf, offsetOf, positionAt } from './axis.js';
import type { SpecModel } from './spec.js';

/**
 * A section's viewport: `[start, end)` in offsets along its axis (fractional; positions on an axis that shows its
 * whole system from 0), its vertical scroll, its width in pixels. `pending` while it cannot be placed (no width yet,
 * or the axis waits for its length): the range is then the intent, applied when it can be.
 */
export interface Viewport {
  readonly start: number;
  /** Past the axis's total only when the range wraps a cyclic origin (`pieces`). */
  readonly end: number;
  readonly top: number;
  readonly width: number;
  readonly pending: boolean;
  /** The section it follows, when it does (a move asked of it moves the leader, `follow.followerMoves`). */
  readonly leader?: string;
  /** A follower's range over a cyclic origin, as two offset pieces in display order: `[start, total, 0, end − total]`. */
  readonly pieces?: readonly number[];
}

interface Group {
  readonly key: string;
  readonly axis: Axis;
  start: number;
  end: number;
  /** It shows the whole axis and stays so across resizes and updates. */
  fitted: boolean;
  /** The last width (px) it was placed at: a resize keeps its zoom against it. */
  placedWidth: number;
  readonly members: string[];
}

interface SectionState {
  readonly section: ResolvedSection;
  /** Its group's key; undefined for an attached follower. */
  group?: string;
  width: number;
  top: number;
}

const DetachedPrefix = '#';

export class ViewModel extends ReactiveModel {
  private readonly _view = new State<string | undefined>(undefined);
  private readonly _viewports = new KeyedState<string, Viewport>();
  private _sections = new Map<string, SectionState>();
  private _groups = new Map<string, Group>();
  /** Sections in an order where each leader comes before its followers. */
  private _order: string[] = [];
  private readonly _spec: SpecModel;
  readonly config: ViewerConfig;
  /** The address of the view shown; undefined: the spec's first view. */
  readonly view: ReadonlyState<string | undefined> = this._view.asReadonly;
  /** The resolved view shown, undefined when the spec has none. */
  readonly resolved: ReadonlyState<ResolvedView | undefined>;

  /** The lengths of source-backed systems the store's coverages give, by system address (`refreshLengths`). */
  private readonly _lengths = new State<ReadonlyMap<string, number>>(new Map(), sameLengths);
  private _lengthOf: ((system: string) => number | undefined) | undefined;

  constructor(spec: SpecModel, config: ViewerConfig = DefaultViewerConfig) {
    super();
    this._spec = spec;
    this.config = config;
    this.resolved = derive([spec.model, this._view, this._lengths], (model, address, lengths) =>
      withLengths(model.spec, address === undefined ? resolveViews(model)[0] : resolveView(model, address), lengths),
    );
    this._update(this.resolved.value);
    this.listenPermanently(this.resolved, (v) => this._update(v));
  }

  /**
   * Where the length of a system the spec does not know comes from (the store's coverage of its source); read by
   * `refreshLengths`, which the viewer calls when the store changes.
   */
  setLengthSource(lengthOf: (system: string) => number | undefined): void {
    this._lengthOf = lengthOf;
    this.refreshLengths();
  }

  /**
   * Reads the lengths of the sections' systems whose axis waits for one: a section on a source-backed system is placed
   * once its source's coverage gives the length (the axis gets its pieces, a kept intent applies, EC-112).
   */
  refreshLengths(): void {
    const lengthOf = this._lengthOf;
    if (!lengthOf) return;
    const model = this._spec.model.value;
    const view = this._view.value === undefined ? resolveViews(model)[0] : resolveView(model, this._view.value);
    const next = new Map<string, number>();
    for (const sec of view?.sections ?? []) {
      for (const axis of sec.follow ? [sec.axis, sec.follow.leaderAxis] : [sec.axis]) {
        if (axis.ranges || next.has(axis.system)) continue;
        const n = lengthOf(axis.system);
        if (n !== undefined && Number.isInteger(n) && n >= 0) next.set(axis.system, n);
      }
    }
    this._lengths.set(next);
  }

  /** Shows another view of the spec (undefined: the first). */
  show(view: string | undefined): void {
    this._view.set(view);
  }

  /** A section's viewport as a state: only its readers are told of its changes. */
  viewport(section: string): ReadonlyState<Viewport | undefined> {
    return this._viewports.key(section);
  }

  get(section: string): Viewport | undefined {
    return this._viewports.get(section);
  }

  /** The sections sharing a section's viewport (itself included); a follower's is its own. */
  groupOf(section: string): readonly string[] {
    const s = this._sections.get(section);
    return s?.group !== undefined ? (this._groups.get(s.group)?.members ?? [section]) : [section];
  }

  /** The section's pixel width (0: hidden for lack of width; its viewport is kept, EC-112). */
  setWidth(section: string, width: number): void {
    const s = this._sections.get(section);
    if (!s || s.width === width || !(width >= 0)) return;
    s.width = width;
    const g = this._groupOf(s);
    if (g) {
      const w = this._groupWidth(g);
      if (
        w > 0 &&
        g.placedWidth > 0 &&
        w !== g.placedWidth &&
        !g.fitted &&
        this.config.viewport.onResize === 'keep-zoom'
      ) {
        const c = (g.start + g.end) / 2;
        const half = ((g.end - g.start) * w) / g.placedWidth / 2;
        g.start = c - half;
        g.end = c + half;
      }
    }
    this._publish();
  }

  setTop(section: string, top: number): void {
    const s = this._sections.get(section);
    if (!s || s.top === top || !Number.isFinite(top)) return;
    s.top = top;
    this._publish();
  }

  /** Shows `[start, end)` (offsets); before a width, the intent (EC-112). False when the section cannot move. */
  setRange(section: string, start: number, end: number): boolean {
    if (!Number.isFinite(start) || !Number.isFinite(end)) return false;
    const lead = this._leading(section);
    if (lead) {
      const v = this.get(section);
      const [s, e] = [Math.min(start, end), Math.max(start, end)];
      return !!v && this._lead(lead, s, e, Math.sign(s - v.start), Math.sign(e - v.end));
    }
    const g = this._movable(section);
    if (!g) return false;
    g.start = Math.min(start, end);
    g.end = Math.max(start, end);
    g.fitted = false;
    this._publish();
    return true;
  }

  /** Moves the range by `delta` offsets. */
  pan(section: string, delta: number): boolean {
    if (!Number.isFinite(delta)) return false;
    if (delta === 0) return this._canMove(section);
    const lead = this._leading(section);
    if (lead) {
      const v = this.get(section);
      const m = Math.sign(delta);
      return !!v && this._lead(lead, v.start + delta, v.end + delta, m, m, delta);
    }
    const g = this._movable(section);
    if (!g) return false;
    g.start += delta;
    g.end += delta;
    g.fitted = false;
    this._publish();
    return true;
  }

  /** Moves the range by `dx` pixels of its width; nothing before a width. */
  panPx(section: string, dx: number): boolean {
    const v = this.get(section);
    if (!v || v.width <= 0) return false;
    return this.pan(section, (dx * (v.end - v.start)) / v.width);
  }

  /** Zooms by `factor` (above 1: in) about `anchor` (an offset; the center when absent), which stays in place. */
  zoom(section: string, factor: number, anchor?: number): boolean {
    if (!(factor > 0) || !Number.isFinite(factor) || (anchor !== undefined && !Number.isFinite(anchor))) return false;
    const lead = this._leading(section);
    if (lead) {
      const v = this.get(section);
      if (!v) return false;
      const a = anchor ?? (v.start + v.end) / 2;
      const s = a - (a - v.start) / factor;
      const e = a + (v.end - a) / factor;
      return this._lead(lead, s, e, Math.sign(s - v.start), Math.sign(e - v.end));
    }
    const g = this._movable(section);
    if (!g) return false;
    const a = anchor ?? (g.start + g.end) / 2;
    g.start = a - (a - g.start) / factor;
    g.end = a + (g.end - a) / factor;
    g.fitted = false;
    this._publish();
    return true;
  }

  /** Shows the whole axis; a detached follower is attached again. */
  fit(section: string): boolean {
    const s = this._sections.get(section);
    if (!s) return false;
    if (s.section.follow) return this.attach(section) || this.isFollowing(section);
    const g = this._groupOf(s);
    if (!g) return false;
    g.start = 0;
    g.end = axisTotal(g.axis) ?? 0;
    g.fitted = true;
    this._publish();
    return true;
  }

  /** Attaches a detached follower to its leader again. */
  attach(section: string): boolean {
    const s = this._sections.get(section);
    if (!s?.section.follow || s.group === undefined) return false;
    this._groups.delete(s.group);
    s.group = undefined;
    this._publish();
    return true;
  }

  /** Whether the section follows another (and is attached to it). */
  isFollowing(section: string): boolean {
    const s = this._sections.get(section);
    return !!s?.section.follow && s.group === undefined;
  }

  /**
   * Whether a move asked of the section would move something: a follower's when it leads (its leader can move), or when
   * it may detach.
   */
  canMove(section: string): boolean {
    return this._canMove(section);
  }

  // ---- internals ----

  private _canMove(section: string, depth = 0): boolean {
    const s = this._sections.get(section);
    if (!s) return false;
    if (s.section.follow && s.group === undefined) {
      const moves = this.config.follow.followerMoves;
      if (moves === 'lead') return depth < MaxFollowDepth && this._canMove(s.section.follow.leader, depth + 1);
      return moves === 'detach';
    }
    return !!this._groupOf(s);
  }

  /**
   * An attached follower whose moves lead (`follow.followerMoves`: `lead`) to a leader that can move, else undefined.
   * `_canMove` bounds the chain (`MaxFollowDepth`), so the move `_lead` passes up it ends.
   */
  private _leading(section: string): SectionState | undefined {
    const s = this._sections.get(section);
    if (!s?.section.follow || s.group !== undefined || this.config.follow.followerMoves !== 'lead') return undefined;
    return this._canMove(section) ? s : undefined;
  }

  /**
   * A move asked of a follower moves its leader (_spec/README.md "ViewModel", Follow): the follower's asked `[start,
   * end)` (offsets) is carried back through the follow, each end to the leader place whose image it is, so the
   * follower shows what was asked and every section moves in step; on a minus strand the follower's ends are the
   * images of the leader's other ends, so its pan moves the leader the other way. `ms`, `me`: which way each end
   * moved. An end that maps nowhere (an intron, a UTR) takes the nearest place that maps the way it moved (`backTo`).
   * A pan (`delta`) moves the follower's unclamped image of the leader, not the window its own `maxPxPerPosition`
   * widened (which would zoom the leader out), and when one end maps nowhere the leader moves by the other's step: the
   * leader keeps its width. A zoom or a range whose end lies past every image outward reaches the leader's axis edge
   * there (a zoom out never stalls); inward, it keeps the leader's end. A leader that follows passes the move on up the
   * chain, bounded by `_leading`.
   */
  private _lead(s: SectionState, start: number, end: number, ms: number, me: number, delta?: number): boolean {
    const follow = s.section.follow!;
    const lead = this.get(follow.leader);
    const v = this.get(s.section.address);
    const axis = s.section.axis;
    const total = axisTotal(axis);
    if (!lead || !v || total === undefined) return false;
    const la = follow.leaderAxis;
    const ps = positionAt(la, lead.start);
    const pe = positionAt(la, lead.end);
    if (ps === undefined || pe === undefined) return false;
    if (delta !== undefined) {
      // the image of the leader's window before the follower's clamp, moved by the pan
      const image = carried(s.section, follow, ps, pe, this.config);
      const a = image && nearestOffsetOf(axis, image[0]);
      const b = image && nearestOffsetOf(axis, image[1]);
      if (a !== undefined && b !== undefined && a < b) [start, end] = [a + delta, b + delta];
    }
    // a cyclic follower's range over its origin runs past the total
    const pos = (o: number) => positionAt(axis, axis.cyclic && o > total ? o - total : o);
    const lo = pos(start);
    const hi = pos(end);
    const inv = inverseOf(follow.path);
    const xLo = lo === undefined ? undefined : backTo(follow, inv, axis, lo, ms);
    const xHi = hi === undefined ? undefined : backTo(follow, inv, axis, hi, me);
    // the follower's start is the image of the leader's start, or of its end on a minus strand
    let a = follow.dir === 1 ? xLo : xHi;
    let b = follow.dir === 1 ? xHi : xLo;
    if (delta !== undefined) {
      if (a === undefined && b === undefined) return false;
      if (a === undefined) a = ps + (b! - pe);
      if (b === undefined) b = pe + (a - ps);
    } else {
      // an end moved outward past every image: the leader's axis edge on that side (its clamp keeps it in)
      const startOut = follow.dir === 1 ? ms < 0 : me > 0;
      const endOut = follow.dir === 1 ? me > 0 : ms < 0;
      const lt = axisTotal(la);
      a ??= startOut && lt !== undefined ? (positionAt(la, 0) ?? ps) : ps;
      b ??= endOut && lt !== undefined ? (positionAt(la, lt) ?? pe) : pe;
    }
    const oa = nearestOffsetOf(la, a);
    const ob = nearestOffsetOf(la, b);
    if (oa === undefined || ob === undefined || oa === ob) return false;
    return this.setRange(follow.leader, Math.min(oa, ob), Math.max(oa, ob));
  }

  private _groupOf(s: SectionState): Group | undefined {
    return s.group === undefined ? undefined : this._groups.get(s.group);
  }

  private _groupWidth(g: Group): number {
    let w = 0;
    for (const m of g.members) w = Math.max(w, this._sections.get(m)?.width ?? 0);
    return w;
  }

  /**
   * The group a move on the section moves, after the move's arguments were checked: none for an attached follower
   * that ignores moves (EC-114); a detached one's own.
   */
  private _movable(section: string): Group | undefined {
    if (!this._canMove(section)) return undefined;
    const s = this._sections.get(section)!;
    if (s.section.follow && s.group === undefined) {
      // detached: its own viewport, from where it is now
      const v = this.get(section);
      const key = `${DetachedPrefix}${section}`;
      const total = axisTotal(s.section.axis) ?? 0;
      this._groups.set(key, {
        key,
        axis: s.section.axis,
        start: v?.start ?? 0,
        end: v?.end ?? total,
        fitted: false,
        placedWidth: s.width,
        members: [section],
      });
      s.group = key;
    }
    return this._groupOf(s);
  }

  /** A new resolved view: sections and groups again, each viewport kept by the rules of EC-111. */
  private _update(view: ResolvedView | undefined): void {
    const root = this._spec.model.value.spec;
    const oldSections = this._sections;
    const oldGroups = this._groups;
    const sections = new Map<string, SectionState>();
    const groups = new Map<string, Group>();
    for (const section of view?.sections ?? []) {
      const old = oldSections.get(section.address);
      const state: SectionState = { section, width: old?.width ?? 0, top: old?.top ?? 0 };
      sections.set(section.address, state);
      if (section.follow) {
        // a detached follower stays detached while its axis is the same
        const kept = old?.group?.startsWith(DetachedPrefix) ? oldGroups.get(old.group) : undefined;
        if (kept && axesEqual(kept.axis, section.axis)) {
          groups.set(kept.key, { ...kept, axis: section.axis });
          state.group = kept.key;
        }
        continue;
      }
      // a pending axis has no pieces: its key holds the section's own range, so only equal ranges group
      const key = this.config.viewport.groupEqualAxes
        ? section.axis.ranges
          ? axisKey(section.axis)
          : `${axisKey(section.axis)}|${rangeParam(root, section.address)}`
        : section.address;
      state.group = key;
      const existing = groups.get(key);
      if (existing) {
        existing.members.push(section.address);
        continue;
      }
      const same = oldGroups.get(key);
      const before = old?.group !== undefined ? oldGroups.get(old.group) : undefined;
      groups.set(key, carry(key, section, same && axesEqual(same.axis, section.axis) ? same : before));
    }
    batch(() => {
      this._sections = sections;
      this._groups = groups;
      this._order = followOrder(view);
      for (const address of [...this._viewports.values.keys()])
        if (!sections.has(address)) this._viewports.delete(address);
      this._publish();
    });
  }

  /** Clamps each group where it can be placed, then writes every section's viewport, leaders before followers. */
  private _publish(): void {
    batch(() => {
      const { maxPxPerPosition } = this.config.viewport;
      for (const g of this._groups.values()) {
        const total = axisTotal(g.axis);
        const width = this._groupWidth(g);
        if (total === undefined || width <= 0) continue;
        g.placedWidth = width;
        if (g.fitted) {
          g.start = 0;
          g.end = total;
          continue;
        }
        const w = Math.min(Math.max(g.end - g.start, Math.min(width / maxPxPerPosition, total)), total);
        const c = (g.start + g.end) / 2;
        g.start = Math.min(Math.max(c - w / 2, 0), total - w);
        g.end = g.start + w;
        if (g.start === 0 && g.end === total) g.fitted = true;
      }
      for (const address of this._order) {
        const s = this._sections.get(address)!;
        const g = this._groupOf(s);
        if (g) {
          const pending = axisTotal(s.section.axis) === undefined || this._groupWidth(g) <= 0;
          this._set(address, g.start, g.end, s, pending, undefined, undefined);
        } else this._follow(s);
      }
    });
  }

  /** Writes a viewport, comparing first: an equal one makes nothing and tells nothing. */
  private _set(
    address: string,
    start: number,
    end: number,
    s: SectionState,
    pending: boolean,
    leader: string | undefined,
    pieces: readonly number[] | undefined,
  ): void {
    const p = this._viewports.get(address);
    if (
      p &&
      p.start === start &&
      p.end === end &&
      p.top === s.top &&
      p.width === s.width &&
      p.pending === pending &&
      p.leader === leader &&
      piecesEqual(p.pieces, pieces)
    )
      return;
    this._viewports.set(address, {
      start,
      end,
      top: s.top,
      width: s.width,
      pending,
      ...(leader !== undefined ? { leader } : {}),
      ...(pieces ? { pieces } : {}),
    });
  }

  /**
   * A follower's viewport: its leader's range carried along the follow path at its scale (1, 3 or 1/3). Each end is
   * the continuous image of the leader's (`point` at the start, the far edge of the last element at a whole end), so
   * it moves smoothly within a run; where an end maps nowhere (an intron, a UTR), the hull of the whole positions'
   * image (`followRange`). Over a cyclic follower's origin the range wraps (`pieces`), or with `follow.cyclicArc:
   * 'hull'` is the hull. Clamped by `maxPxPerPosition` like a leader. A follower whose leader has no viewport yet keeps
   * its own until the leader's arrives.
   */
  private _follow(s: SectionState): void {
    const follow = s.section.follow!;
    const lead = this.get(follow.leader);
    const address = s.section.address;
    const prev = this.get(address);
    const axis = s.section.axis;
    const total = axisTotal(axis);
    if (!lead) {
      if (prev) this._set(address, prev.start, prev.end, s, true, follow.leader, prev.pieces);
      return;
    }
    const pending = lead.pending || total === undefined || s.width <= 0;
    let start = prev?.start ?? 0;
    let end = prev?.end ?? total ?? 0;
    let pieces: readonly number[] | undefined = prev?.pieces;
    const ps = positionAt(follow.leaderAxis, lead.start);
    const pe = positionAt(follow.leaderAxis, lead.end);
    const range = ps === undefined || pe === undefined ? undefined : carried(s.section, follow, ps, pe, this.config);
    if (range && total !== undefined) {
      const [lo, hi] = range;
      const a = nearestOffsetOf(axis, lo);
      const b = nearestOffsetOf(axis, hi);
      if (a !== undefined && b !== undefined) {
        if (a > b && axis.cyclic) {
          // over the origin: on to the total, then from 0
          start = a;
          end = total + b;
          pieces = [a, total, 0, b];
        } else {
          [start, end] = [Math.min(a, b), Math.max(a, b)];
          pieces = undefined;
          if (s.width > 0) [start, end] = clampRange(start, end, total, s.width, this.config.viewport.maxPxPerPosition);
        }
      }
    }
    this._set(address, start, end, s, pending, follow.leader, pieces);
  }
}

/**
 * A group's viewport across an update (EC-111): an equal axis keeps it; a changed one keeps its center and zoom when
 * the center is still shown, else fits; an axis that was pending (no length) keeps its intent where offsets are
 * positions from 0 on the same system. A new group fits.
 */
function carry(key: string, section: ResolvedSection, prev: Group | undefined): Group {
  const axis = section.axis;
  const total = axisTotal(axis) ?? 0;
  const group: Group = { key, axis, start: 0, end: total, fitted: true, placedWidth: 0, members: [section.address] };
  if (!prev) return group;
  group.placedWidth = prev.placedWidth;
  if (prev.fitted) return group;
  if (axesEqual(prev.axis, axis)) return { ...group, start: prev.start, end: prev.end, fitted: false };
  if (prev.axis.system !== axis.system) return group;
  if (!prev.axis.ranges) {
    // a pending axis's offsets are positions from 0: kept where the new axis starts there
    return axis.ranges?.[0] === 0 ? { ...group, start: prev.start, end: prev.end, fitted: false } : group;
  }
  const center = positionAt(prev.axis, (prev.start + prev.end) / 2);
  const c = center === undefined ? undefined : offsetOf(axis, center);
  if (c === undefined) return group;
  const half = (prev.end - prev.start) / 2;
  return { ...group, start: c - half, end: c + half, fitted: false };
}

/** A leader's `[ps, pe)` positions carried to the follower's positions, `[lo, hi]` in its display direction. */
function carried(
  section: ResolvedSection,
  follow: Follow,
  ps: number,
  pe: number,
  config: ViewerConfig,
): [number, number] | undefined {
  let a = follow.point(ps);
  let b = pathEnd(follow.path, pe);
  if (a === undefined || b === undefined) {
    const arc = config.follow.cyclicArc === 'hull' ? 'hull' : 'shortest';
    const hull = followRange(section, [Math.floor(ps), Math.ceil(pe)], { arc });
    if (!hull || hull.length < 2) return undefined;
    const lo = hull[0];
    const hi = hull[hull.length - 1];
    if (a === undefined) a = follow.dir === 1 ? lo : hi;
    if (b === undefined) b = follow.dir === 1 ? hi : lo;
  }
  const [lo, hi] = follow.dir === 1 ? [a, b] : [b, a];
  // ends past each other: over a cyclic follower's origin, the wrap (or its hull when asked); on a linear follower,
  // the leader's range wraps the coding sequence's origin, so the hull of the image
  if (lo > hi || (lo === hi && pe > ps)) {
    if (section.axis.cyclic && config.follow.cyclicArc !== 'hull') return [lo, hi];
    const hull = followRange(section, [Math.floor(ps), Math.ceil(pe)], { arc: 'hull' });
    return hull && hull.length >= 2 ? [hull[0], hull[hull.length - 1]] : undefined;
  }
  return [lo, hi];
}

/** How many follows a move asked of a follower goes up at most (a chain of follows; a cycle is cut). */
const MaxFollowDepth = 16;

/**
 * The continuous image of a position along a path: where element `floor(x)` starts plus the fraction in its run's
 * direction (times a codon step's scale); undefined where the element maps nowhere. `Follow.point` for any path, the
 * follow's inverse too.
 */
function pathPoint(p: Path, x: number): number | undefined {
  const k = Math.floor(x);
  const frac = x - k;
  const fine = p.step === 1 ? 3 * k : k;
  const r = p.linear.forward(fine, fine + 1)[0];
  if (!r) return undefined;
  const edge = r.dir === 1 ? r.to : r.to + 1;
  if (p.step === 1) return edge + frac * 3 * r.dir;
  if (p.step === -1) return (edge + frac * r.dir) / 3;
  return edge + frac * r.dir;
}

/** The image of a range's end: the far edge of its last element's last place, or `pathPoint` within an element. */
function pathEnd(p: Path, x: number): number | undefined {
  if (!Number.isInteger(x)) return pathPoint(p, x);
  const k = x - 1;
  const fine = p.step === 1 ? 3 * k + 2 : k;
  const r = p.linear.forward(fine, fine + 1)[0];
  if (!r) return undefined;
  const edge = r.dir === 1 ? r.to + 1 : r.to;
  return p.step === -1 ? edge / 3 : edge;
}

/** A follow's path read back, made once per path. */
const inverses = new WeakMap<Path, Path>();
function inverseOf(p: Path): Path {
  let inv = inverses.get(p);
  if (!inv) inverses.set(p, (inv = inversePath(p)));
  return inv;
}

/** Just short of an element's far edge: the limit of its continuous image there, not the next element's start. */
const BeforeEdge = 1 - 1e-9;

/**
 * The leader position whose image (`Follow.point`, what the follower's ends are made of) is the follower's position
 * `y`: the inverse of the forward map itself, so a follower moved through it shows what was asked. Each leader
 * element's image is a run of `point(k) + f · slope`, `f` in [0, 1) (a codon split by an intron runs on from its first
 * nucleotide: its image is not its nucleotides); the elements tried are those near `y`'s preimage through `inv` (the
 * follow read back) and near the nearest mapped positions on either side. Where `y` is in no element's image (an
 * intron, a UTR), the nearest image the way the end moved (`m`), so a pan steps over a gap and never sticks. Undefined
 * when nothing lies that way.
 */
function backTo(follow: Follow, inv: Path, axis: Axis, y: number, m: number): number | undefined {
  const tried = new Set<number>();
  let best: number | undefined;
  let distance = Infinity;
  const near = (x: number | undefined) => {
    if (x === undefined) return undefined;
    const k0 = Math.floor(x);
    for (let k = k0 - 1; k <= k0 + 1; k++) {
      if (tried.has(k)) continue;
      tried.add(k);
      const a = follow.point(k);
      const half = follow.point(k + 0.5);
      if (a === undefined || half === undefined || half === a) continue;
      const slope = 2 * (half - a);
      const f = (y - a) / slope;
      if (f >= 0 && f < 1) return k + f;
      if (m === 0) continue;
      // the end of this element's image nearer `y`, when it lies the way the end moved
      const b = a + slope;
      const v = m > 0 ? Math.min(a, b) : Math.max(a, b);
      const d = (v - y) * m;
      if (d > 0 && d < distance) {
        distance = d;
        best = k + (v === a ? 0 : BeforeEdge);
      }
    }
    return undefined;
  };
  return (
    near(pathPoint(inv, y)) ??
    near(nearestMapped(inv, axis, y, 1)) ??
    near(nearestMapped(inv, axis, y, -1)) ??
    (m === 0 ? undefined : best)
  );
}

/**
 * The leader place of the follower's nearest mapped position past `y` the way `m` says, within the follower's axis:
 * where the next one starts (`m` > 0), or the far edge of the last one before (`m` < 0).
 */
function nearestMapped(inv: Path, axis: Axis, y: number, m: number): number | undefined {
  const r = axis.ranges;
  if (!r || r.length < 2) return undefined;
  let lo = Infinity;
  let hi = -Infinity;
  for (let i = 0; i + 1 < r.length; i += 2) {
    lo = Math.min(lo, r[i]);
    hi = Math.max(hi, r[i + 1]);
  }
  // the inverse's source is in thirds of a residue when it steps up to nucleotides
  const scale = inv.step === 1 ? 3 : 1;
  if (m > 0) {
    const from = Math.ceil(y);
    if (from >= hi) return undefined;
    const run = inv.linear.forward(from * scale, hi * scale)[0];
    return run ? pathPoint(inv, Math.ceil(run.from / scale)) : undefined;
  }
  const to = Math.floor(y);
  if (to <= lo) return undefined;
  const runs = inv.linear.forward(lo * scale, to * scale);
  const last = runs[runs.length - 1];
  return last ? pathEnd(inv, Math.floor((last.from + last.length - 1) / scale) + 1) : undefined;
}

/** At least `width / maxPx` wide (within the total), centered, within the axis. */
function clampRange(start: number, end: number, total: number, width: number, maxPx: number): [number, number] {
  const w = Math.min(Math.max(end - start, Math.min(width / maxPx, total)), total);
  if (w === end - start) return [start, end];
  const c = (start + end) / 2;
  const s = Math.min(Math.max(c - w / 2, 0), total - w);
  return [s, s + w];
}

function piecesEqual(a: readonly number[] | undefined, b: readonly number[] | undefined): boolean {
  if (a === b) return true;
  if (!a || !b || a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** A section's own `range` param as written ('all' without one): part of a pending axis's group key. */
function rangeParam(root: Node, address: string): string {
  const p = parseAddress(address).pattern;
  let n: Node | undefined = root;
  for (const seg of p?.segments ?? []) n = typeof seg === 'string' ? n?.children?.find((c) => c.id === seg) : undefined;
  const r = n?.params?.range;
  return Array.isArray(r) ? r.join(',') : 'all';
}

function sameLengths(a: ReadonlyMap<string, number>, b: ReadonlyMap<string, number>): boolean {
  if (a.size !== b.size) return false;
  for (const [k, v] of a) if (b.get(k) !== v) return false;
  return true;
}

/** An axis waiting for its length, given one: the whole system, or a cyclic range over the origin in two pieces. */
function placedAxis(axis: Axis, lengths: ReadonlyMap<string, number>, range?: readonly number[]): Axis {
  if (axis.ranges) return axis;
  const length = lengths.get(axis.system);
  if (length === undefined) return axis;
  const ranges = range && range.length === 2 ? [range[0], length, 0, range[1]] : [0, length];
  return { ...axis, length, ranges };
}

/**
 * The resolved view with the lengths the store gave to axes the spec leaves unknown (_spec/README.md "ViewModel",
 * Placed): its sections' axes and their follows' leader axes. The same view object when nothing changes.
 */
function withLengths(
  root: Node,
  view: ResolvedView | undefined,
  lengths: ReadonlyMap<string, number>,
): ResolvedView | undefined {
  if (!view || !lengths.size) return view;
  let changed = false;
  /** A section's own range: on an axis waiting for its length, only a cyclic range over the origin has one. */
  const range = (address: string) => {
    const r = rangeParam(root, address);
    return r === 'all' ? undefined : r.split(',').map(Number);
  };
  const sections = view.sections.map((sec) => {
    const axis = placedAxis(sec.axis, lengths, range(sec.address));
    const leaderAxis = sec.follow ? placedAxis(sec.follow.leaderAxis, lengths, range(sec.follow.leader)) : undefined;
    if (axis === sec.axis && leaderAxis === sec.follow?.leaderAxis) return sec;
    changed = true;
    return { ...sec, axis, ...(sec.follow ? { follow: { ...sec.follow, leaderAxis: leaderAxis! } } : {}) };
  });
  return changed ? { ...view, sections } : view;
}

/** The view's sections with every leader before its followers (a chain of follows, EC-114); a cycle is cut. */
function followOrder(view: ResolvedView | undefined): string[] {
  const out: string[] = [];
  const done = new Set<string>();
  const byAddress = new Map((view?.sections ?? []).map((s) => [s.address, s] as const));
  const visit = (s: ResolvedSection, seen: Set<string>) => {
    if (done.has(s.address) || seen.has(s.address)) return;
    seen.add(s.address);
    const leader = s.follow && byAddress.get(s.follow.leader);
    if (leader) visit(leader, seen);
    done.add(s.address);
    out.push(s.address);
  };
  for (const s of view?.sections ?? []) visit(s, new Set());
  return out;
}
