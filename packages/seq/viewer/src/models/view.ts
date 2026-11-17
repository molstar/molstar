/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The view's model (_spec/README.md "ViewModel"): the view shown, a viewport per section, scroll groups on equal axes
// (EC-113), follows through the resolved view's follow path (EC-113, EC-114), viewports kept across updates (EC-111)
// and moves before a width kept as intents (EC-112). Session state: never in the spec or its history.

import { derive, KeyedState, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import {
  followRange,
  type ResolvedSection,
  type ResolvedView,
  resolveView,
  resolveViews,
} from '@molstar/seqx-core/view';
import { DefaultViewerConfig, type ViewerConfig } from '../config.js';
import { axesEqual, axisKey, axisTotal, nearestOffsetOf, offsetOf, positionAt } from './axis.js';
import type { SpecModel } from './spec.js';

/**
 * A section's viewport: `[start, end)` in offsets along its axis (fractional; positions on an axis that shows its
 * whole system from 0), its vertical scroll, its width in pixels. `pending` while it cannot be placed (no width yet,
 * or the axis waits for its length): the range is then the intent, applied when it can be.
 */
export interface Viewport {
  readonly start: number;
  readonly end: number;
  readonly top: number;
  readonly width: number;
  readonly pending: boolean;
  /** The section it follows, when it does (followers never lead, EC-114). */
  readonly leader?: string;
}

interface Group {
  readonly key: string;
  readonly axis: ResolvedSection['axis'];
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

const viewportsEqual = (a: Viewport, b: Viewport) =>
  a.start === b.start &&
  a.end === b.end &&
  a.top === b.top &&
  a.width === b.width &&
  a.pending === b.pending &&
  a.leader === b.leader;

export class ViewModel extends ReactiveModel {
  private readonly _view = new State<string | undefined>(undefined);
  private readonly _viewports = new KeyedState<string, Viewport>(viewportsEqual);
  private _sections = new Map<string, SectionState>();
  private _groups = new Map<string, Group>();
  /** Sections in an order where each leader comes before its followers. */
  private _order: string[] = [];
  readonly config: ViewerConfig;
  /** The address of the view shown; undefined: the spec's first view. */
  readonly view: ReadonlyState<string | undefined> = this._view.asReadonly;
  /** The resolved view shown, undefined when the spec has none. */
  readonly resolved: ReadonlyState<ResolvedView | undefined>;

  constructor(spec: SpecModel, config: ViewerConfig = DefaultViewerConfig) {
    super();
    this.config = config;
    this.resolved = derive([spec.model, this._view], (model, address) =>
      address === undefined ? resolveViews(model)[0] : resolveView(model, address),
    );
    this._update(this.resolved.value);
    this.listenPermanently(this.resolved, (v) => this._update(v));
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
    if (!s || s.width === width) return;
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
    if (!s || s.top === top) return;
    s.top = top;
    this._publish();
  }

  /** Shows `[start, end)` (offsets); before a width, the intent (EC-112). False when the section cannot move. */
  setRange(section: string, start: number, end: number): boolean {
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
    const g = this._movable(section);
    if (!g || delta === 0) return !!g;
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
    const g = this._movable(section);
    if (!g || !(factor > 0)) return false;
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
    if (s.section.follow) return this.attach(section);
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

  // ---- internals ----

  private _groupOf(s: SectionState): Group | undefined {
    return s.group === undefined ? undefined : this._groups.get(s.group);
  }

  private _groupWidth(g: Group): number {
    let w = 0;
    for (const m of g.members) w = Math.max(w, this._sections.get(m)?.width ?? 0);
    return w;
  }

  /** The group a move on the section moves: none for an attached follower that ignores moves (EC-114). */
  private _movable(section: string): Group | undefined {
    const s = this._sections.get(section);
    if (!s) return undefined;
    if (s.section.follow && s.group === undefined) {
      if (this.config.follow.followerMoves !== 'detach') return undefined;
      // detached: its own viewport, from where it is now
      const v = this.get(section);
      const key = `#${section}`;
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
        const kept = old?.group !== undefined && old.group.startsWith('#') ? oldGroups.get(old.group) : undefined;
        if (kept && axesEqual(kept.axis, section.axis)) {
          groups.set(kept.key, kept);
          state.group = kept.key;
        }
        continue;
      }
      const key = this.config.viewport.groupEqualAxes ? axisKey(section.axis) : section.address;
      state.group = key;
      const existing = groups.get(key);
      if (existing) {
        existing.members.push(section.address);
        continue;
      }
      groups.set(key, this._carry(key, section, oldGroups.get(key), old && oldGroups.get(old.group ?? '')));
    }
    this._sections = sections;
    this._groups = groups;
    this._order = followOrder(view);
    for (const address of this._viewports.values.keys()) if (!sections.has(address)) this._viewports.delete(address);
    this._publish();
  }

  /**
   * A group's viewport across an update (EC-111): the same axis keeps it; a changed one keeps its center and zoom when
   * the center is still shown, else fits. A new group fits.
   */
  private _carry(key: string, section: ResolvedSection, same: Group | undefined, before: Group | undefined): Group {
    const axis = section.axis;
    const total = axisTotal(axis) ?? 0;
    const group: Group = {
      key,
      axis,
      start: 0,
      end: total,
      fitted: true,
      placedWidth: 0,
      members: [section.address],
    };
    const prev = same && axesEqual(same.axis, axis) ? same : before;
    if (!prev) return group;
    group.placedWidth = prev.placedWidth;
    if (prev.fitted) return group;
    if (prev === same) return { ...group, start: prev.start, end: prev.end, fitted: false };
    const center = positionAt(prev.axis, (prev.start + prev.end) / 2);
    const c = center === undefined || prev.axis.system !== axis.system ? undefined : offsetOf(axis, center);
    if (c === undefined) return group;
    const half = (prev.end - prev.start) / 2;
    return { ...group, start: c - half, end: c + half, fitted: false };
  }

  /** Clamps each group where it can be placed, then writes every section's viewport, leaders before followers. */
  private _publish(): void {
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
      const total = axisTotal(s.section.axis);
      let next: Viewport | undefined;
      if (g) {
        const pending = total === undefined || this._groupWidth(g) <= 0;
        next = { start: g.start, end: g.end, top: s.top, width: s.width, pending };
      } else next = this._followed(s, total);
      if (next) this._set(address, next);
    }
  }

  private _set(address: string, next: Viewport): void {
    const prev = this._viewports.get(address);
    if (!prev || !viewportsEqual(prev, next)) this._viewports.set(address, next);
  }

  /**
   * A follower's viewport: its leader's range carried along the follow path at its scale (1, 3 or 1/3), continuously
   * within a run (`point`), else the hull of the integer image (`followRange`); its own top and width.
   */
  private _followed(s: SectionState, total: number | undefined): Viewport | undefined {
    const follow = s.section.follow!;
    const lead = this.get(follow.leader);
    const prev = this.get(s.section.address);
    if (!lead) return prev;
    const leaderAxis = follow.leaderAxis;
    const ps = positionAt(leaderAxis, lead.start);
    const pe = positionAt(leaderAxis, lead.end);
    let range: [number, number] | undefined;
    if (ps !== undefined && pe !== undefined) {
      const a = follow.point(ps);
      const b = endPoint(follow, pe);
      if (a !== undefined && b !== undefined) range = [Math.min(a, b), Math.max(a, b)];
      else {
        const hull = followRangeOf(s.section, ps, pe);
        if (hull) range = hull;
      }
    }
    const axis = s.section.axis;
    let start = prev?.start ?? 0;
    let end = prev?.end ?? total ?? 0;
    if (range) {
      const a = nearestOffsetOf(axis, range[0]);
      const b = nearestOffsetOf(axis, range[1]);
      if (a !== undefined && b !== undefined) [start, end] = [Math.min(a, b), Math.max(a, b)];
    }
    const pending = lead.pending || total === undefined || s.width <= 0;
    return { start, end, top: s.top, width: s.width, pending, leader: follow.leader };
  }
}

type Follow = NonNullable<ResolvedSection['follow']>;

/** The image of a range's end: the far edge of the element before it when it is whole (`point` reads its start). */
function endPoint(follow: Follow, x: number): number | undefined {
  if (!Number.isInteger(x)) return follow.point(x);
  const p = follow.point(x - 1);
  return p === undefined ? undefined : p + follow.scale * follow.dir;
}

/** The hull of the integer image of `[start, end)` positions, for ends that map nowhere. */
function followRangeOf(section: ResolvedSection, start: number, end: number): [number, number] | undefined {
  const pieces = section.follow ? followRange(section, [Math.floor(start), Math.ceil(end)]) : undefined;
  return pieces && pieces.length >= 2 ? [pieces[0], pieces[pieces.length - 1]] : undefined;
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
