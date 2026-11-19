/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The view's model (_spec/README.md "ViewModel"): the view shown, a viewport per section, scroll groups on equal axes
// of one layout (EC-113), each move's arithmetic its section's layout's (`SectionNav`, through the placement the last
// frame gave it), follows through the resolved view's follow path (EC-113, EC-114), a move asked of a follower carried back
// to its leader (`follow.followerMoves`: `lead`), viewports kept across updates (EC-111) and moves before a width kept
// as intents (EC-112), and the view's history of jumps (back and forward). Session state: never in the spec or its
// history.

import { batch, derive, KeyedState, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import { type Node, parseAddress } from '@molstar/core/repr/tree';
import {
  BuiltInLayouts,
  type NavContext,
  type Placement,
  type Viewport as RenderViewport,
  type SectionLayoutType,
  type SectionNav,
  type StepDir,
} from '@molstar/seqx-core/render';
import type { BandViewport } from '@molstar/seqx-core/render';
import { inversePath, type Path, projectEnd, projectPoint } from '@molstar/seqx-core/resolve';
import {
  type Axis,
  type Follow,
  followWindow,
  type ResolvedSection,
  type ResolvedView,
  resolveView,
  resolveViews,
} from '@molstar/seqx-core/view';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { axesEqual, axisKey, axisTotal, breakHolding, nearestOffsetOf, offsetOf, positionAt } from './axis.js';
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
  /**
   * px a position, where the layout reads one (text: the letter size; none: its preset's). Each layout reads the
   * record its way (_spec/README.md "ViewModel", Layouts): text's `start` is its anchor, its `end` an estimate.
   */
  readonly scale?: number;
}

/**
 * Where the models read what a layout's arithmetic needs (_spec/README.md "ViewModel", Layouts): a section's placement
 * in the last frame (text: its lines and their heights), the layouts by name, the least px a position (text: the
 * letters' threshold).
 */
export interface LayoutSource {
  placement(section: string): Placement | undefined;
  readonly layouts?: Readonly<Record<string, SectionLayoutType<any>>>;
  readonly minScale?: number;
}

interface Group {
  readonly key: string;
  readonly axis: Axis;
  /** Its sections' layout: a group holds one (scroll groups across layouts would read one range two ways). */
  readonly layout: string;
  start: number;
  end: number;
  /** px a position, where the layout reads one (text). */
  scale?: number;
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

/**
 * A group's range in a step of the history (_spec/README.md "ViewModel", History): its offsets with its axis's
 * signature (system, pieces, direction), its centre as a position and its members, so a step over a changed axis is
 * restored by position (the centre and the zoom), never dropped.
 */
interface Kept {
  readonly start: number;
  readonly end: number;
  readonly fitted: boolean;
  readonly sig?: string;
  readonly center?: number;
  readonly members?: readonly string[];
  readonly scale?: number;
}

/**
 * Every band's own viewport by its track's address (_spec/README.md "ViewModel", Bands): a rows band's scroll
 * (`bandTop`, display rows below its pinned band, fractional), row zoom (`rowScale`, px a row) and folded clades
 * (`folds`, session state). One object, replaced on each change.
 */
export type Bands = ReadonlyMap<string, BandViewport>;

const NoBands: Bands = new Map();

/**
 * A band as a step of the history keeps it (_spec/README.md "ViewModel", History): its top row by key (a row's key, and
 * the fraction of it scrolled past) with its row zoom and folds, so a step back lands on the same rows whatever was
 * folded, hidden or reordered since; `bandTop` as it was, where no layout names the rows (`BandRows`).
 */
export interface KeptBand {
  readonly band: BandViewport;
  readonly key?: string;
  readonly offset?: number;
}

/**
 * What names a band's rows (the rows layer's last layout, given by the viewer): the key of the row at a band's top
 * and the fraction into it, and the top (display rows) a key stands at under a band's folds and row zoom, clamped (a
 * folded key at its summary row, a hidden one at the nearest row shown).
 */
export interface BandRows {
  keyAt(track: string, band: BandViewport): { readonly key: string; readonly offset: number } | undefined;
  topOf(track: string, band: BandViewport, key: string, offset: number): number | undefined;
  /** A rows system as an axis in its rows layer's display order (what Shift+click extends along). */
  axisOf(system: string): Axis | undefined;
}

/** The leading viewports at one time: each group's range by key (a detached follower's among them), and the bands. */
interface Snapshot {
  readonly groups: ReadonlyMap<string, Kept>;
  readonly bands: ReadonlyMap<string, KeptBand>;
}

const sameBand = (a: KeptBand | undefined, b: KeptBand | undefined) =>
  !!a &&
  !!b &&
  a.band.rowScale === b.band.rowScale &&
  a.band.folds === b.band.folds &&
  (a.key !== undefined && b.key !== undefined
    ? a.key === b.key && a.offset === b.offset
    : a.band.bandTop === b.band.bandTop);

/** What a jump did: moved the view, found it already there, or could not move. */
export type JumpOutcome = 'moved' | 'same' | 'none';

/** The view's history of jumps: how many steps back and forward there are. */
export interface ViewHistory {
  readonly back: number;
  readonly forward: number;
}

/** Two steps are the same view: the same groups at the same offsets, compared only on equal axes. */
function snapshotsEqual(a: Snapshot, b: Snapshot): boolean {
  if (a.groups.size !== b.groups.size || a.bands.size !== b.bands.size) return false;
  for (const [k, v] of a.groups) {
    const w = b.groups.get(k);
    if (!w || w.sig !== v.sig || w.start !== v.start || w.end !== v.end || w.fitted !== v.fitted || w.scale !== v.scale)
      return false;
  }
  for (const [k, v] of a.bands) if (!sameBand(v, b.bands.get(k))) return false;
  return true;
}

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

  /** The jumps' history (_spec/README.md "ViewModel", History): the viewports each jump left, and those back left. */
  private readonly _back: Snapshot[] = [];
  private readonly _forward: Snapshot[] = [];
  private readonly _history = new State<ViewHistory>(
    { back: 0, forward: 0 },
    (a, b) => a.back === b.back && a.forward === b.forward,
  );
  private readonly _jumps = new State(0);
  private readonly _bands = new State<Bands>(NoBands);
  /** Every band's own viewport by its track's address: a rows band's scroll, row zoom and folds. */
  readonly bands: ReadonlyState<Bands> = this._bands.asReadonly;
  /** How many steps back and forward the view's history holds. */
  readonly history: ReadonlyState<ViewHistory> = this._history.asReadonly;
  /** Moves on each jump, back and forward (and each step zoom): what a viewer animates; a move by hand does not. */
  readonly jumps: ReadonlyState<number> = this._jumps.asReadonly;

  /** Where a section's placement comes from: its moves are then its layout's arithmetic (`setLayoutSource`). */
  private _layouts: LayoutSource | undefined;

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

  /**
   * Where a section's placement in the last frame comes from (the viewer's frames): each move asked of a section is
   * then its layout's arithmetic over it (_spec/README.md "ViewModel", Layouts); without one, a window of offsets.
   */
  setLayoutSource(source: LayoutSource): void {
    this._layouts = source;
  }

  /**
   * A frame placed the sections again: a layout's estimate drawn from its placement (a text's end: the lines in view)
   * is read again. Tells only what changed.
   */
  placed(): void {
    if (this._layouts) this._publish();
  }

  /** A section's layout's name (`tracks` when it has none). */
  layoutOf(section: string): string {
    return this._sections.get(section)?.section.layout ?? 'tracks';
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

  /**
   * Sets a leading section's viewport as given (offsets: a story step's view, an animation left where it was drawn):
   * its window, its letter size and its top, then clamped as any move; no history step (_spec/README.md "ViewModel",
   * Stories). False when the section does not lead (a follower follows) or cannot move.
   */
  setViewport(
    section: string,
    v: { readonly start: number; readonly end: number; readonly top?: number; readonly scale?: number },
  ): boolean {
    if (!Number.isFinite(v.start) || !Number.isFinite(v.end) || this.isFollowing(section)) return false;
    const g = this._movable(section);
    if (!g) return false;
    this._take(g, {
      start: Math.min(v.start, v.end),
      end: Math.max(v.start, v.end),
      ...(v.scale !== undefined ? { scale: v.scale } : {}),
    });
    const s = this._sections.get(section);
    if (s && v.top !== undefined && Number.isFinite(v.top)) s.top = v.top;
    this._publish();
    return true;
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
    const [s, e] = [Math.min(start, end), Math.max(start, end)];
    const c = this._navContext(section, g);
    const focus = c && this._navOf(g.layout).focus;
    if (c && focus) this._take(g, focus(c, rangeOf(g), s, e));
    else {
      g.start = s;
      g.end = e;
      g.fitted = false;
    }
    this._publish();
    return true;
  }

  /**
   * The least move that shows `[start, end)` (offsets): the cursor out of view, a linked selection (the layout's
   * `reveal`: a window pans to its nearer edge, a text scrolls the lines in). False when nothing moved.
   */
  reveal(section: string, start: number, end: number): boolean {
    if (!Number.isFinite(start) || !Number.isFinite(end) || this._leading(section)) return false;
    const g = this._movable(section);
    const c = g && this._navContext(section, g);
    if (!g || !c) return false;
    const next = this._navOf(g.layout).reveal(c, rangeOf(g), Math.min(start, end), Math.max(start, end));
    if (next.start === g.start && next.end === g.end) return false;
    this._take(g, next);
    this._publish();
    return true;
  }

  /**
   * Where a cursor at offset `offset` of a section steps (`up`, `down`: a line, in a layout of lines), by its layout's
   * arithmetic; −1 for none (tracks have one line: up and down leave the section).
   */
  step(section: string, offset: number, dir: StepDir): number {
    const s = this._sections.get(section);
    const c = s && this._navContext(section);
    return s && c ? this._navOf(s.section.layout).step(c, offset, dir) : -1;
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

  /**
   * Moves the view by a scroll's deltas, px (positive shows what is right and below), by its layout's arithmetic: a
   * window pans across, a text scrolls its lines; nothing before a width. A follower that leads moves its leader.
   */
  panPx(section: string, dx: number, dy = 0): boolean {
    const v = this.get(section);
    if (!v || v.width <= 0) return false;
    const g = this._leading(section) ? undefined : this._movable(section);
    const c = g && this._navContext(section, g);
    if (!g || !c) return dx ? this.pan(section, (dx * (v.end - v.start)) / v.width) : this._canMove(section);
    this._take(g, this._navOf(g.layout).pan(c, rangeOf(g), dx, dy));
    this._publish();
    return true;
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
    const c = this._navContext(section, g);
    if (c) this._take(g, this._navOf(g.layout).zoom(c, rangeOf(g), factor, a));
    else {
      // before a width: the intent, a window of offsets
      g.start = a - (a - g.start) / factor;
      g.end = a + (g.end - a) / factor;
      g.fitted = false;
    }
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
    // its layout's fit (text: the preset's letter size)
    const c = this._navContext(section, g);
    if (c) {
      const f = this._navOf(g.layout).fit(c, rangeOf(g));
      if (f.scale !== undefined) g.scale = f.scale;
      else delete g.scale;
    }
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

  // ---- the history of jumps ----

  /**
   * Runs a jump (a focus, a fit, a reveal, a bookmark's or a link's move): when it moves the view, the viewports it left
   * are pushed on the history (`record`, the default; a step zoom animates without one) and the steps forward dropped.
   * `same` when the view already showed what was asked (a double-click on what is focused), `none` when nothing could
   * move. A pan or a zoom by hand pushes nothing: a jump after it pushes where the hand left the view.
   */
  jump(move: () => boolean, record = true): JumpOutcome {
    const before = this._snapshot();
    if (!move()) return 'none';
    if (snapshotsEqual(before, this._snapshot())) return 'same';
    if (record) {
      this._back.push(before);
      const limit = Math.max(0, this.config.view.historyLimit);
      if (this._back.length > limit) this._back.splice(0, this._back.length - limit);
      this._forward.length = 0;
      this._tellHistory();
    }
    this._jumps.set(this._jumps.value + 1);
    return 'moved';
  }

  /** What names the bands' rows for the history (the viewer's rows layers); none: a band kept as its `bandTop`. */
  private _bandRows: BandRows | undefined;

  setBandRows(rows: BandRows | undefined): void {
    this._bandRows = rows;
  }

  /** A rows system as an axis in display order, when a rows layer shows it; undefined else. */
  rowsAxis(system: string): Axis | undefined {
    return this._bandRows?.axisOf(system);
  }

  /** A band's own viewport (its track's address), undefined when it has none yet. */
  band(track: string): BandViewport | undefined {
    return this._bands.value.get(track);
  }

  /**
   * Sets a band's own viewport (a rows band's scroll, row zoom, folds), merged over what it holds; undefined fields are
   * left out (`folds: undefined` drops the session's folds: the layer's `collapsed` shows). Not a jump: a hand's move,
   * kept by the next jump's step (a bulk fold runs through `jump`).
   */
  setBand(track: string, band: Partial<BandViewport>): void {
    const prev = this._bands.value.get(track);
    const next: { bandTop: number; rowScale?: number; folds?: readonly string[] } = {
      bandTop: Math.max(0, band.bandTop ?? prev?.bandTop ?? 0),
    };
    const scale = 'rowScale' in band ? band.rowScale : prev?.rowScale;
    const folds = 'folds' in band ? band.folds : prev?.folds;
    if (scale !== undefined) next.rowScale = scale;
    if (folds !== undefined) next.folds = folds;
    if (prev && prev.bandTop === next.bandTop && prev.rowScale === next.rowScale && prev.folds === next.folds) return;
    const map = new Map(this._bands.value);
    map.set(track, next);
    this._bands.set(map);
  }

  /** The bands of a section's tracks (by track address), as a frame reads them; undefined when none has one. */
  bandsOf(section: string): Readonly<Record<string, BandViewport>> | undefined {
    const all = this._bands.value;
    if (!all.size) return undefined;
    const kept = this._bandsKept.get(section);
    if (kept && kept[0] === all) return kept[1];
    let out: Record<string, BandViewport> | undefined;
    for (const [track, b] of all)
      if (track.startsWith(`${section}/`)) {
        out ??= {};
        out[track] = b;
      }
    this._bandsKept.set(section, [all, out]);
    return out;
  }

  private readonly _bandsKept = new Map<string, [Bands, Readonly<Record<string, BandViewport>> | undefined]>();

  /** Whether there is a step back or forward. */
  get canBack(): boolean {
    return this._back.length > 0;
  }

  get canForward(): boolean {
    return this._forward.length > 0;
  }

  /** Returns to the viewports the last jump left, every leading section's together; false when there is none. */
  back(): boolean {
    return this._step(this._back, this._forward);
  }

  /** Goes again where the last back left. */
  forward(): boolean {
    return this._step(this._forward, this._back);
  }

  /** Forgets the history (a new spec loaded). */
  clearHistory(): void {
    this._back.length = 0;
    this._forward.length = 0;
    this._tellHistory();
  }

  private _step(from: Snapshot[], to: Snapshot[]): boolean {
    // a step to where the view is now (a group gone, a hand that returned) is skipped
    const now = this._snapshot();
    let target: Snapshot | undefined;
    while (from.length) {
      const t = from.pop()!;
      if (!snapshotsEqual(this._restored(t), now)) {
        target = t;
        break;
      }
    }
    if (!target) {
      this._tellHistory();
      return false;
    }
    to.push(now);
    this._restore(target);
    this._tellHistory();
    this._jumps.set(this._jumps.value + 1);
    return true;
  }

  private _tellHistory(): void {
    this._history.set({ back: this._back.length, forward: this._forward.length });
  }

  /** Each group's range now, with its axis's signature, its centre as a position and its members. */
  private _snapshot(): Snapshot {
    const out = new Map<string, Kept>();
    for (const g of this._groups.values()) {
      const center = g.axis.ranges ? positionAt(g.axis, (g.start + g.end) / 2) : undefined;
      out.set(g.key, {
        start: g.start,
        end: g.end,
        fitted: g.fitted,
        sig: sigOf(g),
        ...(center !== undefined ? { center } : {}),
        members: g.members,
        ...(g.scale !== undefined ? { scale: g.scale } : {}),
      });
    }
    // the bands with the windows, by their top rows' keys: back after a fold or a scroll lands on the same rows
    const bands = new Map<string, KeptBand>();
    for (const [track, band] of this._bands.value) {
      const at = this._bandRows?.keyAt(track, band);
      bands.set(track, at ? { band, key: at.key, offset: at.offset } : { band });
    }
    return { groups: out, bands };
  }

  /**
   * A group's kept range in a step: by its key, else (its axis changed since: a reverse, a collapse) by a member; on an
   * axis of another signature, restored by position: the centre and the zoom, the centre clamped to what is drawn (a
   * centre now hidden goes to its break). Undefined when the step has none of it.
   */
  private _keptOf(snap: Snapshot, g: Group): Kept | undefined {
    let v = snap.groups.get(g.key);
    if (!v)
      for (const w of snap.groups.values())
        if (w.members?.some((m) => g.members.includes(m))) {
          v = w;
          break;
        }
    const sig = sigOf(g);
    if (!v || v.sig === undefined || v.sig === sig) return v;
    // the letter size is kept within one layout only (a text's, across a reverse or a collapse)
    const scale = v.scale !== undefined && layoutOfSig(v.sig) === g.layout ? { scale: v.scale } : {};
    const total = axisTotal(g.axis);
    if (v.fitted || total === undefined) return { start: 0, end: total ?? v.end, fitted: true, sig, ...scale };
    const c = v.center === undefined ? undefined : centerOffset(g.axis, v.center);
    if (c === undefined) return { start: 0, end: total, fitted: true, sig, ...scale };
    // the centre and the zoom: a window's middle, a text's lines in view (its anchor and its estimated end) alike
    const half = (v.end - v.start) / 2;
    return { start: c - half, end: c + half, fitted: false, sig, ...scale };
  }

  /**
   * What restoring a snapshot would give: the groups there now, at its ranges (by position onto a changed axis); a
   * follower detached since, attached; one detached then and attached since, detached again (a section gone since is
   * left out).
   */
  private _restored(snap: Snapshot): Snapshot {
    const out = new Map<string, Kept>();
    for (const g of this._groups.values()) {
      const v = this._keptOf(snap, g);
      if (v) out.set(g.key, { ...v, sig: sigOf(g) });
      else if (!g.key.startsWith(DetachedPrefix))
        out.set(g.key, {
          start: g.start,
          end: g.end,
          fitted: g.fitted,
          sig: sigOf(g),
          ...(g.scale !== undefined ? { scale: g.scale } : {}),
        });
    }
    for (const [key, v] of snap.groups)
      if (key.startsWith(DetachedPrefix) && !out.has(key) && this._detachable(key)) out.set(key, v);
    return { groups: out, bands: snap.bands };
  }

  /** The follower a detached group's key names, when it is in the view and follows. */
  private _detachable(key: string): SectionState | undefined {
    const s = this._sections.get(key.slice(DetachedPrefix.length));
    return s?.section.follow ? s : undefined;
  }

  /** The groups take the snapshot's ranges in one publish; a follower detached since it was taken is attached again. */
  private _restore(snap: Snapshot): void {
    // each band at its kept top row's place under its kept folds and zoom, clamped; a band it had none of has none
    const map = new Map<string, BandViewport>();
    for (const [track, k] of snap.bands) {
      const top = k.key !== undefined ? this._bandRows?.topOf(track, k.band, k.key, k.offset ?? 0) : undefined;
      map.set(track, top === undefined ? k.band : { ...k.band, bandTop: top });
    }
    if (map.size || this._bands.value.size) this._bands.set(map.size ? map : NoBands);
    for (const [key, g] of [...this._groups]) {
      const v = this._keptOf(snap, g);
      if (v) {
        g.start = v.start;
        g.end = v.end;
        g.fitted = v.fitted;
        if (v.scale !== undefined) g.scale = v.scale;
        else delete g.scale;
      } else if (key.startsWith(DetachedPrefix)) {
        this._groups.delete(key);
        for (const s of this._sections.values()) if (s.group === key) s.group = undefined;
      }
    }
    // a follower detached when the snapshot was taken is detached again, at its range then
    for (const [key, v] of snap.groups) {
      const s = !this._groups.has(key) && key.startsWith(DetachedPrefix) ? this._detachable(key) : undefined;
      if (!s) continue;
      this._groups.set(key, {
        key,
        axis: s.section.axis,
        layout: s.section.layout,
        start: v.start,
        end: v.end,
        fitted: v.fitted,
        ...(v.scale !== undefined ? { scale: v.scale } : {}),
        placedWidth: s.width,
        members: [s.section.address],
      });
      s.group = key;
    }
    this._publish();
  }

  // ---- internals ----

  /** A layout's navigation by name (tracks' when it is unknown). */
  private _navOf(layout: string): SectionNav {
    return this._layoutOf(layout).nav;
  }

  /** A layout by name: the layout source's (the frame's), else the built-ins; tracks' when it is unknown, as drawn. */
  private _layoutOf(layout: string): SectionLayoutType<any> {
    const layouts = this._layouts?.layouts ?? BuiltInLayouts;
    return layouts[layout] ?? layouts.tracks ?? BuiltInLayouts.tracks;
  }

  /**
   * Whether a section's window wraps its cyclic axis's origin (its layout's `wraps`: tracks and rings; a text's anchor
   * does not): one rule for its navigation, its window's drags and its input.
   */
  wraps(section: string): boolean {
    const s = this._sections.get(section);
    return !!s && s.section.axis.cyclic && this._layoutOf(s.section.layout).wraps === true;
  }

  /**
   * Whether two ends on a section are read as the shorter arc between them: a layout that draws the whole circle (one
   * with a `window`: a ring), or a wrapping window over the origin now (a linear view panned across it).
   */
  arcs(section: string): boolean {
    if (!this.wraps(section)) return false;
    const s = this._sections.get(section)!;
    if (this._layoutOf(s.section.layout).window) return true;
    const v = this.get(section);
    const total = axisTotal(s.section.axis);
    return !!v && total !== undefined && v.end > total;
  }

  /**
   * What a section's layout's arithmetic reads (_spec/README.md "ViewModel", Layouts): its axis's total, its group's
   * width (`g`) or its own, its placement in the last frame; undefined while it cannot be placed (no length or width):
   * a move is then an intent, a window of offsets.
   */
  private _navContext(section: string, g?: Group): NavContext | undefined {
    const s = this._sections.get(section);
    const total = s && axisTotal(s.section.axis);
    const width = g ? this._groupWidth(g) : (s?.width ?? 0);
    if (!s || total === undefined || !(width > 0)) return undefined;
    const placement = this._layouts?.placement(section);
    const minScale = this._layouts?.minScale;
    return {
      total,
      width,
      height: placement?.box.height ?? 0,
      // a cyclic axis's window passes its origin where its layout draws it there (`wraps`: tracks, rings); a text's
      // anchor keeps within the axis
      cyclic: s.section.axis.cyclic && this._layoutOf(g?.layout ?? s.section.layout).wraps === true,
      maxScale: this.config.viewport.maxPxPerPosition,
      ...(minScale !== undefined ? { minScale } : {}),
      ...(placement ? { placement } : {}),
    };
  }

  /** A group takes a viewport its layout's arithmetic gave. */
  private _take(g: Group, v: RenderViewport): void {
    g.start = v.start;
    g.end = v.end;
    if (v.scale !== undefined) g.scale = v.scale;
    else delete g.scale;
    g.fitted = false;
  }

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
    const l0 = positionAt(la, lead.start, 'start');
    const l1 = positionAt(la, lead.end, 'end');
    if (l0 === undefined || l1 === undefined) return false;
    // the leader's window in positions, lowest first (a reversed leader's start is its highest)
    const ps = Math.min(l0, l1);
    const pe = Math.max(l0, l1);
    if (delta !== undefined) {
      // the image of the leader's window before the follower's clamp, moved by the pan
      const w = followWindow(s.section, la, axis, lead.start, lead.end, this._arc());
      if (w && !w.wraps && w.start < w.end) [start, end] = [w.start + delta, w.end + delta];
    }
    // a cyclic follower's range over its origin runs past the total
    const pos = (o: number, side: 'start' | 'end') => positionAt(axis, axis.cyclic && o > total ? o - total : o, side);
    const f0 = pos(start, 'start');
    const f1 = pos(end, 'end');
    // the follower's ends in positions, lowest first, with the way each moved there (navigation is in display order:
    // on a reversed follower an end moved right moved down)
    const rev = axis.dir === -1;
    const lo = rev ? f1 : f0;
    const hi = rev ? f0 : f1;
    const mLo = rev ? -me : ms;
    const mHi = rev ? -ms : me;
    const inv = inverseOf(follow.path);
    const xLo = lo === undefined ? undefined : backTo(follow, inv, axis, lo, mLo);
    const xHi = hi === undefined ? undefined : backTo(follow, inv, axis, hi, mHi);
    // the leader's lowest position is the preimage of the follower's lowest, or of its highest on a minus strand
    let a = follow.dir === 1 ? xLo : xHi;
    let b = follow.dir === 1 ? xHi : xLo;
    if (delta !== undefined) {
      if (a === undefined && b === undefined) return false;
      if (a === undefined) a = ps + (b! - pe);
      if (b === undefined) b = pe + (a - ps);
    } else {
      // an end moved outward past every image: the leader's axis edge on that side (its clamp keeps it in)
      const lowOut = follow.dir === 1 ? mLo < 0 : mHi > 0;
      const highOut = follow.dir === 1 ? mHi > 0 : mLo < 0;
      const lt = axisTotal(la);
      const e0 = lt === undefined ? undefined : positionAt(la, 0, 'start');
      const e1 = lt === undefined ? undefined : positionAt(la, lt, 'end');
      const lowEdge = e0 === undefined || e1 === undefined ? ps : Math.min(e0, e1);
      const highEdge = e0 === undefined || e1 === undefined ? pe : Math.max(e0, e1);
      a ??= lowOut ? lowEdge : ps;
      b ??= highOut ? highEdge : pe;
    }
    const oa = nearestOffsetOf(la, a);
    const ob = nearestOffsetOf(la, b);
    if (oa === undefined || ob === undefined || oa === ob) return false;
    return this.setRange(follow.leader, Math.min(oa, ob), Math.max(oa, ob));
  }

  /** The arc a cyclic follower takes (`follow.cyclicArc`). */
  private _arc(): 'shortest' | 'hull' {
    return this.config.follow.cyclicArc === 'hull' ? 'hull' : 'shortest';
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
        layout: s.section.layout,
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
      // groups only among sections of one layout: a range is read its layout's way (_spec/README.md "ViewModel")
      const axisPart = section.axis.ranges
        ? axisKey(section.axis)
        : `${axisKey(section.axis)}|${rangeParam(root, section.address)}`;
      const key = this.config.viewport.groupEqualAxes ? `${axisPart}|${section.layout}` : section.address;
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
        const c = this._navContext(g.members[0], g);
        // a layout whose window is not its range (a text's lines in view, from its last frame): its end that window's
        // width on from its start, an estimate (lines at full height until laid out)
        const w = c?.placement?.window;
        if (g.fitted) {
          g.start = 0;
          g.end = w ? Math.min(total, Math.max(1, w.end - w.start)) : total;
          continue;
        }
        // the layout's clamp (a window within its axis and the most px a position; a text's anchor and size)
        if (c) {
          const v = this._navOf(g.layout).clamp(c, rangeOf(g));
          g.start = v.start;
          g.end = v.end;
          if (v.scale !== undefined) g.scale = v.scale;
          if (w) g.end = Math.min(total, g.start + Math.max(1, w.end - w.start));
        } else {
          const w = Math.min(Math.max(g.end - g.start, Math.min(width / maxPxPerPosition, total)), total);
          const m = (g.start + g.end) / 2;
          g.start = Math.min(Math.max(m - w / 2, 0), total - w);
          g.end = g.start + w;
        }
        if (g.start === 0 && g.end === total) g.fitted = true;
      }
      for (const address of this._order) {
        const s = this._sections.get(address)!;
        const g = this._groupOf(s);
        if (g) {
          const pending = axisTotal(s.section.axis) === undefined || this._groupWidth(g) <= 0;
          this._set(address, g.start, g.end, s, pending, undefined, undefined, g.scale);
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
    scale?: number,
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
      p.scale === scale &&
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
      ...(scale !== undefined ? { scale } : {}),
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
    const w =
      total === undefined
        ? undefined
        : followWindow(s.section, follow.leaderAxis, axis, lead.start, lead.end, this._arc());
    if (w && total !== undefined) {
      if (w.wraps) {
        // over the origin: on to the total, then from 0
        start = w.start;
        end = total + w.end;
        pieces = [w.start, total, 0, w.end];
      } else {
        start = w.start;
        end = w.end;
        pieces = undefined;
        if (s.width > 0) [start, end] = clampRange(start, end, total, s.width, this.config.viewport.maxPxPerPosition);
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
  const group: Group = {
    key,
    axis,
    layout: section.layout,
    start: 0,
    end: total,
    fitted: true,
    placedWidth: 0,
    members: [section.address],
  };
  if (!prev) return group;
  group.placedWidth = prev.placedWidth;
  // the letter size is kept within one layout (a text across a patch of its params)
  if (prev.scale !== undefined && prev.layout === section.layout) group.scale = prev.scale;
  if (prev.fitted) return group;
  if (axesEqual(prev.axis, axis)) return { ...group, start: prev.start, end: prev.end, fitted: false };
  if (prev.axis.system !== axis.system) return group;
  if (!prev.axis.ranges) {
    // a pending axis's offsets are positions from 0: kept where the new axis starts there
    return axis.ranges?.[0] === 0 ? { ...group, start: prev.start, end: prev.end, fitted: false } : group;
  }
  const center = positionAt(prev.axis, (prev.start + prev.end) / 2);
  const c = center === undefined ? undefined : centerOffset(axis, center);
  if (c === undefined) return group;
  const half = (prev.end - prev.start) / 2;
  return { ...group, start: c - half, end: c + half, fitted: false };
}

/**
 * A centre position's offset on an axis: where it is shown; in a stretch a break hides (a collapse), the break's (the
 * nearest piece edge, zero wide); outside the axis, none.
 */
function centerOffset(axis: Axis, center: number): number | undefined {
  const o = offsetOf(axis, center);
  if (o !== undefined) return o;
  return breakHolding(axis, Math.floor(center)) >= 0 ? nearestOffsetOf(axis, center) : undefined;
}

/** A group's range as its layout reads it. */
function rangeOf(g: Group): RenderViewport {
  return g.scale !== undefined ? { start: g.start, end: g.end, scale: g.scale } : { start: g.start, end: g.end };
}

/** A group's signature in a step of the history: its axis's and its layout (a step over either is restored by position). */
function sigOf(g: Group): string {
  return `${axisKey(g.axis)}#${g.layout}`;
}

function layoutOfSig(sig: string): string {
  return sig.slice(sig.lastIndexOf('#') + 1);
}

/** How many follows a move asked of a follower goes up at most (a chain of follows; a cycle is cut). */
const MaxFollowDepth = 16;

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
 * `y`: the inverse of the forward map itself, so a follower moved through it shows what was asked. The forward map is
 * place by place (a residue's three fine places across a codon step: a codon split by an intron is its nucleotides),
 * each place's image a run of `point(k) + f · slope`, `f` in [0, 1); the places tried are those near `y`'s preimage
 * through `inv` (the follow read back) and near the nearest mapped positions on either side. Where `y` is in no place's
 * image (an intron, a UTR), the nearest image the way the end moved (`m`), so a pan steps over a gap and never sticks.
 * Undefined when nothing lies that way.
 */
function backTo(follow: Follow, inv: Path, axis: Axis, y: number, m: number): number | undefined {
  // a place's width in leader positions: a third of a residue across a codon step
  const u = follow.path.step === 1 ? 1 / 3 : 1;
  const tried = new Set<number>();
  let best: number | undefined;
  let distance = Infinity;
  const near = (x: number | undefined) => {
    if (x === undefined) return undefined;
    const k0 = Math.floor(x / u + 1e-9);
    for (let k = k0 - 1; k <= k0 + 1; k++) {
      if (tried.has(k)) continue;
      tried.add(k);
      const a = follow.point(k * u);
      const half = follow.point((k + 0.5) * u);
      if (a === undefined || half === undefined || half === a) continue;
      const slope = 2 * (half - a);
      const f = (y - a) / slope;
      if (f >= 0 && f < 1) return (k + f) * u;
      if (m === 0) continue;
      // the end of this place's image nearer `y`, when it lies the way the end moved
      const b = a + slope;
      const v = m > 0 ? Math.min(a, b) : Math.max(a, b);
      const d = (v - y) * m;
      if (d > 0 && d < distance) {
        distance = d;
        best = (k + (v === a ? 0 : BeforeEdge)) * u;
      }
    }
    return undefined;
  };
  return (
    near(projectPoint(inv, y)) ??
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
    return run ? projectPoint(inv, Math.ceil(run.from / scale)) : undefined;
  }
  const to = Math.floor(y);
  if (to <= lo) return undefined;
  const runs = inv.linear.forward(lo * scale, to * scale);
  const last = runs[runs.length - 1];
  return last ? projectEnd(inv, Math.floor((last.from + last.length - 1) / scale) + 1) : undefined;
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
