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
// history. This module holds the moves and the publish; beside it, the scroll groups across updates (`view-groups.ts`),
// the bands (`view-bands.ts`), the follows' arithmetic (`view-follow.ts`) and the history's steps (`view-history.ts`).

import { batch, derive, KeyedState, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import {
  type NavContext,
  type Placement,
  type Viewport as RenderViewport,
  type SectionLayoutType,
  type SectionNav,
  type StepDir,
  TracksLayout,
} from '@molstar/seqx-core/render';
import type { Location } from '@molstar/seqx-core/query';
import type { BandViewport } from '@molstar/seqx-core/render';
import {
  type Axis,
  blockSpan,
  BuiltInLayerKinds,
  followWindow,
  type LayerKinds,
  opensComplexes,
  type ResolvedView,
  resolveView,
  resolveViews,
} from '@molstar/seqx-core/view';
import { openedModel, type Overlay, type StaticModel } from '@molstar/seqx-core/resolve';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { axisTotal } from './axis.js';
import type { SpecModel } from './spec.js';
import { type BandRows, type Bands, ViewBands } from './view-bands.js';
import {
  type FollowerRangeOut,
  followerRange,
  followOrder,
  type LeaderRangeOut,
  leaderRange,
  MaxFollowDepth,
} from './view-follow.js';
import { DetachedPrefix, type Group, groupSections, rangeOf, type SectionState } from './view-groups.js';
import {
  type JumpKind,
  type JumpOutcome,
  JumpSteps,
  type KeptBand,
  type LastJump,
  restoredOf,
  restoreInto,
  type Snapshot,
  snapshotOf,
  snapshotsEqual,
  type ViewHistory,
  type ViewSnapshot,
} from './view-history.js';

export type { BandRows, Bands, JumpKind, JumpOutcome, KeptBand, LastJump, ViewHistory, ViewSnapshot };

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
  /** A grid's own y window (`follow_y: none`), offsets down its y axis; absent: the whole axis, or y follows. */
  readonly y?: { readonly start: number; readonly end: number };
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

export class ViewModel extends ReactiveModel {
  private readonly _view = new State<string | undefined>(undefined);
  private readonly _viewports = new KeyedState<string, Viewport>();
  private _sections = new Map<string, SectionState>();
  private _groups = new Map<string, Group>();
  /** Sections in an order where each leader comes before its followers. */
  private _order: string[] = [];
  /** The records the follow functions write into: this model's own, so a move and a publish make none. */
  private readonly _leaderOut: LeaderRangeOut = { start: 0, end: 0 };
  private readonly _followerOut: FollowerRangeOut = { start: 0, end: 0, pieces: undefined };
  private readonly _spec: SpecModel;
  readonly config: ViewerConfig;
  /** The address of the view shown; undefined: the spec's first view. */
  readonly view: ReadonlyState<string | undefined> = this._view.asReadonly;
  /**
   * The model the view is resolved over: the spec's, or with a section of a source's complex the opened model over the
   * store's overlay (core's _spec "Overlay", Opened sources), one a node version.
   */
  readonly model: ReadonlyState<StaticModel>;
  /** The resolved view shown, undefined when the spec has none. */
  readonly resolved: ReadonlyState<ResolvedView | undefined>;

  /** The jumps' history (_spec/README.md "ViewModel", History): the viewports each jump left, and those back left. */
  private readonly _steps = new JumpSteps();
  private readonly _jumps = new State(0);
  private readonly _bands = new ViewBands();
  /** Every band's own viewport by its track's address: a rows band's scroll, row zoom and folds. */
  readonly bands: ReadonlyState<Bands> = this._bands.bands;
  /** How many steps back and forward the view's history holds. */
  readonly history: ReadonlyState<ViewHistory> = this._steps.history;
  /** Moves on each jump, back and forward (and each step zoom): what a viewer animates; a move by hand does not. */
  readonly jumps: ReadonlyState<number> = this._jumps.asReadonly;
  /** Where the last jump left the view, while no hand moved it: the step a jump pushes keeps its Location. */
  private _landed: Snapshot | undefined;
  private _lastJump: LastJump | undefined;

  /** Where a section's placement comes from: its moves are then its layout's arithmetic (`setLayoutSource`). */
  private _layouts: LayoutSource | undefined;

  /** The layer kinds the view is resolved with (the composition's: `setKinds` when an extension comes or goes). */
  private readonly _kinds: State<LayerKinds>;

  /**
   * `overlay`: the store's overlay over the spec's model (the viewer's, one a node version); a section of a source's
   * complex resolves over it, and reflows when it moves. None: the spec's model alone.
   */
  constructor(
    spec: SpecModel,
    config: ViewerConfig = DefaultViewerConfig,
    kinds: LayerKinds = BuiltInLayerKinds,
    overlay?: ReadonlyState<Overlay>,
  ) {
    super();
    this._spec = spec;
    this.config = config;
    this._kinds = new State(kinds);
    // a spec with no section of a source's complex never reads the overlay: no view resolved again a node version
    this.model = overlay
      ? derive([spec.model, overlay], (m, o) => (opensComplexes(m) ? openedModel(m, o) : m))
      : spec.model;
    // a source-backed system's length is the model's once its source opened (`SpecModel.setLengths`)
    this.resolved = derive([this.model, this._view, this._kinds], (model, address, k) =>
      address === undefined ? resolveViews(model, k)[0] : resolveView(model, address, k),
    );
    this._update(this.resolved.value);
    this.listenPermanently(this.resolved, (v) => this._update(v));
  }

  /**
   * The layer kinds the view is resolved with: the view is resolved again (a new `resolved`) when they change.
   * @internal
   */
  setKinds(kinds: LayerKinds): void {
    this._kinds.set(kinds);
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
    this._stopFree = true;
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
    // a range set (a jump, a story step) may move to another block under `navigate.acrossBlocks: stop`
    this._stopFree = true;
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
      if (!v) return false;
      // at the follower's end the pan stops there, its width kept (a leader asked past it would narrow instead)
      const total = lead.section.axis.cyclic ? undefined : axisTotal(lead.section.axis);
      const d = total === undefined ? delta : stopAtEnds(delta, v.start, v.end, total);
      if (!d) return false;
      const m = Math.sign(d);
      return this._lead(lead, v.start + d, v.end + d, m, m, d);
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
   * window pans across (a row of blocks keeping the content under the pointer's `x`, px along the axis), a text scrolls its lines; nothing before a width. A follower that leads moves its leader.
   */
  panPx(section: string, dx: number, dy = 0, x?: number): boolean {
    const v = this.get(section);
    if (!v || v.width <= 0) return false;
    // a grid whose y moves apart from x: dy moves y (its own window, or the section it follows), dx moves x
    if (this._yOf(section)) {
      const moved = dy ? this._panY(section, dy) : false;
      return (dx ? this._panX(section, v, dx, 0, x) : false) || moved;
    }
    return this._panX(section, v, dx, dy, x);
  }

  private _panX(section: string, v: Viewport, dx: number, dy: number, x?: number): boolean {
    const g = this._leading(section) ? undefined : this._movable(section);
    const c = g && this._navContext(section, g);
    if (!g || !c) return dx ? this.pan(section, (dx * (v.end - v.start)) / v.width) : this._canMove(section);
    this._take(g, this._navOf(g.layout).pan(c, rangeOf(g), dx, dy, x));
    this._publish();
    return true;
  }

  /**
   * Zooms by `factor` (above 1: in) about `anchor` (an offset; the center when absent), which stays in place; a grid
   * whose y moves apart from x zooms y too, about `anchorY` (an offset down y).
   */
  zoom(section: string, factor: number, anchor?: number, anchorY?: number): boolean {
    if (!(factor > 0) || !Number.isFinite(factor) || (anchor !== undefined && !Number.isFinite(anchor))) return false;
    if (this._yOf(section)) {
      const y = this._zoomY(section, factor, anchorY);
      return this._zoomX(section, factor, anchor) || y;
    }
    return this._zoomX(section, factor, anchor);
  }

  private _zoomX(section: string, factor: number, anchor?: number): boolean {
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

  /**
   * Shows the whole axis; a detached follower is attached again; a grid's own y window shows all of y; a grid that
   * follows (its moves leading) fits its leaders, x's and y's.
   */
  fit(section: string): boolean {
    const s = this._sections.get(section);
    if (!s) return false;
    if (s.y) {
      delete s.y;
      this._publish();
    }
    if (s.section.follow) {
      const attached = this.attach(section) || this.isFollowing(section);
      // a grid whose moves lead fits what it shows: its x's leader and its y's (`0` over the plot fits the plot)
      const y = s.section.y;
      if (!y || !this._leading(section)) return attached;
      let moved = this.fit(s.section.follow.leader);
      if (!y.followsX && y.follow) moved = this.fit(y.follow.leader) || moved;
      return moved || attached;
    }
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
   * The last jump (_spec/README.md "ViewModel", History, The last jump): its count (`jumps`), kind, the Location it named
   * and its duration; set before `jumps` moves, so a listener of `jumps` reads it. Undefined before the first.
   */
  get lastJump(): LastJump | undefined {
    return this._lastJump;
  }

  /**
   * Runs a jump (a focus, a fit, a reveal, a bookmark's or a link's move): when it moves the view, the viewports it left
   * are pushed on the history (`record`, the default; a step zoom animates without one) and the steps forward dropped.
   * `same` when the view already showed what was asked (a double-click on what is focused), `none` when nothing could
   * move. A pan or a zoom by hand pushes nothing: a jump after it pushes where the hand left the view.
   */
  jump(move: () => boolean, record = true, kind: JumpKind = 'move', location?: Location): JumpOutcome {
    const before = this._snapshot();
    if (!move()) return 'none';
    const after = this._snapshot();
    if (snapshotsEqual(before, after)) return 'same';
    if (record) this.record(before as unknown as ViewSnapshot);
    this._jumped(kind, location, after);
    return 'moved';
  }

  /** A jump landed (`after`: where it left the view): the last jump told, then the count moved; `after` kept with its Location. */
  private _jumped(kind: JumpKind, location: Location | undefined, after: Snapshot): void {
    const n = this._jumps.value + 1;
    this._landed = location && { ...after, location };
    this._lastJump = { n, kind, location, durationMs: this.config.view.jumpMs };
    this._jumps.set(n);
  }

  /** A step as the history keeps it: where the last jump left the view, with its Location, unless a hand moved it. */
  private _named(snapshot: Snapshot): Snapshot {
    const l = this._landed;
    return l && snapshotsEqual(l, snapshot) ? l : snapshot;
  }

  /**
   * Every group's range and every band's scroll now, as `jump` keeps them (_spec/README.md "ViewModel", History): what
   * a gesture that moves by hand reads at its press, to `record` on its release when it went far (a navigator's drag).
   */
  snapshot(): ViewSnapshot {
    return this._snapshot() as ViewSnapshot;
  }

  /**
   * Pushes a snapshot taken earlier as a step back (the forward steps dropped), unless the view is there now: a gesture's
   * one step, pushed when it ends. False when nothing was pushed.
   */
  record(snapshot: ViewSnapshot): boolean {
    const kept = snapshot as unknown as Snapshot;
    if (snapshotsEqual(kept, this._snapshot())) return false;
    this._steps.push(this._named(kept), this.config.view.historyLimit);
    return true;
  }

  /** What names the bands' rows for the history (the viewer's rows layers); none: a band kept as its `bandTop`. */
  setBandRows(rows: BandRows | undefined): void {
    this._bands.rows = rows;
  }

  /** A rows system as an axis in display order, when a rows layer shows it; undefined else. */
  rowsAxis(system: string): Axis | undefined {
    return this._bands.rows?.axisOf(system);
  }

  /** A band's own viewport (its track's address), undefined when it has none yet. */
  band(track: string): BandViewport | undefined {
    return this._bands.get(track);
  }

  /**
   * Sets a band's own viewport (a rows band's scroll, row zoom, folds), merged over what it holds; undefined fields are
   * left out (`folds: undefined` drops the session's folds: the layer's `collapsed` shows). Not a jump: a hand's move,
   * kept by the next jump's step (a bulk fold runs through `jump`).
   */
  setBand(track: string, band: Partial<BandViewport>): void {
    this._bands.set(track, band);
  }

  /** A band as it was before a gesture (its cancel): its viewport, or none when it had none. */
  restoreBand(track: string, band: BandViewport | undefined): void {
    this._bands.restore(track, band);
  }

  /** The bands of a section's tracks (by track address), as a frame reads them; undefined when none has one. */
  bandsOf(section: string): Readonly<Record<string, BandViewport>> | undefined {
    return this._bands.of(section);
  }

  /** Whether there is a step back or forward. */
  get canBack(): boolean {
    return this._steps.canBack;
  }

  get canForward(): boolean {
    return this._steps.canForward;
  }

  /** Returns to the viewports the last jump left, every leading section's together; false when there is none. */
  back(): boolean {
    return this._step(false);
  }

  /** Goes again where the last back left. */
  forward(): boolean {
    return this._step(true);
  }

  /** Forgets the history (a new spec loaded). */
  clearHistory(): void {
    this._steps.clear();
    this._landed = undefined;
  }

  private _step(forward: boolean): boolean {
    const own = (address: string) => this._ownY(address);
    let location: Location | undefined;
    const moved = this._steps.step(
      forward,
      this._named(this._snapshot()),
      (t) => restoredOf(t, this._groups, this._sections, own),
      (t) => {
        location = t.location;
        restoreInto(t, this._groups, this._sections, this._bands, own);
        this._publish();
      },
    );
    if (moved) this._jumped(forward ? 'forward' : 'back', location, this._snapshot());
    return moved;
  }

  /** Each group's range now, with its axis's signature, its centre as a position and its members. */
  private _snapshot(): Snapshot {
    return snapshotOf(this._groups, this._sections, this._bands);
  }

  // ---- internals ----

  /** A layout's navigation by name (tracks' when it is unknown). */
  private _navOf(layout: string): SectionNav {
    return this._layoutOf(layout).nav;
  }

  /**
   * A layout by name: the layout source's (the frame's, the viewer's registry); tracks' when it is unknown or there is
   * no source, as drawn.
   */
  private _layoutOf(layout: string): SectionLayoutType<any> {
    return this._layouts?.layouts?.[layout] ?? TracksLayout;
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
   * end)` carried back through the follow (`leaderRange`; `ms`, `me`: which way each end moved, `delta`: a pan's). A
   * leader that follows passes the move on up the chain, bounded by `_leading`.
   */
  private _lead(s: SectionState, start: number, end: number, ms: number, me: number, delta?: number): boolean {
    const follow = s.section.follow!;
    const lead = this.get(follow.leader);
    if (!lead || !this.get(s.section.address)) return false;
    const r = leaderRange(s.section, lead, start, end, ms, me, delta, this._arc(), this._leaderOut);
    return !!r && this.setRange(follow.leader, r.start, r.end);
  }

  /** Whether a section has a y that moves apart from x (a grid's own window, or a section it follows). */
  movesY(section: string): boolean {
    return this._yOf(section) !== undefined;
  }

  /**
   * A grid's y that moves apart from x (_spec/README.md "ViewModel", Grids): its own window (`follow_y: none`) or a
   * section it follows (`follow_y: <section>`); undefined when y shows x's window (`follow_y: x`) or there is no y.
   */
  private _yOf(section: string): { s: SectionState; total: number; lead?: SectionState } | undefined {
    const s = this._sections.get(section);
    const y = s?.section.y;
    if (!s || !y || y.followsX) return undefined;
    const total = axisTotal(y.axis);
    if (total === undefined) return undefined;
    // following a section: a state whose follow is y's, so a move asked of y is carried to y's leader (`_lead`)
    const lead = y.follow ? { ...s, section: { ...s.section, axis: y.axis, follow: y.follow } } : undefined;
    return { s, total, ...(lead ? { lead } : {}) };
  }

  /** Whether a section is a grid with a y window of its own (`follow_y: none`). */
  private _ownY(section: string): boolean {
    const r = this._yOf(section);
    return !!r && !r.lead;
  }

  /** A grid's y window now, offsets: its own (all of y by default), or its leader's image. */
  private _yWindow(r: NonNullable<ReturnType<ViewModel['_yOf']>>): { start: number; end: number } | undefined {
    if (!r.lead) return r.s.y ?? { start: 0, end: r.total };
    const follow = r.lead.section.follow!;
    const lead = this.get(follow.leader);
    if (!lead) return undefined;
    const w = followWindow(r.lead.section, follow.leaderAxis, r.lead.section.axis, lead.start, lead.end, this._arc());
    return w && !w.wraps ? w : undefined;
  }

  /** A grid's y set to `[start, end)` (clamped to y), or carried to its leader. */
  private _moveY(r: NonNullable<ReturnType<ViewModel['_yOf']>>, start: number, end: number, delta?: number): boolean {
    const w = this._yWindow(r);
    if (!w) return false;
    if (r.lead) return this._lead(r.lead, start, end, Math.sign(start - w.start), Math.sign(end - w.end), delta);
    const width = Math.min(Math.max(end - start, 1), r.total);
    const from = Math.min(Math.max(start, 0), r.total - width);
    r.s.y = from === 0 && width === r.total ? undefined : { start: from, end: from + width };
    if (!r.s.y) delete r.s.y;
    this._publish();
    return true;
  }

  /** A grid's y panned by `dy` px of its plot's height. */
  private _panY(section: string, dy: number): boolean {
    const r = this._yOf(section);
    const w = r && this._yWindow(r);
    if (!r || !w) return false;
    const height = this._layouts?.placement(section)?.box.height ?? 0;
    if (!(height > 0)) return false;
    // at y's end the pan stops there, its window's width kept (a leader asked past its end would narrow instead)
    const d = stopAtEnds((dy * (w.end - w.start)) / height, w.start, w.end, r.total);
    if (!d) return false;
    return this._moveY(r, w.start + d, w.end + d, d);
  }

  /** A grid's y zoomed by `factor` about `anchor` (an offset down y; the centre when absent). */
  private _zoomY(section: string, factor: number, anchor?: number): boolean {
    const r = this._yOf(section);
    const w = r && this._yWindow(r);
    if (!r || !w) return false;
    const a = anchor !== undefined && Number.isFinite(anchor) ? anchor : (w.start + w.end) / 2;
    return this._moveY(r, a - (a - w.start) / factor, a + (w.end - a) / factor);
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
    const { sections, groups } = groupSections(
      view,
      this._spec.model.value.spec,
      this._sections,
      this._groups,
      this.config.viewport.groupEqualAxes,
    );
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
          // one block at a time: the fit is its block's
          if (this._stopIn(g, true)) g.fitted = false;
          continue;
        }
        // the layout's clamp (a window within its axis and the most px a position; a text's anchor and size)
        if (c) {
          const v = this._navOf(g.layout).clamp(c, rangeOf(g));
          g.start = v.start;
          g.end = v.end;
          if (v.scale !== undefined) g.scale = v.scale;
          if (w) g.end = Math.min(total, g.start + Math.max(1, w.end - w.start));
          this._stopIn(g, false);
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
      this._stopFree = false;
    });
  }

  /**
   * Each group's block under `navigate.acrossBlocks: stop` (the one its window was last kept in), by its group's first section (a group made again over a new block list keeps it): its copy
   * key (system, range, form) and its index then.
   */
  private readonly _stopBlock = new Map<string, { readonly copy: string; readonly index: number }>();
  /** The next publish follows a range set (a jump): under `stop` its window may move to the block holding it. */
  private _stopFree = false;

  /**
   * `navigate.acrossBlocks: stop` on a row of blocks (_spec/README.md "ViewModel", Blocks): the window kept in one
   * block: after a jump the block holding its centre, else the one it was kept in, else its start's; its
   * width at most the block's, its ends inside it (a fit: the whole block). True when it applies.
   */
  private _stopIn(g: Group, fit: boolean): boolean {
    if (this.config.navigate.acrossBlocks !== 'stop') return false;
    const section = this._sections.get(g.members[0])?.section;
    const blocks = section?.blocks;
    if (!blocks || blocks.length < 2 || section.params.arrange === 'stack') return false;
    const at = (u: number) => {
      let lo = 0;
      let hi = blocks.length - 1;
      while (lo < hi) {
        const mid = (lo + hi + 1) >> 1;
        if ((blocks[mid].offset ?? 0) <= u) lo = mid;
        else hi = mid - 1;
      }
      return lo;
    };
    // the kept block found again by its copy key (at its index, else wherever it now is: `section.include` and a
    // reverse change the list under it), never whichever block now sits at its old index
    const k = this._stopBlock.get(g.members[0]);
    let kept = k && blocks[k.index]?.copy === k.copy ? k.index : -1;
    if (k && kept < 0) kept = blocks.findIndex((x) => x.copy === k.copy);
    // a jump (a click on another chromosome, a range set) moves to the block holding its centre; a pan, a zoom and a
    // fit keep the block the window was in; the first time, the block at its start
    const b = !fit && this._stopFree ? at((g.start + g.end) / 2) : kept >= 0 ? kept : at(g.start);
    const lo = blocks[b].offset ?? 0;
    const hi = lo + blockSpan(blocks[b]);
    const width = fit ? hi - lo : Math.min(g.end - g.start, hi - lo);
    const start = fit ? lo : Math.min(Math.max(g.start, lo), hi - width);
    g.start = start;
    g.end = start + width;
    this._stopBlock.set(g.members[0], { copy: blocks[b].copy, index: b });
    return true;
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
      p.y === s.y &&
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
      ...(s.y ? { y: s.y } : {}),
    });
  }

  /**
   * A follower's viewport: its leader's range carried along the follow path (`followerRange`), clamped by
   * `maxPxPerPosition` like a leader. A follower whose leader has no viewport yet keeps its own until the leader's
   * arrives.
   */
  private _follow(s: SectionState): void {
    const follow = s.section.follow!;
    const lead = this.get(follow.leader);
    const address = s.section.address;
    const prev = this.get(address);
    if (!lead) {
      if (prev) this._set(address, prev.start, prev.end, s, true, follow.leader, prev.pieces);
      return;
    }
    const pending = lead.pending || axisTotal(s.section.axis) === undefined || s.width <= 0;
    const r = followerRange(
      s.section,
      lead,
      prev,
      s.width,
      this._arc(),
      this.config.viewport.maxPxPerPosition,
      this._followerOut,
    );
    this._set(address, r.start, r.end, s, pending, follow.leader, r.pieces);
  }
}

function piecesEqual(a: readonly number[] | undefined, b: readonly number[] | undefined): boolean {
  if (a === b) return true;
  if (!a || !b || a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/**
 * A pan by `delta` of window `[start, end)` stopped at the axis's ends `[0, total)` (the window kept whole: never past
 * them); 0 when it cannot move that way. A window that does not fit in them (a follower shown past its own ends) pans
 * as asked.
 */
function stopAtEnds(delta: number, start: number, end: number, total: number): number {
  if (start < -1e-9 || end > total + 1e-9) return delta;
  if (delta > 0) return Math.max(0, Math.min(delta, total - end));
  if (delta < 0) return Math.min(0, Math.max(delta, -start));
  return 0;
}
