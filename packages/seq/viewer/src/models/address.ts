/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Addresses in the viewer (_spec/README.md "AddressModel"): what an address names in the spec shown and the store's
// nodes, the actions on it (focus, reveal, select, highlight on a named channel, hover on the links' own channel), the
// user-defined refs with what each names now (a ref to a node not held waits, and resolves when it arrives), adding one
// as one patch, links in both directions, and the session (the selection, the focus) read and set as addresses. A
// link, a key, a test and an agent reach these through the same commands (the `addressing` behavior).

import { derive, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import { type Issue, issueOf, parseAddress, printAddress, type Segment } from '@molstar/core/repr/tree';
import { addRefPatch, RefsId, type RefInput, refAddresses } from '@molstar/seqx-core/lang';
import {
  type AddressResolution,
  AddressIssues,
  type AddressStatus,
  addressFor,
  EmptyKeys,
  EmptyLocation,
  keySet,
  isEmptyLocation,
  type LinkActionName,
  type Location,
  location,
  locationOfResult,
  parseLink,
  printLink,
  type QueryModel,
  type QueryResult,
  rangesLocation,
  resolveAddressIn,
  unionLocations,
} from '@molstar/seqx-core/query';
import type { Want, Wanter } from '@molstar/seqx-core/data';
import {
  DefaultRenderConfig,
  type HighlightStyle,
  HighlightStyles,
  withSourceNumberings,
} from '@molstar/seqx-core/render';
import { ViewerIssues } from '../issues.js';
import { overlayQueryModel, type RowsSource } from '@molstar/seqx-core/resolve';
import type { Axis } from '@molstar/seqx-core/view';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { breakHolding, nearestOffsetOf, positionAt } from './axis.js';
import type { HoverModel } from './hover.js';
import { type Emphasized, type HighlightModel, isNothing, type SelectionModel } from './selection.js';
import type { EditOptions, EditOutcome, SpecModel } from './spec.js';
import type { StoreModel } from './store.js';
import type { ViewModel } from './view.js';

/**
 * On a cyclic axis over its whole system (offsets the positions), the shortest arc holding sorted, merged ranges: the
 * circle less its widest gap (EC-12), `[start, end]` in offsets with `end` past the length when it wraps the origin;
 * undefined elsewhere (the hull is then the focus).
 */
function shortestArc(axis: Axis, ranges: ArrayLike<number>): [number, number] | undefined {
  const n = axis.length;
  const r = axis.ranges;
  const k = ranges.length;
  if (!axis.cyclic || n === undefined || axis.dir === -1 || r?.length !== 2 || r[0] !== 0 || r[1] !== n || k < 4)
    return undefined;
  let gap = n - ranges[k - 1] + ranges[0];
  let arc: [number, number] | undefined;
  for (let i = 2; i < k; i += 2)
    if (ranges[i] - ranges[i - 1] > gap) {
      gap = ranges[i] - ranges[i - 1];
      arc = [ranges[i], ranges[i - 1] + n];
    }
  return arc;
}

/** A user-defined ref of the spec shown, with what it names now. */
export interface RefEntry {
  /** `/refs/<id>`. */
  readonly address: string;
  readonly id: string;
  /** Its label, else its id. */
  readonly label: string;
  readonly color?: string;
  /** What it names, as written. */
  readonly to: readonly string[];
  readonly status: AddressStatus;
  readonly issues: readonly Issue[];
}

/** What adding a ref did: true when applied, else why not. */
export type AddRefResult = true | readonly Issue[];

/** What an action on an address did. */
export interface AddressOutcome {
  /** Whether it did what was asked (a waiting focus is not done yet: it runs when the address resolves). */
  readonly done: boolean;
  readonly status: AddressStatus;
  readonly message: string;
  readonly issues: readonly Issue[];
}

function sameRefs(a: readonly RefEntry[], b: readonly RefEntry[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    const x = a[i];
    const y = b[i];
    if (x.address !== y.address || x.status !== y.status || x.label !== y.label || x.color !== y.color) return false;
    if (x.to.join(' ') !== y.to.join(' ') || x.issues.length !== y.issues.length) return false;
  }
  return true;
}

/** A focus's outcome when what it names is wholly in a hidden stretch: at its break, and how to show it. */
export const HiddenMessage = 'At its break: hidden by collapsed introns (section.collapse shows them)';

export class AddressModel extends ReactiveModel {
  readonly config: ViewerConfig;
  /** The refs of the spec shown, each with its status; moves with the spec and with the store's nodes. */
  readonly refs: ReadonlyState<readonly RefEntry[]>;
  /** A focus asked of an address that waits: run when it resolves (a ref on a node that has not arrived). */
  private _pendingFocus: { address: string; reveal: boolean; section?: string } | undefined;
  /** The address last focused and the viewports it left, so the focus reads as that address until the view moves. */
  private _focused: { address: string; viewports: string } | undefined;
  /** A focus kept until its sections are placed: records `_focused` when the kept move lands. */
  private _landing: { address: string; stop(): void } | undefined;
  /** The sections the last move moved: what a kept move waits on. */
  private _moved: readonly string[] = [];
  /** The query model last made, and what it was made of (the model, the overlay, the store's version). */
  private _query: { model: unknown; overlay: unknown; version: number; query: QueryModel } | undefined;
  /** The tiles of sources' numberings a form named, asked for through one wanter (core "Numberings from a source"). */
  private _numberingWanter: Wanter | undefined;
  private readonly _numberingWants = new Map<string, Want>();
  private _numberingModel: unknown;

  constructor(
    private readonly spec: SpecModel,
    private readonly view: ViewModel,
    private readonly store: StoreModel,
    private readonly hover: HoverModel,
    readonly selection: SelectionModel,
    readonly highlight: HighlightModel,
    config: ViewerConfig = DefaultViewerConfig,
    /** Positions a source's numbering may cover (`RenderConfig.maxNumberingPositions`). */
    private readonly maxNumberingPositions = DefaultRenderConfig.maxNumberingPositions,
  ) {
    super();
    this.config = config;
    this.refs = derive(
      [spec.model, store.store.nodeVersion],
      () => this._refEntries(),
      // a new list only when an entry changed: what lists the bookmarks redraws once a change
      sameRefs,
    );
    // a waiting focus runs when its address resolves: a node arriving, or the data it reads (refs or not)
    this.listenPermanently(this.refs, () => this._retryFocus());
    this.listenPermanently(store.store.version, () => this._retryFocus());
    this.listenPermanently(store.store.nodeVersion, () => this._retryFocus());
    // the session follows the spec: through a patch's moves before its root is set, then resolved again on the new
    // one, and dropped with an issue when it names nothing
    this.registerPermanentDispose(spec.beforeCommit((changes) => this.followMoves(changes.moved)));
    this.listenPermanently(spec.root, () => this._respec());
    this.registerPermanentDispose(() => this._landing?.stop());
    this.registerPermanentDispose(() => this._numberingWanter?.dispose());
  }

  private readonly _issues = new State<readonly Issue[]>([]);
  /**
   * The session's issues: a selection, a highlight, a focus or a waiting focus dropped because a change of the spec
   * left its address naming nothing (each `seq.address.missing`, naming the address).
   */
  readonly issues: ReadonlyState<readonly Issue[]> = this._issues.asReadonly;

  /**
   * Session addresses after a patch moved or renamed nodes (`Changes.moved` of `repr/patch`): each address whose
   * leading ids name a moved node, or a node below it, is rewritten to the new place, as the spec's refs are.
   */
  followMoves(moved: readonly { readonly from: readonly string[]; readonly to: readonly string[] }[]): void {
    // a gesture's selection (no address) takes its shortest address on the spec it was made on, so the new spec
    // resolves it again (or drops it with an issue) as it does an address's
    let sel = this.selection.value.value;
    if (sel.address === undefined && !isNothing(sel)) {
      const a = this.selectionAddress();
      if (a !== undefined) {
        sel = { ...sel, address: a };
        this.selection.set(sel, this.selection.anchor);
      }
    }
    if (!moved.length) return;
    const rewrite = (address: string) => rewriteMoved(address, moved);
    // the anchor and the cursor name systems: rewritten as the addresses are
    const an = this.selection.anchor;
    const anchor = an ? { ...an, system: rewrite(an.system) } : undefined;
    if (sel.address !== undefined) {
      const a = rewrite(sel.address);
      if (a !== sel.address) this.selection.set({ ...sel, address: a }, anchor);
      else this.selection.set(sel, anchor);
    } else this.selection.set(sel, anchor);
    const c = this.selection.cursor.value;
    if (c) this.selection.setCursor({ ...c, section: rewrite(c.section), system: rewrite(c.system) });
    for (const name of this.highlight.names.value) {
      const h = this.highlight.get(name);
      if (h?.address === undefined) continue;
      const a = rewrite(h.address);
      if (a !== h.address) this.highlight.set(name, { ...h, address: a });
    }
    if (this._pendingFocus)
      this._pendingFocus = { ...this._pendingFocus, address: rewrite(this._pendingFocus.address) };
    if (this._focused) this._focused = { ...this._focused, address: rewrite(this._focused.address) };
    if (this._landing) this._landing.address = rewrite(this._landing.address);
  }

  /** The spec changed: each session address resolved again; one that names nothing now is dropped with an issue. */
  private _respec(): void {
    const dropped: Issue[] = [];
    const check = (address: string): AddressResolution | undefined => {
      const r = this.resolve(address);
      if (r.status === 'resolved' || r.status === 'waiting') return r;
      dropped.push(
        issueOf(AddressIssues, 'seq.address.missing', 'A session address names nothing after the spec changed', {
          address,
          hint: AddressIssues.codes['seq.address.missing'].hint,
        }),
      );
      return undefined;
    };
    const sel = this.selection.value.value;
    if (sel.address !== undefined) {
      if (!check(sel.address)) this.selection.clear();
      else {
        const v = this._emphasized(sel.address).value;
        if (v) this.selection.set(v, this.selection.anchor);
      }
    }
    for (const name of [...this.highlight.names.value]) {
      const h = this.highlight.get(name);
      if (h?.address === undefined) continue;
      if (!check(h.address)) this.highlight.clear(name);
      else {
        const v = this._emphasized(h.address).value;
        if (v) this.highlight.set(name, h.style ? { ...v, style: h.style } : v);
      }
    }
    if (this._pendingFocus && !check(this._pendingFocus.address)) this._pendingFocus = undefined;
    if (this._focused && !check(this._focused.address)) this._focused = undefined;
    if (this._landing && !check(this._landing.address)) this._landing.stop();
    this._issues.set(dropped);
  }

  // ---- resolving ----

  /**
   * The query model addresses resolve over: the static model of what is shown, the store's nodes and held rows now,
   * and the labels of the numberings read from sources that the store holds.
   */
  queryModel(): QueryModel {
    const model = this.spec.model.value;
    const overlay = this.store.currentOverlay(model);
    const version = this.store.store.version.value;
    const q = this._query;
    if (q && q.model === model && q.overlay === overlay && q.version === version) return q.query;
    // one rows source a model (the query models are kept by it): rows read from a source asked for as numberings are
    if (this._rows?.model !== model)
      this._rows = { model, source: { data: this.store.store, want: (w) => this._wantTiles(model, w) } };
    const query = withSourceNumberings(overlayQueryModel(model, overlay, this._rows.source), model, this.store.store, {
      maxPositions: this.maxNumberingPositions,
      want: (w) => this._wantTiles(model, w),
    });
    this._query = { model, overlay, version, query };
    return query;
  }

  /** The rows source of the model shown (`queryModel`). */
  private _rows: { model: unknown; source: RowsSource } | undefined;

  /**
   * Asks for the tiles a form named, a source numbering's labels or an alignment's rows (whoever resolves asks, so a
   * link in a numbering no ruler shows, or to a row no layer read, resolves): the wants of the model shown, all at
   * once, after the resolve that named them (never inside it).
   */
  private _wantTiles(model: unknown, w: Want) {
    if (this._numberingModel !== model) {
      this._numberingModel = model;
      this._numberingWants.clear();
    }
    const key = `${w.source}|${w.system}|${w.end}`;
    if (this._numberingWants.has(key)) return;
    this._numberingWants.set(key, w);
    const wants = [...this._numberingWants.values()];
    queueMicrotask(() => {
      if (this._numberingModel !== model) return;
      this._numberingWanter ??= this.store.store.wanter();
      this._numberingWanter.want(wants);
    });
  }

  /** What an address names in the spec shown (_spec/README.md "AddressModel"; the core's `resolveAddressIn`). */
  resolve(address: string): AddressResolution {
    // no spec set yet: no address names anything, but one will (a link followed before the document arrives waits)
    if (!this.spec.loaded.value && parseAddress(address).pattern) {
      const waiting = issueOf(AddressIssues, 'seq.address.waiting', 'The address waits for its data', { address });
      return { address, status: 'waiting', via: [], issues: [waiting], complete: false, pending: EmptyLocation };
    }
    const { limit, maxDepth, refDepth } = this.config.address;
    return resolveAddressIn(this.spec.shown.value, this.queryModel(), address, { limit, maxDepth, refDepth });
  }

  /** What an address names as positions (elements by their extents, a system whole), with its resolution. */
  locate(address: string): { readonly resolution: AddressResolution; readonly location: Location } {
    const resolution = this.resolve(address);
    if (!resolution.result) return { resolution, location: EmptyLocation };
    return {
      resolution,
      location: locationOfResult(this.spec.shown.value, this.queryModel(), resolution.result).location,
    };
  }

  /** The shortest stable address of a result: a ref naming exactly it, else ids and forms (`copyLink` settings). */
  addressOf(result: QueryResult): string | undefined {
    const { prefer, numbering } = this.config.copyLink;
    return addressFor(this.spec.shown.value, this.queryModel(), result, { prefer, numbering, ...this.config.address });
  }

  private _refEntries(): RefEntry[] {
    const root = this.spec.shown.value;
    const container = root.children?.find((c) => c.id === RefsId && c.kind === 'refs');
    const out: RefEntry[] = [];
    for (const r of container?.children ?? []) {
      if (r.kind !== 'ref' || r.id === undefined) continue;
      const address = printAddress([RefsId, r.id]);
      const res = this.resolve(address);
      const label = r.params?.label;
      const color = r.params?.color;
      out.push({
        address,
        id: r.id,
        label: typeof label === 'string' ? label : r.id,
        ...(typeof color === 'string' ? { color } : {}),
        to: refAddresses(r),
        status: res.status,
        issues: res.issues,
      });
    }
    return out;
  }

  // ---- refs as patches ----

  /**
   * Adds a ref (a bookmark) as one patch over the root (`addRefPatch`, through `SpecModel.patch`: a ref that would
   * dangle is refused), one step of the history when one is bound. True when applied, else the issues; a promise when
   * the history runs it.
   */
  addRef(ref: RefInput, options: EditOptions = {}): AddRefResult | Promise<AddRefResult> {
    // made from the root the edit runs on: bookmarks added one after another each apply
    const outcome = this.spec.patch('Add a bookmark', (root) => addRefPatch(root, ref), options);
    const result = (o: EditOutcome): AddRefResult => (o.ok ? true : o.issues);
    return outcome instanceof Promise ? outcome.then(result) : result(outcome);
  }

  /**
   * Bookmarks the selection as one patch (_spec/README.md "AddressModel"): a ref to the selection's address, labeled by
   * it (what a link would say) unless a label is given, under the id given or the next free `bookmark-<n>` of the root
   * the edit runs on. Undefined when nothing is selected, else as `addRef`.
   */
  bookmarkSelection(
    options: EditOptions & { readonly id?: string; readonly label?: string } = {},
  ): AddRefResult | Promise<AddRefResult> | undefined {
    const to = this.selectionAddress();
    if (to === undefined) return undefined;
    const { id, label, ...edit } = options;
    const outcome = this.spec.patch(
      'Add a bookmark',
      (root) => {
        const taken = new Set(root.children?.find((c) => c.id === RefsId)?.children?.map((c) => c.id));
        let n = 1;
        while (id === undefined && taken.has(`bookmark-${n}`)) n++;
        return addRefPatch(root, { id: id ?? `bookmark-${n}`, to, label: label ?? to });
      },
      edit,
    );
    const result = (o: EditOutcome): AddRefResult => (o.ok ? true : o.issues);
    return outcome instanceof Promise ? outcome.then(result) : result(outcome);
  }

  // ---- actions ----

  private _outcome(resolution: AddressResolution, done: boolean, message: string): AddressOutcome {
    return { done, status: resolution.status, message, issues: resolution.issues };
  }

  /** Whether the last move found what it moved to wholly in a hidden stretch (and moved to its break). */
  private _hidden = false;

  /** Whether the last focus or reveal moved to a break: what it names is hidden by collapsed pieces. */
  get movedToBreak(): boolean {
    return this._hidden;
  }

  /** Moves the view to what an address names: each leading section showing it (or a path to it), its whole extent. */
  focus(address: string, options: { section?: string } = {}): AddressOutcome {
    return this._move(address, false, options.section);
  }

  /** Brings what an address names into view, panning only when it is not shown; zooming out only when it must. */
  reveal(address: string, options: { section?: string } = {}): AddressOutcome {
    return this._move(address, true, options.section);
  }

  private _move(address: string, reveal: boolean, section?: string): AddressOutcome {
    const { resolution, location: loc } = this.locate(address);
    if (resolution.status === 'waiting') {
      this._pendingFocus = { address, reveal, ...(section !== undefined ? { section } : {}) };
      return this._outcome(resolution, false, 'Waiting for its data');
    }
    if (this._pendingFocus?.address !== address) this._pendingFocus = undefined;
    this._landing?.stop();
    if (resolution.status !== 'resolved') return this._outcome(resolution, false, 'Nothing to focus');
    this._pendingFocus = undefined;
    const r = this._moveTo(loc, reveal, section);
    // a section not placed yet (no width, an axis waiting for its length) keeps the move as its intent: not done yet
    if (r === 'pending') {
      this._keepUntilPlaced(address);
      return this._outcome(resolution, false, 'Kept until the view is placed');
    }
    if (r === 'moved') this._focused = { address, viewports: this._viewportsKey() };
    return this._outcome(
      resolution,
      r === 'moved',
      r === 'moved' ? (this._hidden ? HiddenMessage : reveal ? 'Revealed' : 'Focused') : 'Not in this view',
    );
  }

  /**
   * Records the focus as the address once the sections the kept move moved are placed (a width, a length): not before,
   * and never after another jump (a focus, back, forward) or a change of the spec that leaves the address naming
   * nothing. A section that stays pending without being moved (hidden) does not hold it.
   */
  private _keepUntilPlaced(address: string): void {
    const waiting = this._moved.filter((a) => this.view.get(a)?.pending);
    const offs: (() => void)[] = [];
    const landing = {
      address,
      stop: () => {
        for (const off of offs.splice(0)) off();
        if (this._landing === landing) this._landing = undefined;
      },
    };
    const land = () => {
      if (waiting.some((a) => this.view.get(a)?.pending)) return;
      landing.stop();
      this._focused = { address: landing.address, viewports: this._viewportsKey() };
    };
    this._landing = landing;
    for (const a of waiting) offs.push(this.view.viewport(a).listen(land));
    // a jump that is not this move's (the view's back and forward, another focus) takes the record away
    // (told after this move's own jump, so the count it started at is compared)
    const jumps = this.view.jumps.value;
    offs.push(
      this.view.jumps.listen((n) => {
        if (n !== jumps) landing.stop();
      }),
    );
  }

  private _retryFocus(): void {
    const p = this._pendingFocus;
    if (!p || this.resolve(p.address).status === 'waiting') return;
    this._move(p.address, p.reveal, p.section);
  }

  /**
   * Moves the view to a Location as `focus` moves it to an address's: each leading, movable section showing it (or a
   * static path to it) shows its hull; followers follow. What the cursor's Enter runs.
   */
  focusLocation(loc: Location): 'moved' | 'pending' | 'none' {
    return this._moveTo(loc, false);
  }

  /** Brings a Location into view as `reveal` brings an address's (_spec/README.md "AddressModel", Actions). */
  revealLocation(loc: Location): 'moved' | 'pending' | 'none' {
    return this._moveTo(loc, true);
  }

  /**
   * A focus or a reveal is a jump (_spec/README.md "ViewModel", History): the viewports it leaves are a step of the
   * view's history.
   */
  private _moveTo(loc: Location, reveal: boolean, only?: string): 'moved' | 'pending' | 'none' {
    let r: 'moved' | 'pending' | 'none' = 'none';
    this._landing?.stop();
    this.view.jump(
      () => {
        r = this._moveEach(loc, reveal, only);
        return r !== 'none';
      },
      true,
      reveal ? 'reveal' : 'focus',
      loc,
    );
    return r;
  }

  /** Each leading section on a system the location has positions on (or reaches by a static path) shows its hull. */
  private _moveEach(loc: Location, reveal: boolean, only?: string): 'moved' | 'pending' | 'none' {
    const resolved = this.view.resolved.value;
    this._hidden = false;
    if (!resolved || isEmptyLocation(loc)) return 'none';
    let hidden = false;
    // a rect's sides are positions on its two systems
    let positions = location({ parts: loc.parts });
    for (const r of loc.rects) positions = unionLocations(positions, location({ parts: [r.a, r.b] }));
    const model = this.spec.model.value;
    let moved = false;
    const done = new Set<string>();
    for (const s of resolved.sections) {
      if (only !== undefined && s.address !== only) continue;
      if (this.view.isFollowing(s.address) || !this.view.canMove(s.address)) continue;
      // one move a scroll group
      const group = this.view.groupOf(s.address);
      if (group.some((g) => done.has(g))) continue;
      let axis = s.axis;
      let part = positions.parts.find((p) => p.system === axis.system);
      // a section of blocks on several systems (its axis its own units, on no system): the first block holding some
      // of the positions, on its own axis from its offset
      let base = 0;
      if (s.blocks && !model.systems.has(axis.system)) {
        for (const b of s.blocks) {
          const q = positions.parts.find((p) => p.system === b.axis.system);
          if (!q?.ranges.length || !b.axis.ranges) continue;
          [axis, part, base] = [b.axis, q, b.offset ?? 0];
          break;
        }
      }
      if (!part) {
        const projected = model.systems.has(axis.system)
          ? model.project(positions, axis.system, { via: s.via }).location.parts.find((p) => p.system === axis.system)
          : undefined;
        part = projected;
      }
      if (!part?.ranges.length) continue;
      // the hull in display offsets (_spec/README.md "AddressModel", Focus): hidden stretches take no width, so what
      // lies wholly in one is its break, zero wide; an axis waiting for its length reads offsets as positions from 0
      // (the intent is kept, EC-112)
      let a = Infinity;
      let b = -Infinity;
      let wide = false;
      for (let i = 0; i + 1 < part.ranges.length; i += 2) {
        const x = axis.ranges ? (nearestOffsetOf(axis, part.ranges[i]) ?? part.ranges[i]) : part.ranges[i];
        const y = axis.ranges ? (nearestOffsetOf(axis, part.ranges[i + 1]) ?? part.ranges[i + 1]) : part.ranges[i + 1];
        a = Math.min(a, x + base, y + base);
        b = Math.max(b, x + base, y + base);
        if (x !== y) wide = true;
      }
      // on a ring's whole system, an item across the origin is its shortest arc, not the hull of its runs (EC-12)
      const arc = base === 0 ? shortestArc(axis, part.ranges) : undefined;
      if (arc) [a, b] = arc;
      const v = this.view.get(s.address);
      // wholly hidden by a break: to its break, centered, the zoom kept; outside the section's range is not hidden by a
      // break (the focus as before: what the view can show of it)
      if (!wide && axis.ranges && breakHolding(axis, part.ranges[0]) >= 0) {
        hidden = true;
        moved = (v && !v.pending ? this.view.pan(s.address, a - (v.start + v.end) / 2) : false) || moved;
        done.add(s.address);
        continue;
      }
      if (reveal && v && !v.pending) {
        if (a >= v.start && b <= v.end) {
          moved = true;
          done.add(s.address);
          continue;
        }
        if (b - a <= v.end - v.start) {
          moved = this.view.pan(s.address, (a + b) / 2 - (v.start + v.end) / 2) || moved;
          done.add(s.address);
          continue;
        }
      }
      moved = this.view.setRange(s.address, Math.min(a, b), Math.max(a, b)) || moved;
      done.add(s.address);
    }
    this._moved = [...done];
    this._hidden = hidden;
    if (!moved) return 'none';
    for (const a of done) if (this.view.get(a)?.pending) return 'pending';
    return 'moved';
  }

  /** What an address names as an emphasis: its elements and positions, and what is drawn of them. */
  private _emphasized(address: string): { resolution: AddressResolution; value?: Emphasized } {
    const resolution = this.resolve(address);
    const r = resolution.result;
    if (!r) return { resolution };
    const shown = locationOfResult(this.spec.shown.value, this.queryModel(), r).location;
    return {
      resolution,
      value: {
        address,
        elements: r.type === 'keys' ? r.keys : EmptyKeys,
        location: r.type === 'location' ? r.location : r.type === 'nodes' ? shown : EmptyLocation,
        shown,
      },
    };
  }

  /** Sets the selection to what an address names (session state). */
  select(address: string): AddressOutcome {
    const { resolution, value } = this._emphasized(address);
    if (!value) return this._outcome(resolution, false, 'Nothing to select');
    this.selection.set(value);
    return this._outcome(resolution, true, 'Selected');
  }

  /** Highlights what an address names on a channel (the links' channel by default). */
  highlightAddress(
    address: string,
    channel = this.config.links.highlightChannel,
    style?: HighlightStyle,
  ): AddressOutcome {
    const { resolution, value } = this._emphasized(address);
    if (!value) return this._outcome(resolution, false, 'Nothing to highlight');
    if (style !== undefined && !HighlightStyles.includes(style)) {
      this._badStyle();
      style = 'tint';
    }
    this.highlight.set(channel, style ? { ...value, style } : value);
    return this._outcome(resolution, true, 'Highlighted');
  }

  /** An unknown style is a tint, reported once: the issue stays while the session's issues do. */
  private _badStyle(): void {
    if (this._issues.value.some((i) => i.code === 'seq.viewer.highlight.style')) return;
    this._issues.set([
      ...this._issues.value,
      issueOf(ViewerIssues, 'seq.viewer.highlight.style', 'A highlight style is unknown: drawn as a tint', {
        hint: ViewerIssues.codes['seq.viewer.highlight.style'].hint,
      }),
    ]);
  }

  /** Clears a highlight channel, or every one: `highlightAddress`'s inverse. The channels cleared. */
  clearHighlight(channel?: string): number {
    const names =
      channel === undefined ? [...this.highlight.names.value] : this.highlight.get(channel) ? [channel] : [];
    for (const n of names) this.highlight.clear(n);
    return names.length;
  }

  /** A hover over a link: its own channel (EC-182), so the pointer's hover is left alone. */
  hoverAddress(address: string): AddressOutcome {
    return this.highlightAddress(address, this.config.links.hoverChannel);
  }

  /** Ends a link's hover. */
  unhover(): void {
    this.highlight.clear(this.config.links.hoverChannel);
  }

  /** Runs an action by its link name. */
  run(action: LinkActionName, address: string, channel?: string): AddressOutcome {
    switch (action) {
      case 'focus':
        return this.focus(address);
      case 'reveal':
        return this.reveal(address);
      case 'select':
        return this.select(address);
      case 'highlight':
        return this.highlightAddress(address, channel);
      case 'hover':
        return this.hoverAddress(address);
    }
  }

  // ---- links ----

  /** A link of one action on an address, in this viewer's prefix. */
  link(action: LinkActionName, address: string, channel?: string): string {
    return printLink(
      { actions: [{ action, address }], ...(channel !== undefined ? { channel } : {}) },
      { prefix: this.config.links.prefix },
    );
  }

  /** Parses a link (undefined when the text is none) with this viewer's prefix and bounds. */
  parseLink(text: string) {
    const { prefix, maxLength, maxActions } = this.config.links;
    return parseLink(text, { prefix, maxLength, maxActions });
  }

  /**
   * The link of what is selected, else of what the pointer is over, as the shortest stable address (`copyLink`):
   * undefined when there is nothing.
   */
  copyLink(): string | undefined {
    const a = this.selectionAddress() ?? this.hoverAddressNow();
    return a === undefined ? undefined : this.link(this.config.copyLink.action, a);
  }

  /** The address of what the pointer is over: its element, else its position. */
  hoverAddressNow(): string | undefined {
    const h = this.hover.current;
    if (h.owner === undefined || h.system === undefined) return undefined;
    if (h.element) return this.addressOf({ type: 'keys', keys: keySet([h.element]) });
    if (Number.isNaN(h.position)) return undefined;
    const p = Math.floor(h.position);
    return this.addressOf({ type: 'location', location: rangesLocation(h.system, [p, p + 1]) });
  }

  // ---- the session as addresses ----

  /** The selection's address: the address it was set from, else the shortest of what it holds. */
  selectionAddress(): string | undefined {
    const s = this.selection.value.value;
    if (isNothing(s)) return undefined;
    if (s.address !== undefined) return s.address;
    return s.elements.keys.length
      ? this.addressOf({ type: 'keys', keys: s.elements })
      : this.addressOf({ type: 'location', location: s.location });
  }

  /** Sets the selection from an address (none: clears it). */
  setSelection(address: string | undefined): AddressOutcome | undefined {
    if (address === undefined) {
      this.selection.clear();
      return undefined;
    }
    return this.select(address);
  }

  /**
   * The focus as an address: the address last focused while the view has not moved since, else the region the first
   * leading section shows (positions, in its system's ordinal numbering).
   */
  focusAddress(): string | undefined {
    if (this._focused && this._focused.viewports === this._viewportsKey()) return this._focused.address;
    return this.shownAddress();
  }

  /**
   * The region the first leading section shows, as an address of ids and forms (positions in `copyLink.numbering`):
   * what a bookmark of the view stores. Undefined while the view is not placed.
   */
  shownAddress(): string | undefined {
    const resolved = this.view.resolved.value;
    for (const s of resolved?.sections ?? []) {
      if (this.view.isFollowing(s.address)) continue;
      const v = this.view.get(s.address);
      if (!v || v.pending) return undefined;
      // a section of blocks on several systems shows no one system's positions
      if (!this.spec.model.value.systems.has(s.axis.system)) continue;
      const a = positionAt(s.axis, v.start, 'start');
      const b = positionAt(s.axis, v.end, 'end');
      if (a === undefined || b === undefined) return undefined;
      // positions, lowest first: a reversed axis's window starts at its highest
      const from = Math.floor(Math.min(a, b));
      const to = Math.max(from + 1, Math.ceil(Math.max(a, b)));
      // ids and forms only: a bookmark of the view never names another ref, so bookmarks never chain
      return addressFor(
        this.spec.shown.value,
        this.queryModel(),
        { type: 'location', location: rangesLocation(s.axis.system, [from, to]) },
        { ...this.config.address, prefer: ['ids'], numbering: this.config.copyLink.numbering },
      );
    }
    return undefined;
  }

  /** Sets the focus from an address: the `focus` action. */
  setFocus(address: string): AddressOutcome {
    return this.focus(address);
  }

  private _viewportsKey(): string {
    let key = '';
    for (const s of this.view.resolved.value?.sections ?? []) {
      const v = this.view.get(s.address);
      key += v ? `${s.address}:${v.start},${v.end};` : `${s.address};`;
    }
    return key;
  }
}

/** An address with its leading ids rewritten through moves (a node moved, or one below it), forms after them kept. */
function rewriteMoved(
  address: string,
  moved: readonly { readonly from: readonly string[]; readonly to: readonly string[] }[],
): string {
  const p = parseAddress(address).pattern;
  if (!p) return address;
  let segments: readonly Segment[] = p.segments;
  for (const m of moved) {
    const n = m.from.length;
    if (!n || segments.length < n) continue;
    let prefix = true;
    for (let i = 0; i < n && prefix; i++) prefix = segments[i] === m.from[i];
    if (prefix) segments = [...m.to, ...segments.slice(n)];
  }
  return segments === p.segments ? address : printAddress(segments, p.subtree);
}
