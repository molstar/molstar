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
  EmptyElements,
  EmptyLocation,
  elements,
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
import { overlayQueryModel } from '@molstar/seqx-core/resolve';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { nearestOffsetOf, positionAt } from './axis.js';
import type { HoverModel } from './hover.js';
import { type Emphasized, type HighlightModel, isNothing, type SelectionModel } from './selection.js';
import type { EditOptions, EditOutcome, SpecModel } from './spec.js';
import type { StoreModel } from './store.js';
import type { ViewModel } from './view.js';

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

export class AddressModel extends ReactiveModel {
  readonly config: ViewerConfig;
  /** The refs of the spec shown, each with its status; moves with the spec and with the store's nodes. */
  readonly refs: ReadonlyState<readonly RefEntry[]>;
  /** A focus asked of an address that waits: run when it resolves (a ref on a node that has not arrived). */
  private _pendingFocus: { address: string; reveal: boolean; section?: string } | undefined;
  /** The address last focused and the viewports it left, so the focus reads as that address until the view moves. */
  private _focused: { address: string; viewports: string } | undefined;

  constructor(
    private readonly spec: SpecModel,
    private readonly view: ViewModel,
    private readonly store: StoreModel,
    private readonly hover: HoverModel,
    readonly selection: SelectionModel,
    readonly highlight: HighlightModel,
    config: ViewerConfig = DefaultViewerConfig,
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
    if (!moved.length) return;
    const rewrite = (address: string) => rewriteMoved(address, moved);
    const sel = this.selection.value.value;
    if (sel.address !== undefined) {
      const a = rewrite(sel.address);
      if (a !== sel.address) this.selection.set({ ...sel, address: a });
    }
    for (const name of this.highlight.names.value) {
      const h = this.highlight.get(name);
      if (h?.address === undefined) continue;
      const a = rewrite(h.address);
      if (a !== h.address) this.highlight.set(name, { ...h, address: a });
    }
    if (this._pendingFocus)
      this._pendingFocus = { ...this._pendingFocus, address: rewrite(this._pendingFocus.address) };
    if (this._focused) this._focused = { ...this._focused, address: rewrite(this._focused.address) };
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
        if (v) this.selection.set(v);
      }
    }
    for (const name of [...this.highlight.names.value]) {
      const h = this.highlight.get(name);
      if (h?.address === undefined) continue;
      if (!check(h.address)) this.highlight.clear(name);
      else {
        const v = this._emphasized(h.address).value;
        if (v) this.highlight.set(name, v);
      }
    }
    if (this._pendingFocus && !check(this._pendingFocus.address)) this._pendingFocus = undefined;
    if (this._focused && !check(this._focused.address)) this._focused = undefined;
    this._issues.set(dropped);
  }

  // ---- resolving ----

  /** The query model addresses resolve over: the static model of what is shown and the store's nodes now. */
  queryModel(): QueryModel {
    const model = this.spec.model.value;
    return overlayQueryModel(model, this.store.currentOverlay(model));
  }

  /** What an address names in the spec shown (_spec/README.md "AddressModel"; the core's `resolveAddressIn`). */
  resolve(address: string): AddressResolution {
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

  // ---- actions ----

  private _outcome(resolution: AddressResolution, done: boolean, message: string): AddressOutcome {
    return { done, status: resolution.status, message, issues: resolution.issues };
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
    if (resolution.status !== 'resolved') return this._outcome(resolution, false, 'Nothing to focus');
    this._pendingFocus = undefined;
    const r = this._moveTo(loc, reveal, section);
    // a section not placed yet (no width, an axis waiting for its length) keeps the move as its intent: not done yet
    if (r === 'pending') return this._outcome(resolution, false, 'Kept until the view is placed');
    if (r === 'moved') this._focused = { address, viewports: this._viewportsKey() };
    return this._outcome(
      resolution,
      r === 'moved',
      r === 'moved' ? (reveal ? 'Revealed' : 'Focused') : 'Not in this view',
    );
  }

  private _retryFocus(): void {
    const p = this._pendingFocus;
    if (!p || this.resolve(p.address).status === 'waiting') return;
    this._move(p.address, p.reveal, p.section);
  }

  /** Each leading section on a system the location has positions on (or reaches by a static path) shows its hull. */
  private _moveTo(loc: Location, reveal: boolean, only?: string): 'moved' | 'pending' | 'none' {
    const resolved = this.view.resolved.value;
    if (!resolved || isEmptyLocation(loc)) return 'none';
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
      const axis = s.axis;
      let part = positions.parts.find((p) => p.system === axis.system);
      if (!part) {
        const projected = model.systems.has(axis.system)
          ? model.project(positions, axis.system).location.parts.find((p) => p.system === axis.system)
          : undefined;
        part = projected;
      }
      if (!part?.ranges.length) continue;
      const first = part.ranges[0];
      const last = part.ranges[part.ranges.length - 1];
      // an axis waiting for its length reads offsets as positions from 0 (the intent is kept, EC-112)
      const a = axis.ranges ? (nearestOffsetOf(axis, first) ?? first) : first;
      const b = axis.ranges ? (nearestOffsetOf(axis, last) ?? last) : last;
      const v = this.view.get(s.address);
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
        elements: r.type === 'elements' ? r.elements : EmptyElements,
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
  highlightAddress(address: string, channel = this.config.links.highlightChannel): AddressOutcome {
    const { resolution, value } = this._emphasized(address);
    if (!value) return this._outcome(resolution, false, 'Nothing to highlight');
    this.highlight.set(channel, value);
    return this._outcome(resolution, true, 'Highlighted');
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
    if (h.element) return this.addressOf({ type: 'elements', elements: elements([h.element]) });
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
      ? this.addressOf({ type: 'elements', elements: s.elements })
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
      const start = positionAt(s.axis, v.start);
      const end = positionAt(s.axis, v.end);
      if (start === undefined || end === undefined) return undefined;
      const from = Math.floor(start);
      const to = Math.max(from + 1, Math.ceil(end));
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
