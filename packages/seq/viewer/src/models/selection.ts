/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The selection and the highlight channels (_spec/README.md "SelectionModel and HighlightModel"): session state, as
// addresses and Locations, never in the spec. The selection holds element keys and a Location (what is drawn: the
// Location and the elements' extents), the anchor an extend starts from, and the keyboard cursor; the grammar's gestures
// (the interaction grammar, by the hit: an item is an element, a cell a position) are methods here, so a binding, a
// key and a test act alike.

import { KeyedState, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import {
  type ElementKey,
  type Elements,
  EmptyElements,
  EmptyLocation,
  elements,
  elementsEqual,
  isEmptyLocation,
  type Location,
  locationsEqual,
  partOf,
  rangesLocation,
  subtractLocations,
  unionLocations,
} from '@molstar/seqx-core/query';
import type { Hit } from './hover.js';

/** What is selected or highlighted: element keys, positions, and the positions drawn (theirs and the elements' extents). */
export interface Emphasized {
  readonly elements: Elements;
  readonly location: Location;
  /** What is drawn: the location and the elements' extents. */
  readonly shown: Location;
  /** The address it was set from, if it was. */
  readonly address?: string;
}

export const Nothing: Emphasized = Object.freeze({
  elements: EmptyElements,
  location: EmptyLocation,
  shown: EmptyLocation,
});

const same = (a: Emphasized, b: Emphasized) =>
  a === b ||
  (a.address === b.address &&
    elementsEqual(a.elements, b.elements) &&
    locationsEqual(a.location, b.location) &&
    locationsEqual(a.shown, b.shown));

export function isNothing(e: Emphasized): boolean {
  return !e.elements.keys.length && !e.location.parts.length && !e.location.rects.length;
}

/** Where an extend starts: a range of positions on one system (a cell's, or an item's extent there). */
export interface Anchor {
  readonly system: string;
  readonly start: number;
  readonly end: number;
}

/** The keyboard cursor: a position on a section's system (drawn on every section of that system). */
export interface Cursor {
  readonly section: string;
  readonly system: string;
  readonly position: number;
}

const sameCursor = (a: Cursor | undefined, b: Cursor | undefined) =>
  a === b || (!!a && !!b && a.section === b.section && a.system === b.system && a.position === b.position);

const keyId = (k: ElementKey) => `${typeof k.key === 'number' ? 'n' : 's'}${k.key}\u0000${k.owner}`;

/** The positions a hit covers on its own system: an item's extent there, else its cell; undefined for neither. */
export function hitRange(hit: Hit): Anchor | undefined {
  const part = hit.extent ? partOf(hit.extent, hit.system) : undefined;
  if (part?.ranges.length)
    return { system: hit.system, start: part.ranges[0], end: part.ranges[part.ranges.length - 1] };
  if (hit.position === undefined || Number.isNaN(hit.position)) return undefined;
  const p = Math.floor(hit.position);
  return { system: hit.system, start: p, end: p + 1 };
}

export class SelectionModel extends ReactiveModel {
  private readonly _value = new State<Emphasized>(Nothing, same);
  private readonly _cursor = new State<Cursor | undefined>(undefined, sameCursor);
  /** The extents of elements selected by a hit, by key: what is drawn of them without resolving them again. */
  private readonly _extents = new Map<string, Location>();
  private _anchor: Anchor | undefined;
  private _drag: { anchor: Anchor; snap: boolean } | undefined;
  /** The selection: an equal value tells nothing. */
  readonly value: ReadonlyState<Emphasized> = this._value.asReadonly;
  /** The keyboard cursor; an equal one tells nothing. */
  readonly cursor: ReadonlyState<Cursor | undefined> = this._cursor.asReadonly;

  /**
   * `extentOf`: an element's extent when no hit gave it (one selected by an address, then toggled): the viewer's
   * addresses resolve it. Without one, an element with no known extent is held and not drawn.
   */
  constructor(private readonly extentOf?: (elements: Elements) => Location) {
    super();
  }

  /** Where an extend starts: the last element or cell set, toggled or begun from (none after a clear). */
  get anchor(): Anchor | undefined {
    return this._anchor;
  }

  set(value: Emphasized, anchor?: Anchor): void {
    this._anchor = anchor;
    this._value.set(value);
  }

  clear(): void {
    this._anchor = undefined;
    this._drag = undefined;
    this._value.set(Nothing);
  }

  // ---- the grammar's gestures, by the hit (the interaction grammar) ----

  /**
   * A click: the selection becomes what was hit (an item: its key; a cell: its position); clicking the only selected
   * element clears it. Nothing hit clears it. True when the selection changed.
   */
  select(hit: Hit | undefined): boolean {
    const before = this._value.value;
    const range = hit ? hitRange(hit) : undefined;
    if (!hit || (!hit.element && !range)) {
      this.clear();
      return before !== this._value.value;
    }
    const next = this._ofHit(hit, range);
    if (same(next, before)) this.clear();
    else this.set(next, range);
    return before !== this._value.value;
  }

  /**
   * Shift+click: the positions from the anchor to what was hit, on the anchor's system (a range in one system); with
   * no anchor there, a click. The anchor stays, so extends from it replace one another.
   */
  extend(hit: Hit | undefined): boolean {
    const a = this._anchor;
    const range = hit ? hitRange(hit) : undefined;
    if (!a || !range || range.system !== a.system) return hit ? this.select(hit) : false;
    const before = this._value.value;
    const loc = rangesLocation(a.system, [Math.min(a.start, range.start), Math.max(a.end, range.end)]);
    this._value.set({ elements: EmptyElements, location: loc, shown: loc });
    return before !== this._value.value;
  }

  /** Cmd/Ctrl+click: what was hit in or out of the selection, the rest kept; it becomes the anchor. */
  toggle(hit: Hit | undefined): boolean {
    const range = hit ? hitRange(hit) : undefined;
    if (!hit || (!hit.element && !range)) return false;
    const before = this._value.value;
    const cur = before.address === undefined ? before : { ...before, address: undefined };
    let els = cur.elements;
    let loc = cur.location;
    if (hit.element) {
      const has = els.keys.some((k) => k.owner === hit.element!.owner && k.key === hit.element!.key);
      if (hit.extent) this._extents.set(keyId(hit.element), hit.extent);
      els = elements(
        has
          ? els.keys.filter((k) => k.owner !== hit.element!.owner || k.key !== hit.element!.key)
          : [...els.keys, hit.element],
      );
    } else if (range) {
      const cell = rangesLocation(range.system, [range.start, range.end]);
      const inside = locationsEqual(subtractLocations(cell, loc), EmptyLocation);
      loc = inside ? subtractLocations(loc, cell) : unionLocations(loc, cell);
    }
    this._anchor = range;
    this._value.set(this._compose(els, loc));
    return before !== this._value.value;
  }

  /**
   * A drag-select begins on what was hit: the range grows from it as the pointer moves (`dragTo`), snapped to whole
   * items when it began on one (the interaction grammar). False when nothing selectable was hit.
   */
  beginDrag(hit: Hit | undefined): boolean {
    const range = hit ? hitRange(hit) : undefined;
    if (!range) return false;
    this._drag = { anchor: range, snap: !!hit?.element };
    return true;
  }

  /**
   * The drag's other end: a position on the drag's system (the pointer's, read from the view even past the data), or
   * an item's extent when the drag began on items and the pointer is over one. The selection becomes the range between
   * the two ends; it is the anchor's for a later extend. True when the selection changed.
   */
  dragTo(position: number, hit?: Hit): boolean {
    const d = this._drag;
    if (!d || Number.isNaN(position)) return false;
    const a = d.anchor;
    let start = Math.floor(position);
    let end = start + 1;
    if (d.snap && hit?.element && hit.system === a.system) {
      const r = hitRange(hit);
      if (r) {
        start = r.start;
        end = r.end;
      }
    }
    const before = this._value.value;
    const s = Math.min(a.start, start);
    const e = Math.max(a.end, end);
    const loc = rangesLocation(a.system, [s, e]);
    this._anchor = a;
    this._value.set({ elements: EmptyElements, location: loc, shown: loc });
    return before !== this._value.value;
  }

  /** Ends a drag-select (what it selected stays). */
  endDrag(): void {
    this._drag = undefined;
  }

  get dragging(): boolean {
    return this._drag !== undefined;
  }

  // ---- the cursor ----

  setCursor(cursor: Cursor | undefined): void {
    this._cursor.set(cursor);
  }

  // ---- internals ----

  private _ofHit(hit: Hit, range: Anchor | undefined): Emphasized {
    if (hit.element) {
      if (hit.extent) this._extents.set(keyId(hit.element), hit.extent);
      return this._compose(elements([hit.element]), EmptyLocation);
    }
    const loc = rangesLocation(range!.system, [range!.start, range!.end]);
    return { elements: EmptyElements, location: loc, shown: loc };
  }

  /** A value of elements and positions: what is drawn is the positions and each element's extent. */
  private _compose(els: Elements, loc: Location): Emphasized {
    let shown = loc;
    const unknown: ElementKey[] = [];
    for (const k of els.keys) {
      const e = this._extents.get(keyId(k));
      if (e) shown = unionLocations(shown, e);
      else unknown.push(k);
    }
    if (unknown.length && this.extentOf) shown = unionLocations(shown, this.extentOf(elements(unknown)));
    // the extents kept are those still selected
    if (this._extents.size > els.keys.length) {
      const keep = new Set(els.keys.map(keyId));
      for (const id of [...this._extents.keys()]) if (!keep.has(id)) this._extents.delete(id);
    }
    return { elements: els, location: loc, shown: isEmptyLocation(shown) ? EmptyLocation : shown };
  }
}

export class HighlightModel extends ReactiveModel {
  private readonly _channels = new KeyedState<string, Emphasized>();
  private readonly _names = new State<readonly string[]>([]);
  private readonly _version = new State(0);
  /** The channels that hold something, in the order they were first set. */
  readonly names: ReadonlyState<readonly string[]> = this._names.asReadonly;
  /** Moves on every change of any channel: what draws them listens to this. */
  readonly version: ReadonlyState<number> = this._version.asReadonly;

  /** A channel as a state: only its readers are told of its changes. */
  channel(name: string): ReadonlyState<Emphasized | undefined> {
    return this._channels.key(name);
  }

  get(name: string): Emphasized | undefined {
    return this._channels.get(name);
  }

  set(name: string, value: Emphasized): void {
    const old = this._channels.get(name);
    if (old && same(old, value)) return;
    this._channels.set(name, value);
    if (!this._names.value.includes(name)) this._names.set([...this._names.value, name]);
    this._version.set(this._version.value + 1);
  }

  clear(name: string): void {
    if (!this._channels.get(name)) return;
    this._channels.delete(name);
    this._names.set(this._names.value.filter((n) => n !== name));
    this._version.set(this._version.value + 1);
  }
}
