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
  location,
  locationsEqual,
  partOf,
  rangesLocation,
  subtractLocations,
  unionLocations,
} from '@molstar/seqx-core/query';
import type { Axis } from '@molstar/seqx-core/view';
import { axisTotal, cellOffset, positionAt } from './axis.js';
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
  /** The keyboard moved it: drawn (focus-visible); a click places it unseen. */
  readonly keyboard?: boolean;
}

const sameCursor = (a: Cursor | undefined, b: Cursor | undefined) =>
  a === b ||
  (!!a &&
    !!b &&
    a.section === b.section &&
    a.system === b.system &&
    a.position === b.position &&
    !!a.keyboard === !!b.keyboard);

const keyId = (k: ElementKey) => `${typeof k.key === 'number' ? 'n' : 's'}${k.key}\u0000${k.owner}`;

/**
 * The positions a hit covers on its own system, in display order on `axis` (an item over a cyclic origin starts at its
 * first piece shown): an item's extent there (`extent`, else the hit's), else its cell; undefined for neither. `start`
 * is the first position, `end` one past the last, so `end` is below `start` for a span over the origin.
 */
export function hitRange(hit: Hit, axis?: Axis, extent: Location | undefined = hit.extent): Anchor | undefined {
  const part = extent ? partOf(extent, hit.system) : undefined;
  // display order on the axis only when it is the hit's system's
  if (part?.ranges.length)
    return { system: hit.system, ...spanOf(part.ranges, axis?.system === hit.system ? axis : undefined) };
  // a letter through a path (a residue a codon wide): its element on its own system, which the selection projects
  // onto the section (its codon, by closure); extends by elements there
  const other = !hit.element && extent?.parts.length === 1 ? extent.parts[0] : undefined;
  if (other?.ranges.length)
    return { system: other.system, start: other.ranges[0], end: other.ranges[other.ranges.length - 1] };
  if (hit.position === undefined || Number.isNaN(hit.position)) return undefined;
  const p = Math.floor(hit.position);
  return { system: hit.system, start: p, end: p + 1 };
}

/** A part's ranges as one span in display order on `axis`: its first position shown and one past its last. */
function spanOf(ranges: ArrayLike<number>, axis: Axis | undefined): { start: number; end: number } {
  const linear = { start: ranges[0], end: ranges[ranges.length - 1] };
  if (!axis?.ranges) return linear;
  let lo = Infinity;
  let hi = -Infinity;
  let start = linear.start;
  let end = linear.end;
  for (let i = 0; i + 1 < ranges.length; i += 2) {
    const a = cellOffset(axis, ranges[i]);
    const b = cellOffset(axis, ranges[i + 1] - 1);
    if (a !== undefined && a < lo) [lo, start] = [a, ranges[i]];
    if (b !== undefined && b > hi) [hi, end] = [b, ranges[i + 1]];
  }
  return lo === Infinity || hi === -Infinity ? linear : { start, end };
}

/**
 * The positions from one span to another, both included, on their system: in display offsets on `axis` and split into
 * its pieces (a range over a cyclic origin is two ranges; on a reversed axis each piece read down), so a drag or an
 * extend across a break selects the positions shown between its ends (`selection.acrossBreaks: shown`), or with `span`
 * the hidden stretches between them too; without an axis, or a span it does not show, the hull.
 */
function between(
  system: string,
  a: Anchor,
  b: Anchor,
  axis: Axis | undefined,
  across: 'shown' | 'span' = 'shown',
  arc = false,
): Location {
  // an axis of another system (a residue's range on its gene's section) reads as no axis: the hull on the system
  const r = axis?.system === system ? axis.ranges : undefined;
  if (r) {
    const ends = [
      cellOffset(axis!, a.start),
      cellOffset(axis!, a.end - 1),
      cellOffset(axis!, b.start),
      cellOffset(axis!, b.end - 1),
    ];
    if (ends.every((o) => o !== undefined)) {
      const lo = Math.min(...(ends as number[]));
      const hi = Math.max(...(ends as number[])) + 1;
      const total = axisTotal(axis!) ?? 0;
      // round a ring (`arc`): the shorter way between the ends, over the origin when that is shorter
      if (arc && axis!.cyclic && axis!.dir === 1 && total > 0 && hi - lo > total / 2) {
        const p0 = positionAt(axis!, hi - 1 + 0.5);
        const p1 = positionAt(axis!, lo + 0.5);
        if (p0 !== undefined && p1 !== undefined) {
          const from = Math.floor(p0);
          const to = Math.floor(p1) + 1;
          const len = axis!.length ?? total;
          return rangesLocation(system, from < to ? [from, to] : [from, len, 0, to]);
        }
      }
      const out: number[] = [];
      let acc = 0;
      for (let i = 0; i + 1 < r.length; i += 2) {
        const w = r[i + 1] - r[i];
        const s = Math.max(lo, acc);
        const e = Math.min(hi, acc + w);
        if (e > s) {
          if (axis!.dir === -1) out.push(r[i + 1] - (e - acc), r[i + 1] - (s - acc));
          else out.push(r[i] + s - acc, r[i] + e - acc);
        }
        acc += w;
      }
      if (out.length) {
        const loc = rangesLocation(system, out);
        if (across === 'shown' || axis!.cyclic) return loc;
        // the hidden stretches between the ends too: the hull of what is shown between them
        const part = loc.parts[0].ranges;
        return rangesLocation(system, [part[0], part[part.length - 1]]);
      }
    }
  }
  const lin = (x: Anchor) =>
    x.end > x.start ? [x.start, x.end] : [Math.min(x.start, x.end - 1), Math.max(x.start + 1, x.end)];
  const [as, ae] = lin(a);
  const [bs, be] = lin(b);
  return rangesLocation(system, [Math.min(as, bs), Math.max(ae, be)]);
}

/** A span as a Location: two ranges when it runs over a cyclic origin (`end` below `start`, `length` known). */
function spanLocation(a: Anchor, axis: Axis | undefined): Location {
  return between(a.system, a, a, axis);
}

/** A cell of two axes (a grid's): a position on each, x's system and y's. */
export interface GridCell {
  readonly system: string;
  readonly position: number;
  readonly ySystem: string;
  readonly yPosition: number;
}

/** The rect from one cell to another, both included: x's positions between theirs × y's. */
export function cellsRect(a: GridCell, b: GridCell): Location {
  const x0 = Math.floor(Math.min(a.position, b.position));
  const x1 = Math.floor(Math.max(a.position, b.position)) + 1;
  const y0 = Math.floor(Math.min(a.yPosition, b.yPosition));
  const y1 = Math.floor(Math.max(a.yPosition, b.yPosition)) + 1;
  return location({
    rects: [{ a: { system: a.system, ranges: [x0, x1] }, b: { system: a.ySystem, ranges: [y0, y1] } }],
  });
}

/** Whether every cell of `r` is in `loc`'s rects. */
function rectInside(r: Location, loc: Location): boolean {
  return subtractLocations(location({ rects: r.rects }), location({ rects: loc.rects })).rects.length === 0;
}

/**
 * How a drag on data changes the selection (the interaction grammar): it becomes the range (`replace`), the range is
 * added to it (`add`) or taken out of it (`subtract`), or the range the anchor is in reaches to the pointer (`extend`).
 */
export type DragMode = 'replace' | 'add' | 'subtract' | 'extend';

/**
 * What an extend (Shift+click, Shift+drag, Shift+arrows) changes: the range the anchor is in, the others kept (`range`),
 * or the whole selection, which becomes that range (`selection`).
 */
export type ExtendScope = 'range' | 'selection';

export interface SelectionOptions {
  /** `range` by default (`selection.extendScope`). */
  readonly extendScope?: ExtendScope;
  /** What a drag or an extend across a break selects: the positions shown (`shown`, the default) or the span. */
  readonly acrossBreaks?: 'shown' | 'span';
}

/** What an extend keeps: the selection less the range the anchor is in, set when the anchor is. */
interface Rest {
  readonly elements: Elements;
  readonly location: Location;
}

const NoRest: Rest = Object.freeze({ elements: EmptyElements, location: EmptyLocation });

export class SelectionModel extends ReactiveModel {
  private readonly _value = new State<Emphasized>(Nothing, same);
  private readonly _cursor = new State<Cursor | undefined>(undefined, sameCursor);
  /** The extents of elements selected by a hit, by key: what is drawn of them without resolving them again. */
  private readonly _extents = new Map<string, Location>();
  private _anchor: Anchor | undefined;
  /** What an extend from the anchor keeps (the other ranges and elements); none after a click sets the selection. */
  private _rest: Rest = NoRest;
  private _drag:
    | {
        anchor: Anchor;
        snap: boolean;
        mode: DragMode;
        before: Emphasized;
        beforeAnchor: Anchor | undefined;
        beforeRest: Rest;
      }
    | undefined;
  /** Where a rect's extend starts (a grid's cell), and a rect drag's first cell and the selection when it began. */
  private _cellAnchor: GridCell | undefined;
  private _rectDrag: { from: GridCell; mode: DragMode; before: Emphasized } | undefined;
  private readonly _scope: ExtendScope;
  private readonly _across: 'shown' | 'span';
  /** The selection: an equal value tells nothing. */
  readonly value: ReadonlyState<Emphasized> = this._value.asReadonly;
  /** The keyboard cursor; an equal one tells nothing. */
  readonly cursor: ReadonlyState<Cursor | undefined> = this._cursor.asReadonly;

  /**
   * `extentOf`: an element's extent when no hit gave it (one selected by an address, then toggled): the viewer's
   * addresses resolve it. Without one, an element with no known extent is held and not drawn.
   */
  constructor(
    private readonly extentOf?: (elements: Elements) => Location,
    options: SelectionOptions = {},
  ) {
    super();
    this._scope = options.extendScope ?? 'range';
    this._across = options.acrossBreaks ?? 'shown';
  }

  /** Where an extend starts: the last element or cell set, toggled or begun from (none after a clear). */
  get anchor(): Anchor | undefined {
    return this._anchor;
  }

  /** Sets the selection; with an anchor, an extend from it keeps the rest of the value (less the anchor's span). */
  set(value: Emphasized, anchor?: Anchor): void {
    this._anchor = anchor;
    this._rest = anchor
      ? { elements: value.elements, location: subtractLocations(value.location, spanLocation(anchor, undefined)) }
      : NoRest;
    this._value.set(value);
  }

  clear(): void {
    this._anchor = undefined;
    this._cellAnchor = undefined;
    this._rectDrag = undefined;
    this._rest = NoRest;
    this._drag = undefined;
    this._extents.clear();
    this._value.set(Nothing);
  }

  // ---- the grammar's gestures, by the hit (the interaction grammar) ----

  /**
   * A click: the selection becomes what was hit (an item: its key; a cell: its position); clicking the only selected
   * element clears it. Nothing hit clears it. True when the selection changed.
   */
  select(hit: Hit | undefined, axis?: Axis): boolean {
    const before = this._value.value;
    const range = hit ? hitRange(hit, axis, this._extentOfHit(hit)) : undefined;
    if (!hit || (!hit.element && !range)) {
      this.clear();
      return before !== this._value.value;
    }
    const next = this._ofHit(hit, range, axis);
    // what is selected, wherever it was set from (a link's DBD clears on the first click on it)
    if (elementsEqual(next.elements, before.elements) && locationsEqual(next.location, before.location)) this.clear();
    else {
      this._anchor = range;
      this._rest = NoRest;
      this._value.set(next);
    }
    return before !== this._value.value;
  }

  /**
   * Shift+click: the range the anchor is in reaches from the anchor to what was hit, on the anchor's system, the other
   * ranges and elements kept (`selection.extendScope`: `selection` makes it the whole selection); with no anchor
   * there, a click. The anchor stays, so extends from it replace one another.
   */
  extend(hit: Hit | undefined, axis?: Axis): boolean {
    const a = this._anchor;
    const range = hit ? hitRange(hit, axis, this._extentOfHit(hit)) : undefined;
    if (!a || !range || range.system !== a.system) return hit ? this.select(hit, axis) : false;
    const before = this._value.value;
    this._value.set(this._extended(between(a.system, a, range, axis, this._across)));
    return before !== this._value.value;
  }

  /** Cmd/Ctrl+click: what was hit in or out of the selection, the rest kept; it becomes the anchor. */
  toggle(hit: Hit | undefined, axis?: Axis): boolean {
    const range = hit ? hitRange(hit, axis, this._extentOfHit(hit)) : undefined;
    if (!hit || (!hit.element && !range)) return false;
    const before = this._value.value;
    let els = before.elements;
    let loc = before.location;
    let cell = EmptyLocation;
    if (hit.element) {
      const has = els.keys.some((k) => k.owner === hit.element!.owner && k.key === hit.element!.key);
      els = elements(
        has
          ? els.keys.filter((k) => k.owner !== hit.element!.owner || k.key !== hit.element!.key)
          : [...els.keys, hit.element],
      );
    } else if (range) {
      cell = spanLocation(range, axis);
      const inside = locationsEqual(subtractLocations(cell, loc), EmptyLocation);
      loc = inside ? subtractLocations(loc, cell) : unionLocations(loc, cell);
    }
    this._anchor = range;
    // an extend from here keeps the rest: all but what was hit
    const el = hit.element;
    this._rest = {
      elements: el ? elements(els.keys.filter((k) => k.owner !== el.owner || k.key !== el.key)) : els,
      location: subtractLocations(loc, cell),
    };
    this._value.set(this._compose(els, loc));
    return before !== this._value.value;
  }

  /**
   * Adds positions to the selection, the rest kept (the `selection.add` command). With `anchor` (the added range's), a
   * later extend starts there, keeping what was there before; without one (an address's positions) the anchor stays
   * and what an extend keeps gains what was added. True when the selection changed.
   */
  add(loc: Location, anchor?: Anchor): boolean {
    const before = this._value.value;
    if (anchor) {
      this._anchor = anchor;
      this._rest = { elements: before.elements, location: before.location };
    } else if (this._anchor) {
      this._rest = { elements: this._rest.elements, location: unionLocations(this._rest.location, loc) };
    }
    this._value.set(this._compose(before.elements, unionLocations(before.location, loc)));
    return before !== this._value.value;
  }

  /**
   * Takes positions out of the selection (the `selection.subtract` command): an element they cut becomes the positions
   * of its extent left; one they miss stays an element. The anchor stays; what an extend keeps is what is left but the
   * run the anchor is in, so a part split off stays.
   */
  subtract(loc: Location): boolean {
    const before = this._value.value;
    const next = this._without(before, loc);
    this._rest = this._restAround(next);
    this._value.set(next);
    return before !== this._value.value;
  }

  /**
   * A drag-select begins on what was hit: the range grows from it as the pointer moves (`dragTo`), snapped to whole
   * items when it began on one (the interaction grammar); `mode` says what the range does to the selection, and an
   * extend grows from the anchor (from what was hit when there is no anchor on its system: a replace). False when
   * nothing selectable was hit.
   */
  beginDrag(hit: Hit | undefined, axis?: Axis, mode: DragMode = 'replace'): boolean {
    const range = hit ? hitRange(hit, axis, this._extentOfHit(hit)) : undefined;
    if (!range) return false;
    const a = this._anchor;
    const extend = mode === 'extend' && !!a && a.system === range.system;
    this._drag = {
      anchor: extend ? a! : range,
      snap: !!hit?.element,
      mode: mode === 'extend' && !extend ? 'replace' : mode,
      before: this._value.value,
      beforeAnchor: this._anchor,
      beforeRest: this._rest,
    };
    return true;
  }

  /**
   * The drag's other end: a position on the drag's system (the pointer's, read from the view even past the data), or
   * an item's extent when the drag began on items and the pointer is over one. The range between the two ends replaces
   * the selection, is added to or taken out of what it was when the drag began, or is the anchor's range (its mode).
   * The range's start is the anchor for a later extend (a subtraction keeps the anchor). True when the selection changed.
   */
  dragTo(position: number, hit?: Hit, axis?: Axis, arc = false): boolean {
    const d = this._drag;
    if (!d || Number.isNaN(position)) return false;
    const a = d.anchor;
    const p = Math.floor(position);
    let cur: Anchor = { system: a.system, start: p, end: p + 1 };
    if (d.snap && hit?.element && hit.system === a.system) cur = hitRange(hit, axis, this._extentOfHit(hit)) ?? cur;
    if (axis && axis.system !== a.system) {
      // a drag begun on a residue of another system (a residue row): whole residues, from the one under the pointer
      const r = hit ? hitRange(hit, axis, this._extentOfHit(hit)) : undefined;
      if (!r || r.system !== a.system) return false;
      cur = r;
    }
    const before = this._value.value;
    const loc = between(a.system, a, cur, axis, this._across, arc);
    switch (d.mode) {
      case 'replace':
        this._anchor = a;
        this._rest = NoRest;
        this._value.set({ elements: EmptyElements, location: loc, shown: loc });
        break;
      case 'add':
        this._anchor = a;
        this._rest = { elements: d.before.elements, location: d.before.location };
        this._value.set(this._compose(d.before.elements, unionLocations(d.before.location, loc)));
        break;
      case 'subtract': {
        const next = this._without(d.before, loc);
        this._rest = this._restAround(next);
        this._value.set(next);
        break;
      }
      case 'extend':
        this._value.set(this._extended(loc));
        break;
    }
    return before !== this._value.value;
  }

  /** Cancels a drag-select (Esc): the selection, the anchor and what an extend keeps, as the drag began. */
  cancelDrag(): void {
    const d = this._drag;
    if (!d) return;
    this._drag = undefined;
    this._anchor = d.beforeAnchor;
    this._rest = d.beforeRest;
    this._value.set(d.before);
  }

  /** Ends a drag-select (what it selected stays). */
  endDrag(): void {
    this._drag = undefined;
  }

  get dragging(): boolean {
    return this._drag !== undefined;
  }

  // ---- rects: a grid's cells (_spec/README.md "SelectionModel and HighlightModel", Rects) ----

  /**
   * A click on a cell of two axes: the selection becomes the cell's rect (`replace`), the rect from the cell anchor to
   * it (`extend`), or the cell is toggled in the selection's rects (`toggle`: taken out when all of it is in, else
   * added). The cell becomes the anchor but for an extend. True when the selection changed.
   */
  selectCell(cell: GridCell, how: 'replace' | 'extend' | 'toggle' = 'replace'): boolean {
    const before = this._value.value;
    const a = this._cellAnchor;
    if (how === 'extend' && a && a.system === cell.system && a.ySystem === cell.ySystem) {
      const r = cellsRect(a, cell);
      this._value.set({ elements: EmptyElements, location: r, shown: r });
      return before !== this._value.value;
    }
    const r = cellsRect(cell, cell);
    this._cellAnchor = cell;
    this._anchor = undefined;
    this._rest = NoRest;
    if (how === 'toggle') {
      const loc = before.location;
      this._value.set(
        rectInside(r, loc) ? this._without(before, r) : this._compose(before.elements, unionLocations(loc, r)),
      );
    } else if (locationsEqual(before.location, r) && !before.elements.keys.length) this.clear();
    else this._value.set({ elements: EmptyElements, location: r, shown: r });
    return before !== this._value.value;
  }

  /** A rect drag begins on a cell: it replaces, adds to, takes out of the selection, or extends from the cell anchor. */
  beginRectDrag(cell: GridCell, mode: DragMode = 'replace'): void {
    const a = this._cellAnchor;
    const extend = mode === 'extend' && !!a && a.system === cell.system && a.ySystem === cell.ySystem;
    this._rectDrag = {
      from: extend ? a! : cell,
      mode: mode === 'extend' && !extend ? 'replace' : mode,
      before: this._value.value,
    };
    if (!extend) this._cellAnchor = cell;
  }

  /**
   * The rect drag's other cell: the rect between the two replaces the selection, is added to or taken out of it as it
   * was when the drag began (rects subtract: at most two rects, the cells exact). True when the selection changed.
   */
  rectDragTo(position: number, yPosition: number): boolean {
    const d = this._rectDrag;
    if (!d || Number.isNaN(position) || Number.isNaN(yPosition)) return false;
    const before = this._value.value;
    const r = cellsRect(d.from, { ...d.from, position, yPosition });
    switch (d.mode) {
      case 'replace':
      case 'extend':
        this._value.set({ elements: EmptyElements, location: r, shown: r });
        break;
      case 'add':
        this._value.set(this._compose(d.before.elements, unionLocations(d.before.location, r)));
        break;
      case 'subtract':
        this._value.set(this._without(d.before, r));
        break;
    }
    this._anchor = undefined;
    this._rest = NoRest;
    return before !== this._value.value;
  }

  /** Ends a rect drag (what it selected stays); with `cancel`, the selection as it began. */
  endRectDrag(cancel = false): void {
    const d = this._rectDrag;
    this._rectDrag = undefined;
    if (cancel && d) this._value.set(d.before);
  }

  get rectDragging(): boolean {
    return this._rectDrag !== undefined;
  }

  // ---- the cursor ----

  setCursor(cursor: Cursor | undefined): void {
    this._cursor.set(cursor);
  }

  // ---- internals ----

  /**
   * What an extend keeps after a subtraction: the value less the run of positions the anchor is in on its system (all
   * of it when the anchor is in none; nothing without an anchor).
   */
  private _restAround(value: Emphasized): Rest {
    const a = this._anchor;
    if (!a) return NoRest;
    const r = partOf(value.location, a.system)?.ranges;
    if (r)
      for (let i = 0; i + 1 < r.length; i += 2)
        if (a.start >= r[i] && a.start < r[i + 1])
          return {
            elements: value.elements,
            location: subtractLocations(value.location, rangesLocation(a.system, [r[i], r[i + 1]])),
          };
    return { elements: value.elements, location: value.location };
  }

  /** An extend's value: the anchor's range `loc` with what the extend keeps (`selection.extendScope`). */
  private _extended(loc: Location): Emphasized {
    const r = this._rest;
    if (this._scope === 'selection' || (!r.elements.keys.length && isEmptyLocation(r.location)))
      return { elements: EmptyElements, location: loc, shown: loc };
    return this._compose(r.elements, unionLocations(r.location, loc));
  }

  /** A value less some positions: an element they cut becomes what is left of its extent; one they miss stays. */
  private _without(value: Emphasized, loc: Location): Emphasized {
    let positions = subtractLocations(value.location, loc);
    const keep: ElementKey[] = [];
    for (const k of value.elements.keys) {
      const e = this._extents.get(keyId(k)) ?? this.extentOf?.(elements([k]));
      const left = e ? subtractLocations(e, loc) : undefined;
      if (!e || !left || locationsEqual(left, e)) keep.push(k);
      else positions = unionLocations(positions, left);
    }
    return this._compose(keep.length === value.elements.keys.length ? value.elements : elements(keep), positions);
  }

  private _ofHit(hit: Hit, range: Anchor | undefined, axis: Axis | undefined): Emphasized {
    if (hit.element) return this._compose(elements([hit.element]), EmptyLocation);
    const loc = spanLocation(range!, axis);
    return { elements: EmptyElements, location: loc, shown: loc };
  }

  /** An item's whole extent for a hit (`_extentOfHit`), kept while it is selected; undefined for a cell. */
  extentOfHit(hit: Hit): Location | undefined {
    return this._extentOfHit(hit);
  }

  /**
   * An item's whole extent, kept while it is selected: the hit's (a layer gives the piece hit, one of a row wrapped over
   * a cyclic origin) with what `extentOf` resolves of its key. Undefined for a cell.
   */
  private _extentOfHit(hit: Hit): Location | undefined {
    const el = hit.element;
    if (!el) return undefined;
    const id = keyId(el);
    const kept = this._extents.get(id);
    if (
      kept &&
      (!hit.extent || kept === hit.extent || locationsEqual(subtractLocations(hit.extent, kept), EmptyLocation))
    )
      return kept;
    let e = hit.extent ?? EmptyLocation;
    if (this.extentOf) e = unionLocations(e, this.extentOf(elements([el])));
    if (isEmptyLocation(e)) return hit.extent;
    this._extents.set(id, e);
    return e;
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
