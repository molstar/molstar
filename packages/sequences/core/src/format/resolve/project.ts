/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The projection core (seqspec:interaction#channels; seqspec:model#projection, seqspec:model#closure; CP1b),
// platform-free: both viewers call it. Extents in coordinate systems are closed by containment over the whole view's
// span maps (`ResolvedView.systems.closure`), widened by a selection's granularity, and projected onto a block as a
// plain union of block-system runs through coordinate maps, along the first route kind that applies: direct, the view's
// links, or the model's shortest paths. Routes are composed once per (source system, target) on the instance; a row
// source goes through its own row's binding to its bound unit and reuses that unit's routes. After a route's first use
// nothing on the per-move paths allocates: the walks read the maps (`Coords`) with module visitors and contexts on the
// instance.

import {
  type CoordinateMap,
  type CoordinateSystem,
  Coords,
  Extents,
  CoordinateSystemKind as K,
  type MapHit,
  Mappings,
  type SeqSet,
} from '@molstar/seq-core/model';
import { reach } from './reach.js';
import type { BlockRow } from './section.js';
import type { ResolvedView } from './view.js';

/** A channel's closure option (seqspec:interaction#hover): `containing` adds the elements that contain each element. */
export type Closure = 'element' | 'containing';
/** How a selection widens after closure, each entry in its own coordinate system (seqspec:interaction#selection). */
export const SelectionGranularity = { position: 0, range: 1, unit: 2, copyGroup: 3, column: 4, row: 5 } as const;

/** Where projected runs go, in any order (a `MarkSet` is one). */
export interface RunSink {
  add(start: number, stop: number): void;
  /**
   * A boundary lit (seqspec:model#projection): the place before index `at`, an insertion's, between two positions. A
   * sink without it (a focus, a host's) lights the position after it.
   */
  boundary?(at: number): void;
}

/** The projection over one resolved view and its set; routes cached on it. Not re-entrant. */
export interface Projection {
  readonly view: ResolvedView;
  readonly set: SeqSet;
  /** Step 2: appends, for each entry, the elements of the view's other systems that contain its elements. One step. */
  close(ext: Extents): void;
  /** Step 3, the selection's: widens each entry in its own system (`SelectionGranularity`; `range` is the viewer's). */
  widen(ext: Extents, granularity: number): void;
  /** Step 4: every entry onto a block, as block-system runs; a boundary, or an insert's place, as a sink's caret. */
  project(ext: Extents, block: number, out: RunSink): void;
  /** Routes composed so far (a counter for tests: row sources add none). */
  readonly routes: number;
}

export function projectionOf(view: ResolvedView, set: SeqSet): Projection {
  return new ViewProjection(view, set);
}

/** A route's leg: the input clipped to `[lo, hi)`, then through `map` (`read` 1: from → to; -1: the inverse image). */
interface Leg {
  readonly map: CoordinateMap | null;
  readonly read: number;
  readonly lo: number;
  readonly hi: number;
}
/** The legs of the first route kind that applies (none: lights nothing), and the target's domain. */
interface Route {
  readonly legs: readonly Leg[];
  readonly lo: number;
  readonly hi: number;
}
/** A block, or a closure map's `to` system no block shows (a merged member). */
interface Target {
  readonly unit: number;
  readonly alignment: number;
  readonly members: Int32Array | null;
  readonly block: number;
  readonly lo: number;
  readonly hi: number;
}

/**
 * A walk's output: clipped to the target's domain, `spanned` once a leg went through a span map from its residues; `n`
 * counts what it passed on.
 */
class Emit {
  sink: RunSink = NoSink;
  lo = 0;
  hi = 0;
  spanned = false;
  n = 0;
}
const NoSink: RunSink = { add() {} };

/** Reused runs, merged when they touch in order: a closure's scratch. */
class Runs implements RunSink {
  count = 0;
  start = new Int32Array(8);
  stop = new Int32Array(8);
  add(s: number, t: number) {
    const n = this.count;
    if (n > 0 && s <= this.stop[n - 1] && t >= this.start[n - 1]) {
      if (s < this.start[n - 1]) this.start[n - 1] = s;
      if (t > this.stop[n - 1]) this.stop[n - 1] = t;
      return;
    }
    if (n === this.start.length) {
      const [a, b] = [new Int32Array(2 * n), new Int32Array(2 * n)];
      a.set(this.start);
      b.set(this.stop);
      [this.start, this.stop] = [a, b];
    }
    this.start[n] = s;
    this.stop[n] = t;
    this.count = n + 1;
  }
}

/**
 * A boundary (s === t) is the sink's caret, within the target's domain; for a sink without carets the position after
 * it, the last at the end. Runs are clipped to the target's domain.
 */
function emit(s: number, t: number, c: Emit) {
  if (s === t) {
    const sink = c.sink;
    if (sink.boundary !== undefined) {
      if (s >= c.lo && s <= c.hi) {
        c.n++;
        sink.boundary(s);
      }
      return;
    }
    if (s >= c.hi) s = c.hi - 1;
    t = s + 1;
  }
  if (s < c.lo) s = c.lo;
  if (t > c.hi) t = c.hi;
  if (s < t) {
    c.n++;
    c.sink.add(s, t);
  }
}
/** A bound row's run, landed on its unit: clipped to the unit's range, then that unit's routes. */
function viaUnit(s: number, t: number, p: ViewProjection) {
  const r = p.row;
  const a = s < r.lo ? r.lo : s;
  const b = t > r.hi ? r.hi : t;
  if (a < b) p.walk(K.unit, r.unit, a, b, r.target, r.emit);
}
/** A closure map's inverse image: the elements of its `from` system that contain the runs, appended. */
function addContaining(s: number, t: number, p: ViewProjection) {
  p.append(s, t);
}
/** A candidate containing element: kept when the element's whole extent lies within the candidate's. */
function addIfContains(h: MapHit, p: ViewProjection) {
  const c = p.cand;
  const m = c.map as CoordinateMap;
  const runs = p.scratch;
  for (let r = 0; r < runs.count; r++) {
    for (let x = runs.start[r]; x < runs.stop[r]; x++) {
      let held = false;
      for (let q = 0; q < m.span && !held; q++) held = Coords.toAt(m, h.element, q) === x;
      if (!held) return;
    }
  }
  if (c.self && h.element === c.element) return;
  p.append(h.element, h.element + 1);
}

/**
 * Whether a target shows a source system: its own, a merged member (`members`; not for links, whose paths start from
 * the block's own system, so a copy never reaches another through a merged block), or its unit's molecule.
 */
const shows = (set: SeqSet, t: Target, kind: number, id: number, members: boolean) =>
  kind === K.columns
    ? t.alignment === id
    : t.alignment < 0 &&
      (kind === K.unit
        ? t.unit === id || (members && t.members?.includes(id) === true)
        : kind === K.molecule && set.units.molecule[t.unit] === id);
/** A path read from its `from` to its `to` (1), or back (-1), as a leg's `read` on its composed map. */
const readOf = (map: CoordinateMap, from: CoordinateSystem, forward: boolean) =>
  Coords.sameSystem(map.from, from) === forward ? 1 : -1;

/** `widen`'s scratch systems (CF-15): `Extents.add` copies what it is given. */
const WidenSystem = { kind: 0, id: 0, row: -1 };
const WidenColumns = { kind: K.columns as number, id: 0, row: -1 };

class ViewProjection implements Projection {
  readonly view: ResolvedView;
  readonly set: SeqSet;
  routes = 0;
  private readonly targets: Target[] = [];
  /** Per source kind (molecule, unit, columns), per id: a route per target, composed on first use. */
  private readonly cache = [new Map<number, Route[]>(), new Map<number, Route[]>(), null, new Map<number, Route[]>()];
  /** The closure maps by `to` system: its target, and its rows `[start, end)` of `systems.closure`. */
  private readonly closureTo: { target: number; start: number; end: number }[] = [];
  private readonly maps: CoordinateMap[];
  private readonly emitOut = new Emit();
  private readonly emitScratch = new Emit();
  readonly scratch = new Runs();
  readonly row = { unit: 0, target: 0, lo: 0, hi: 0, emit: this.emitOut };
  /** The closure in progress: the extents appended to, from which entry on, the system added, a candidate's test. */
  readonly cand = {
    ext: null as Extents | null,
    from: 0,
    system: { kind: 0, id: 0, row: -1 } as CoordinateSystem,
    map: null as CoordinateMap | null,
    self: false,
    element: 0,
  };
  private readonly hit: MapHit = { element: 0, phase: 0 };

  constructor(view: ResolvedView, set: SeqSet) {
    this.view = view;
    this.set = set;
    const b = view.blocks;
    for (let k = 0; k < b.count; k++) this.target(b.unit[k], b.alignment[k], b.members[k], k);
    const c = view.systems.closure;
    const s = view.systems;
    this.maps = Array.from(c.path, (p) => view.paths[p].map);
    for (let k = 0; k < c.count; k++) {
      const to = c.to[k];
      if (k > 0 && to === c.to[k - 1]) {
        (this.closureTo.at(-1) as { end: number }).end = k + 1;
        continue;
      }
      const [kind, id] = [s.kind[to], s.id[to]];
      let t = -1;
      for (let x = 0; x < b.count && t < 0; x++) {
        const unit = b.alignment[x] < 0 && b.unit[x] === id;
        if (kind === K.unit ? unit : kind === K.columns && b.alignment[x] === id) t = x;
      }
      if (t < 0 && kind === K.unit) t = this.target(id, -1, null, -1);
      if (t >= 0) this.closureTo.push({ target: t, start: k, end: k + 1 });
    }
  }

  private target(unit: number, alignment: number, members: Int32Array | null, block: number): number {
    const system = { kind: alignment >= 0 ? K.columns : K.unit, id: alignment >= 0 ? alignment : unit, row: -1 };
    const [lo, hi] = Mappings.domainOf(this.set, system);
    return this.targets.push({ unit, alignment, members, block, lo, hi }) - 1;
  }

  project(ext: Extents, block: number, out: RunSink): void {
    const c = this.emitOut;
    c.sink = out;
    for (let i = 0; i < ext.count; i++)
      this.entry(ext.kind[i], ext.id[i], ext.row[i], ext.start[i], ext.stop[i], block, c);
    c.sink = NoSink;
  }

  /** One entry onto a target: a row through its own row only, a molecule, unit or columns along its route. */
  private entry(kind: number, id: number, row: number, s: number, t: number, target: number, c: Emit) {
    if (kind === K.row) this.walkRow(id, row, s, t, target, c);
    else if (kind === K.molecule || kind === K.unit || kind === K.columns) this.walk(kind, id, s, t, target, c);
  }

  walk(kind: number, id: number, s: number, t: number, target: number, c: Emit) {
    const r = this.route(kind, id, target);
    for (let i = 0; i < r.legs.length; i++) {
      const leg = r.legs[i];
      const a = s < leg.lo ? leg.lo : s;
      const b = t > leg.hi ? leg.hi : t;
      if (s === t ? s < leg.lo || s > leg.hi : a >= b) continue;
      c.lo = r.lo;
      c.hi = r.hi;
      const m = leg.map;
      const before = c.n;
      const caret = c.sink.boundary !== undefined;
      if (m === null) emit(a, b, c);
      else if (leg.read > 0) {
        if (m.span > 1) c.spanned = true;
        Coords.forEachTo(m, a, b, emit, c);
      } else if (a === b && caret) Coords.boundaryFrom(m, a, emit, c);
      else Coords.forEachFrom(m, a, a === b ? a + 1 : b, emit, c);
      // An insert (unpaired between pairs: a design's grafted loop, a stretch a design deletes, an intron), or a
      // boundary inside one, lights the caret before the next pair on the other side (seqspec:model#projection).
      if (m === null || !caret || c.n !== before) continue;
      const e = Coords.pairedAfter(m, a, b, leg.read > 0);
      if (e >= 0 && leg.read > 0) Coords.forEachTo(m, e, e, emit, c);
      else if (e >= 0) Coords.boundaryFrom(m, e, emit, c);
    }
  }

  /**
   * A row source: on its own alignment's columns through the row's map; elsewhere through its binding to its bound
   * unit, then that unit's routes. An unbound row, or one bound more than once, lights nothing elsewhere.
   */
  private walkRow(alignment: number, row: number, s: number, t: number, target: number, c: Emit) {
    const T = this.targets[target];
    if (T.alignment === alignment) {
      const rows = this.set.alignments[alignment]?.rows;
      // A few residues (a hover) read their columns one by one, allocation-free; more, the row's cached map.
      const cols = t - s > 32 ? (rows?.posMap(row as never) ?? null) : null;
      if (!rows || (t - s > 32 && !cols)) return;
      c.lo = T.lo;
      c.hi = T.hi;
      let a = -1;
      let b = -1;
      for (let j = s < 0 ? 0 : s; j < t && (!cols || j < cols.length); j++) {
        const x = cols ? cols[j] : rows.colOf(row as never, j);
        if (x < 0) continue;
        if (x === b) b++;
        else {
          if (a >= 0) emit(a, b, c);
          a = x;
          b = x + 1;
        }
      }
      if (a >= 0) emit(a, b, c);
      return;
    }
    const bind = Mappings.bindingsOf(this.set, alignment);
    const k = bind ? Mappings.boundRow(bind, row) : -1;
    if (!bind || k < 0 || bind.bindings[k] !== 1) return;
    const r = this.row;
    const units = this.set.units;
    for (let e = bind.start[k]; e < bind.start[k + 1]; e++) {
      r.unit = bind.unit[e];
      r.lo = units.rangeStart[r.unit];
      r.hi = units.rangeStop[r.unit];
      r.target = target;
      r.emit = c;
      if (bind.rowIsFrom[e]) Coords.forEachTo(bind.map[e], s, t, viaUnit, this);
      else Coords.forEachFrom(bind.map[e], s, t, viaUnit, this);
    }
  }

  private route(kind: number, id: number, target: number): Route {
    const byId = this.cache[kind] as Map<number, Route[]>;
    let list = byId.get(id);
    if (!list) {
      list = new Array<Route>(this.targets.length);
      byId.set(id, list);
    }
    let r = list[target];
    if (r === undefined) {
      r = this.compose(kind, id, target);
      list[target] = r;
    }
    return r;
  }

  /** Direct; else every link whose other block shows the source; else the shortest paths to the block or its members. */
  private compose(kind: number, id: number, target: number): Route {
    this.routes++;
    const set = this.set;
    const T = this.targets[target];
    const from = { kind, id, row: -1 };
    const [lo, hi] = Mappings.domainOf(set, from);
    const legs: Leg[] = [];
    if (shows(set, T, kind, id, true)) return { legs: [{ map: null, read: 1, lo, hi }], lo: T.lo, hi: T.hi };
    if (this.ownShown(T, kind, id)) return { legs, lo: T.lo, hi: T.hi };
    const l = this.view.links;
    for (let k = 0; T.block >= 0 && k < l.count; k++) {
      const into = l.toBlock[k] === T.block;
      const other = into ? l.fromBlock[k] : l.fromBlock[k] === T.block ? l.toBlock[k] : -1;
      const O = this.targets[other];
      if (other < 0 || !shows(set, O, kind, id, false)) continue;
      const p = this.view.paths[l.path[k]];
      legs.push({ map: p.map, read: readOf(p.map, p.from, into), lo: Math.max(lo, O.lo), hi: Math.min(hi, O.hi) });
    }
    if (!legs.length) {
      const block = { unit: T.unit, alignment: T.alignment, members: T.members } as unknown as BlockRow;
      reach(set, from, block, (_m, p) =>
        legs.push(p ? { map: p.map, read: readOf(p.map, p.from, true), lo, hi } : { map: null, read: 1, lo, hi }),
      );
    }
    return { legs, lo: T.lo, hi: T.hi };
  }

  /**
   * Direct wins (seqspec:format#placement, as placement): a unit or molecule whose molecule a block of the target's
   * section shows lights only that molecule's blocks there, never another molecule's through an alignment or a path.
   */
  private ownShown(T: Target, kind: number, id: number): boolean {
    const mol = kind === K.unit ? this.set.units.molecule[id] : kind === K.molecule ? id : -1;
    if (mol < 0 || T.block < 0) return false;
    const { blocks, sections } = this.view;
    const section = blocks.section[T.block];
    const end = sections.blockStart[section] + sections.blockCount[section];
    for (let b = sections.blockStart[section]; b < end; b++)
      if (blocks.alignment[b] < 0 && this.set.units.molecule[blocks.unit[b]] === mol) return true;
    return false;
  }

  close(ext: Extents): void {
    const n = ext.count;
    const cand = this.cand;
    cand.ext = ext;
    cand.from = n;
    for (let i = 0; i < n; i++) {
      if (ext.start[i] >= ext.stop[i]) continue; // a boundary is no element
      for (let f = 0; f < this.closureTo.length; f++) this.closeIn(ext, i, f);
    }
    cand.ext = null;
    cand.map = null;
  }

  /**
   * Entry `i`'s containing elements through the closure maps into one `to` system: its extent there by the routes; a
   * span-1 extent is held by the inverse image, a residue's extent (a codon) only by an element whose extent holds it.
   */
  private closeIn(ext: Extents, i: number, f: number) {
    const kind = ext.kind[i];
    const id = ext.id[i];
    const row = ext.row[i];
    const s = ext.start[i];
    const t = ext.stop[i];
    const { target, start: c0, end: c1 } = this.closureTo[f];
    const c = this.emitScratch;
    const runs = this.scratch;
    const C = this.view.systems.closure;
    runs.count = 0;
    c.sink = runs;
    c.spanned = false;
    this.entry(kind, id, row, s, t, target, c);
    for (let r = 0; !c.spanned && r < runs.count; r++) {
      const a = runs.start[r];
      const b = runs.stop[r];
      for (let k = this.lastBefore(c0, c1, b); k >= c0 && C.toMaxStop[k] > a; k--) {
        if (C.toStop[k] <= a) continue;
        this.candidate(k, -1, false);
        Coords.forEachFrom(this.maps[k], a, b, addContaining, this);
      }
    }
    for (let x = s; c.spanned && x < t; x++) {
      runs.count = 0;
      this.entry(kind, id, row, x, x + 1, target, c);
      let a = 0x7fffffff;
      let b = -1;
      for (let r = 0; r < runs.count; r++) {
        if (runs.start[r] < a) a = runs.start[r];
        if (runs.stop[r] > b) b = runs.stop[r];
      }
      for (let k = this.lastBefore(c0, c1, b); b > a && k >= c0 && C.toMaxStop[k] > a; k--) {
        if (C.toStop[k] <= a) continue;
        const g = this.candidate(k, x, false);
        this.cand.self = g.kind === kind && g.id === id && g.row === row;
        Coords.forEachFromAt(this.maps[k], a, this.hit, addIfContains, this);
      }
    }
    c.sink = NoSink;
  }

  /** The closure map k as the candidate: its `from` system is what is added; `element`, the closed element. */
  private candidate(k: number, element: number, self: boolean): CoordinateSystem {
    const s = this.view.systems;
    const from = s.closure.from[k];
    const cand = this.cand;
    cand.system.kind = s.kind[from];
    cand.system.id = s.id[from];
    cand.system.row = s.row[from];
    cand.map = this.maps[k];
    cand.self = self;
    cand.element = element;
    return cand.system;
  }

  /** The last closure row in `[c0, c1)` whose hull starts before `x`; c0 - 1: none. */
  private lastBefore(c0: number, c1: number, x: number): number {
    const toStart = this.view.systems.closure.toStart;
    let [lo, hi] = [c0, c1];
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (toStart[mid] < x) lo = mid + 1;
      else hi = mid;
    }
    return lo - 1;
  }

  /** Appends a run in the candidate's system, merged into an appended entry of that system it touches. */
  append(s: number, t: number) {
    const { ext, system: g, from } = this.cand;
    if (!ext) return;
    for (let n = ext.count - 1; n >= from; n--) {
      if (ext.kind[n] !== g.kind || ext.id[n] !== g.id || ext.row[n] !== g.row) continue;
      if (s > ext.stop[n] || t < ext.start[n]) continue;
      if (s < ext.start[n]) ext.start[n] = s;
      if (t > ext.stop[n]) ext.stop[n] = t;
      return;
    }
    Extents.add(ext, g, s, t);
  }

  widen(ext: Extents, granularity: number): void {
    const G = SelectionGranularity;
    const { set } = this;
    const n = ext.count;
    const sys = WidenSystem;
    let removed = 0;
    for (let i = 0; i < n; i++) {
      [sys.kind, sys.id, sys.row] = [ext.kind[i], ext.id[i], ext.row[i]];
      if (granularity === G.unit || (granularity === G.row && sys.kind === K.row)) {
        [ext.start[i], ext.stop[i]] = Mappings.domainOf(set, sys);
      } else if (granularity === G.copyGroup && sys.kind === K.unit) {
        const members = set.groups.members[set.units.group[sys.id]] ?? [];
        for (const u of members) {
          if (u !== sys.id) Extents.add(ext, { kind: K.unit, id: u, row: -1 }, ext.start[i], ext.stop[i]);
        }
      } else if (granularity === G.column && sys.kind === K.row) {
        const cols = set.alignments[sys.id]?.rows.posMap(sys.row as never);
        const columns = WidenColumns;
        columns.id = sys.id;
        // CF-15: a stretch of adjacent columns is one run, not one extent a column.
        let from = -1;
        let to = -1;
        for (let j = ext.start[i]; cols && j < ext.stop[i] && j < cols.length; j++) {
          const c = cols[j];
          if (c < 0) continue;
          if (c === to) {
            to++;
            continue;
          }
          if (to > from) Extents.add(ext, columns, from, to);
          from = c;
          to = c + 1;
        }
        if (to > from) Extents.add(ext, columns, from, to);
        ext.kind[i] = 255;
        removed++;
      }
    }
    if (!removed) return;
    let w = 0;
    for (let i = 0; i < ext.count; i++) {
      if (ext.kind[i] === 255) continue;
      [ext.kind[w], ext.id[w], ext.row[w], ext.start[w], ext.stop[w]] = [
        ext.kind[i],
        ext.id[i],
        ext.row[i],
        ext.start[i],
        ext.stop[i],
      ];
      w++;
    }
    ext.count = w;
  }
}
