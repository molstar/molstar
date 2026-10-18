/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's projector over the projection core (seqspec:interaction#channels; seqspec:model#projection,
// seqspec:libraries#translate; CP1b): one per surface, plus one for edge-rate callers (focus, `toLoci`,
// `selection.close`). It holds one core `Projection` per (view, set), where routes live, so a new view or set makes a
// new one; channels' loci and `Coordinates` become `Extents` in Seq* systems (a host's external system through the maps
// registered with it in `ctx.systems`, then on the usual route), close, and project into a channel's `MarkSet` block by
// block. A ligand in the `position` form lights its block's one shown position for every atom. Nothing of it is module
// state, so a disposed viewer frees its set; it is not re-entrant (its sinks never call out).

import { type Closure, type Projection, projectionOf, type ResolvedView, type RunSink } from '@molstar/seq-core/format';
import {
  type CoordinateSystem,
  Coords,
  Extents,
  CoordinateSystemKind as K,
  type Loci,
  type SeqLoci,
  type SeqSet,
} from '@molstar/seq-core/model';
import { Coordinates, type Translator } from '@molstar/seq-core/translate';
import { HoverKind } from './hit.js';
import type { HoverState } from './hover.js';
import { MarkBit, type MarkSet, type RowMarks } from './marks.js';

/** A ligand in the `position` form: every atom lights the block's one position. */
class LigandSink implements RunSink {
  out: RunSink | null = null;
  at = 0;
  add() {
    this.out?.add(this.at, this.at + 1);
  }
}
/**
 * What a block's runs are clipped to: its piece, or across a cyclic origin (stop <= start: two runs, or the whole from
 * an `origin`) its domain, whose positions outside the piece no draw shows (the run converters map none of them).
 */
const pieceLo = (v: ResolvedView, b: number) =>
  v.blocks.stop[b] <= v.blocks.start[b] ? v.blocks.base[b] : v.blocks.start[b];
const pieceHi = (v: ResolvedView, b: number) =>
  v.blocks.stop[b] <= v.blocks.start[b] ? v.blocks.base[b] + v.blocks.length[b] : v.blocks.stop[b];

/** The hull of what lights one block, clipped to its piece. */
class HullSink implements RunSink {
  lo = 0;
  hi = 0;
  min = 0;
  max = 0;
  /** On a whole cyclic unit, the runs too, for their shortest arc (`arc`). */
  cyclic = false;
  n = 0;
  runs = new Int32Array(16);
  private readonly arcOut = new Int32Array(4);
  add(s: number, t: number) {
    if (s < this.min) s = this.min;
    if (t > this.max) t = this.max;
    if (s >= t) return;
    if (s < this.lo) this.lo = s;
    if (t > this.hi) this.hi = t;
    if (!this.cyclic) return;
    if (2 * this.n + 2 > this.runs.length) {
      const r = new Int32Array(2 * this.runs.length);
      r.set(this.runs);
      this.runs = r;
    }
    // Kept sorted by start (few runs: a focus target's).
    let i = this.n++;
    for (; i > 0 && this.runs[2 * i - 2] > s; i--) {
      this.runs[2 * i] = this.runs[2 * i - 2];
      this.runs[2 * i + 1] = this.runs[2 * i - 1];
    }
    this.runs[2 * i] = s;
    this.runs[2 * i + 1] = t;
  }
  /**
   * The runs' shortest arc on a cyclic system of `length` (seqspec:model#cyclic-systems) as lo, hi: lo > hi across the
   * origin.
   */
  arc(length: number) {
    const o = this.arcOut;
    if (Coords.cyclicHull(this.runs, this.n, length, o) !== 2) return;
    this.lo = o[0];
    this.hi = o[3];
  }
}
/** Merged runs of one block: a row's residues on its columns, before they become rectangles. */
class RowSink implements RunSink {
  count = 0;
  start = new Int32Array(8);
  stop = new Int32Array(8);
  add(s: number, t: number) {
    const n = this.count;
    if (n > 0 && s <= this.stop[n - 1] && t >= this.start[n - 1]) {
      this.start[n - 1] = Math.min(s, this.start[n - 1]);
      this.stop[n - 1] = Math.max(t, this.stop[n - 1]);
      return;
    }
    if (n === this.start.length) {
      const a = new Int32Array(2 * n);
      const b = new Int32Array(2 * n);
      a.set(this.start);
      b.set(this.stop);
      this.start = a;
      this.stop = b;
    }
    this.start[n] = s;
    this.stop[n] = t;
    this.count = n + 1;
  }
}

export const isCoordinates = (t: Loci | Coordinates): t is Coordinates => Coordinates.is(t);
export const ofRuns = Coordinates.ofRuns;
export type { Coordinates };

/** Per block lit: (block in the section, start, stop) triples in `data`, `count` of them; reused, grown when full. */
export interface BlockHulls {
  count: number;
  data: Int32Array;
}
export const BlockHulls = { create: (): BlockHulls => ({ count: 0, data: new Int32Array(12) }) };

/** Read once: a member read in a hover-rate path allocates under the test runner's module wrapper. */
const HEADER = HoverKind.header;
const SEQ = HoverKind.seq;
const UNIT = K.unit;

export class Projector {
  /** The host's systems (`ctx.systems`): `Coordinates` keys and external systems' maps; null: loci only. */
  systems: Translator | null = null;
  private view: ResolvedView | null = null;
  private set: SeqSet | null = null;
  private core: Projection | null = null;
  private readonly ligand = new LigandSink();
  private readonly hull = new HullSink();
  private readonly row = new RowSink();
  private readonly one = Extents.create(1);
  private readonly single = Extents.create(1);
  /** A row entry's system, passed to `Extents.add` (which copies it). */
  private readonly sys: CoordinateSystem = { kind: 0, id: -1, row: -1 };

  /** Bumped by each new view or set (`update`): what a reader resolved through the projector before is stale. */
  epoch = 0;

  /** The view and set channels project onto: a different one drops the core and its routes. */
  update(view: ResolvedView | null, set: SeqSet | null): void {
    if (view === this.view && set === this.set) return;
    this.view = view;
    this.set = set;
    this.core = null;
    this.epoch++;
  }

  /** The core over the current view and set, made on first use; null: none shown. */
  get projection(): Projection | null {
    if (!this.core && this.view && this.set) this.core = projectionOf(this.view, this.set);
    return this.core;
  }

  /**
   * A channel's target appended to `out` as entries in Seq* systems (the pipeline's first step); the entries added.
   * Loci of another set light nothing until carried; `Coordinates` read their keys through `systems` (an unregistered
   * key: nothing meanwhile), an external system's runs going through each map registered with it. Allocation-free once
   * `out` has room.
   */
  fromLoci(target: Loci | Coordinates, out: Extents): number {
    if (!isCoordinates(target)) return target.set === this.set && this.set ? Extents.addLoci(out, target) : 0;
    return this.systems ? Coordinates.toExtents(this.systems, target, out, this.one) : 0;
  }

  /**
   * The hovered element appended to `out`: its location's index in `h.system` (a residue in its unit or molecule, a
   * row's residue, a gap's or a column's column), an item's rows; nothing over chrome or a set not shown.
   * Allocation-free once `out` has room.
   */
  fromHover(h: HoverState, out: Extents): void {
    const set = this.set;
    let at = -1;
    if (!set) return;
    if (h.kind === HoverKind.seq && h.seq.set === set) at = h.seq.index;
    else if (h.kind === HoverKind.row && h.row.set === set) at = h.row.index;
    else if (h.kind === HoverKind.column && h.column.set === set) at = h.column.column;
    else if (h.kind === HoverKind.aligned && h.aligned.set === set) {
      const { alignment, row, column } = h.aligned;
      at = h.system.kind === K.row ? (set.alignments[alignment]?.rows.colToPos(row, column) ?? -1) : column;
    } else if (h.kind === HoverKind.item && h.item.set === set)
      Extents.addItem(out, set, h.item.annotation, h.item.item);
    // A header: its runs, when the `headerHover` behavior lights it (header-hover.ts).
    else if (h.kind === HEADER && h.header.light) {
      const c = h.header.extents;
      const s = this.sys;
      for (let i = 0; i < c.count; i++) {
        s.kind = c.kind[i];
        s.id = c.id[i];
        s.row = c.row[i];
        Extents.add(out, s, c.start[i], c.stop[i]);
      }
    }
    // A streamed item (ST4b): its range on the unit, as an item's extent lights.
    if (h.kind === HoverKind.seq && h.span[1] > h.span[0]) Extents.add(out, h.system, h.span[0], h.span[1]);
    else if (at >= 0) Extents.add(out, h.system, at, at + 1);
    // A dot plot's cell: its partner under y too (seqspec:interaction#hover).
    if (h.kind === SEQ && h.pairUnit >= 0 && h.seq.set === set) {
      const s = this.sys;
      s.kind = UNIT;
      s.id = h.pairUnit;
      s.row = -1;
      Extents.add(out, s, h.pairIndex, h.pairIndex + 1);
    }
  }

  /** Closure, appending in place: the elements of the view's other systems that contain each element. */
  close(ext: Extents): void {
    this.projection?.close(ext);
  }

  /**
   * Projection and marks: the extents onto every block of section `section`, into `out` (reset first): per block its
   * runs merged and clipped to the block's piece, a ligand in the `position` form at its one position.
   */
  project(ext: Extents, section: number, out: MarkSet): void {
    const v = this.view;
    const p = this.projection;
    const n = v && section >= 0 && section < v.sections.count ? v.sections.blockCount[section] : 0;
    out.reset(n);
    if (!v || !p) return;
    const b0 = v.sections.blockStart[section];
    for (let b = 0; b < n; b++) {
      out.begin(b);
      this.onBlock(p, ext, b0 + b, out);
      out.close(pieceLo(v, b0 + b), pieceHi(v, b0 + b));
    }
  }

  /**
   * Focus and scroll-to: the hull of the extents' projection (no closure; links and rows included) on every block of
   * the section that it lights, in axis order, into `out` as (block in the section, start, stop) in each block's
   * coordinate system, clipped to its piece; the count of blocks lit.
   */
  extents(ext: Extents, section: number, out: BlockHulls): number {
    out.count = 0;
    const v = this.view;
    const p = this.projection;
    if (!v || !p || section < 0 || section >= v.sections.count) return 0;
    const h = this.hull;
    const b0 = v.sections.blockStart[section];
    for (let b = 0; b < v.sections.blockCount[section]; b++) {
      h.lo = 0x7fffffff;
      h.hi = -0x80000000;
      h.min = pieceLo(v, b0 + b);
      h.max = pieceHi(v, b0 + b);
      h.cyclic = this.wholeCyclic(v, b0 + b);
      h.n = 0;
      this.onBlock(p, ext, b0 + b, h);
      if (h.hi <= h.lo) continue;
      // A whole cyclic unit: the shortest arc, so a feature across the origin focuses as one (start > stop).
      if (h.cyclic) h.arc(v.blocks.length[b0 + b]);
      if (3 * out.count + 3 > out.data.length) {
        const d = new Int32Array(2 * out.data.length + 3);
        d.set(out.data);
        out.data = d;
      }
      const i = 3 * out.count++;
      out.data[i] = b;
      out.data[i + 1] = h.lo;
      out.data[i + 2] = h.hi;
    }
    return out.count;
  }

  /**
   * In a section of an alignment's rows: each row entry's residues as a rectangle of its row by the hull of the columns
   * they land on (gaps included), into `out` under `bit | MarkBit.row` (row marks: the row is lit whole). Other entries
   * light positions, through `project`.
   */
  rows(ext: Extents, section: number, bit: number, out: RowMarks): void {
    const v = this.view;
    const p = this.projection;
    if (!v || !p || section < 0 || section >= v.sections.count) return;
    const b0 = v.sections.blockStart[section];
    const one = this.single;
    for (let b = 0; b < v.sections.blockCount[section]; b++) {
      const a = v.blocks.alignment[b0 + b];
      for (let i = 0; a >= 0 && i < ext.count; i++) {
        if (ext.kind[i] !== K.row || ext.id[i] !== a) continue;
        one.count = 0;
        this.sys.kind = K.row;
        this.sys.id = a;
        this.sys.row = ext.row[i];
        Extents.add(one, this.sys, ext.start[i], ext.stop[i]);
        this.row.count = 0;
        p.project(one, b0 + b, this.row);
        // One rectangle over the runs' hull, gaps included, flagged as a row lit whole.
        let s0 = 0x7fffffff;
        let t0 = -1;
        for (let r = 0; r < this.row.count; r++) {
          s0 = Math.min(s0, this.row.start[r]);
          t0 = Math.max(t0, this.row.stop[r]);
        }
        out.add(bit | MarkBit.row, b, ext.row[i], s0, t0);
      }
    }
  }

  /**
   * Edge rate, allocates (`selection.close`, the select behavior): seq loci closed by `closure`, each molecule element
   * in every alive unit of its molecule. Other kinds are returned as given: rows, cells and columns close only through
   * span-1 maps, which add nothing, and an item's closure is none (its rows; the channels' table).
   */
  closeLoci(loci: Loci, closure: Closure): Loci {
    if (closure === 'element' || loci.kind !== 'seq-loci' || loci.set !== this.set || !this.projection) return loci;
    const ext = Extents.create(8);
    Extents.addLoci(ext, loci);
    this.projection.close(ext);
    return this.toSeqLoci(ext) ?? loci;
  }

  /**
   * The unit and molecule entries of `ext` as one seq loci, a molecule entry in every alive unit of its molecule (what
   * `hover.toLoci` and a click give: every copy); null: none. Edge rate: per unit its runs, never a Set of positions,
   * so a chromosome-long drag closes in O(runs).
   */
  toSeqLoci(ext: Extents): SeqLoci | null {
    return this.set ? Coordinates.toSeqLoci(this.set, ext) : null;
  }

  /**
   * `Coordinates` as loci in the shown set (edge rate, allocates; the selection's verbs): its units and molecules as
   * one seq loci, rows as one row loci and columns as one column loci per alignment; unregistered keys give nothing.
   */
  lociOf(c: Coordinates): Loci[] {
    const set = this.set;
    const ext = Extents.create(8);
    if (!set || this.fromLoci(c, ext) === 0) return [];
    return Coordinates.toLoci(set, ext);
  }

  /** A block over a unit over all of a cyclic molecule. */
  private wholeCyclic(v: ResolvedView, b: number): boolean {
    const u = v.blocks.unit[b];
    const m = u >= 0 && this.set ? this.set.molecules[this.set.units.molecule[u]] : null;
    return !!m?.cyclic && v.blocks.base[b] === 0 && v.blocks.length[b] === m.length;
  }

  private onBlock(p: Projection, ext: Extents, block: number, out: RunSink) {
    const v = this.view as ResolvedView;
    if (!v.blocks.ligand[block]) return p.project(ext, block, out);
    this.ligand.out = out;
    this.ligand.at = v.blocks.start[block];
    p.project(ext, block, this.ligand);
    this.ligand.out = null;
  }
}
