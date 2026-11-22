/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge laid out (_spec/README.md "Render", Bridges): its band, the room between its two sides in the stack, and
// each side as an end reads it: the y of its edge facing the band and the page x of a range on a system it shows
// (every place of a block showing it, a reversed block's first position the larger x; past the window toward its true
// x, clamped one data box past the side's edge; a part in a hidden stretch at its break's hairline). Part of the bridge
// chunk: loaded when a view has a bridge.

import { blockUnitOf, type ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from '../axis.js';
import type { Box } from '../layer-type.js';
import type { SectionLayout } from '../tracks-types.js';
import { BlockPlacement } from './blocks.js';

/** Each system's blocks (indices, in order) of a resolved section: built once a section. */
const SystemBlocks = new WeakMap<ResolvedSection, ReadonlyMap<string, Int32Array>>();

function systemBlocks(section: ResolvedSection): ReadonlyMap<string, Int32Array> {
  let m = SystemBlocks.get(section);
  if (m) return m;
  const lists = new Map<string, number[]>();
  const blocks = section.blocks ?? [];
  for (let b = 0; b < blocks.length; b++)
    for (const s of new Set(blocks[b].systems.length ? blocks[b].systems : [blocks[b].axis.system])) {
      let l = lists.get(s);
      if (!l) lists.set(s, (l = []));
      l.push(b);
    }
  m = new Map([...lists].map(([s, l]) => [s, Int32Array.from(l)]));
  SystemBlocks.set(section, m);
  return m;
}

const Assemblies = new WeakMap<ResolvedSection, ReadonlySet<string>>();

/** The assemblies a section's blocks' sequences are of: built once a section. */
function assembliesOf(section: ResolvedSection, blocks: ReadonlyMap<string, Int32Array>): ReadonlySet<string> {
  let a = Assemblies.get(section);
  if (a) return a;
  a = new Set([...blocks.keys()].filter((s) => s.lastIndexOf('/(') > 0).map((s) => s.slice(0, s.lastIndexOf('/('))));
  Assemblies.set(section, a);
  return a;
}

/** The x a unit row writes into (a typed array's slot: a number returned from a call not inlined is boxed). */
const RowX = new Float64Array(2);

/**
 * One side of a bridge as its ribbons read it (_spec/README.md "Render", Bridges: `BridgeSide`): its section laid out,
 * the y of its edge facing the band (its bottom when above, its top when below), its data box.
 */
export class BridgeSide {
  layout!: SectionLayout;
  section!: ResolvedSection;
  /** The y of its edge facing the band. */
  edge = 0;
  /** Whether it is the side above the band. */
  above = false;
  /** Its data box's x and width (an end past it is clamped one width further). */
  x0 = 0;
  width = 1;
  /** Its device pixels a CSS px. */
  dpr = 1;
  private blocks: ReadonlyMap<string, Int32Array> | undefined;
  /** The assemblies its blocks' sequences are of (by the address before a sequence's query). */
  private assemblies: ReadonlySet<string> | undefined;
  private row: BlockPlacement | undefined;
  private axis: AxisTransform | undefined;

  set(layout: SectionLayout, above: boolean, dpr: number): this {
    this.layout = layout;
    this.section = layout.section;
    this.above = above;
    this.edge = above ? layout.y + layout.height : layout.y;
    this.dpr = dpr;
    const a = layout.axis;
    this.axis = a;
    const p = layout.placement;
    this.row = p instanceof BlockPlacement && !p.stacked ? p : undefined;
    this.blocks = layout.section.blocks ? systemBlocks(layout.section) : undefined;
    this.assemblies = this.blocks ? assembliesOf(layout.section, this.blocks) : undefined;
    const box = this.row ? { x: this.row.row.x0, width: this.row.row.width } : a ? { x: a.x0, width: a.width } : p.box;
    this.x0 = box.x;
    this.width = box.width;
    return this;
  }

  /** Whether the side shows a system (a block of it, or its own axis's). */
  shows(system: string): boolean {
    return this.blocks ? this.blocks.has(system) : this.section.axis.system === system;
  }

  /**
   * Whether an end on a system is the side's: shown, or a sequence of an assembly the section shows (a block of a
   * sequence beside it, `<assembly>/(sequence …)`) that its `include` leaves out (Defaults, 4: a stub, never unplaced).
   */
  holds(system: string): boolean {
    if (this.shows(system)) return true;
    const at = system.lastIndexOf('/(');
    return at > 0 && !!this.assemblies?.has(system.slice(0, at));
  }

  /** Whether its axis is laid out (a pending section has none: nothing is placed on it). */
  get ready(): boolean {
    return this.axis !== undefined;
  }

  /**
   * The places of `[start, end)` on `system` (`shows` it): each place's x of its first position and of its end into
   * `out` (`2k`, `2k + 1`; a reversed block's first the larger) and its block into `blocks` (−1 without blocks), at most
   * `max`; their count, 0 when no block shows any of it (a stub's end). With `pos`, each place's positions `[s, e)`
   * (the part of `[start, end)` it shows).
   */
  extentInto(
    system: string,
    start: number,
    end: number,
    out: Float64Array,
    blocks: Int32Array,
    max: number,
    pos?: Float64Array,
  ): number {
    if (!this.axis) return 0;
    if (!this.blocks) {
      const a = this.axis;
      const x0 = xAt(a, start, true);
      const x1 = xAt(a, end, false);
      if (Number.isNaN(x0) || Number.isNaN(x1)) return 0;
      out[0] = this.clamp(x0 + this.dx());
      out[1] = this.clamp(x1 + this.dx());
      blocks[0] = -1;
      if (pos) {
        pos[0] = start;
        pos[1] = end;
      }
      return 1;
    }
    const list = this.blocks.get(system);
    if (!list) return 0;
    const bs = this.section.blocks!;
    let n = 0;
    for (let i = 0; i < list.length && n < max; i++) {
      const b = bs[list[i]];
      const r = b.axis.ranges;
      if (!r || b.offset === undefined) continue;
      for (let k = 0; k + 1 < r.length && n < max; k += 2) {
        const s = start > r[k] ? start : r[k];
        const e = end < r[k + 1] ? end : r[k + 1];
        if (e <= s) continue;
        const u0 = blockUnitOf(b, s);
        const u1 = blockUnitOf(b, e);
        if (u0 === undefined || u1 === undefined) continue;
        out[2 * n] = this.clamp(this.unitX(u0));
        out[2 * n + 1] = this.clamp(this.unitX(u1));
        if (pos) {
          pos[2 * n] = s;
          pos[2 * n + 1] = e;
        }
        blocks[n] = list[i];
        n++;
      }
    }
    return n;
  }

  /** The page x of a section unit (a row of blocks: as last laid, past the window on at its scale). */
  private unitX(u: number): number {
    if (this.row) {
      this.row.row.xInto(u, RowX, 0);
      return RowX[0];
    }
    return this.axis!.xOfUnit(u) + this.dx();
  }

  private dx(): number {
    const s = this.layout.placement.strips;
    return s.n ? s.dx[0] : 0;
  }

  /** An end past the window toward its true x, clamped one data box past the side's edge (EC-129: continuous on pan). */
  clamp(x: number): number {
    const lo = this.x0 - this.width;
    const hi = this.x0 + 2 * this.width;
    return x < lo ? lo : x > hi ? hi : x;
  }
}

/**
 * The x of the boundary before position `p` on a section's transform; in a hidden stretch, its break's hairline (the
 * next shown position for a start, the last for an end); NaN when no piece is on that side of it.
 */
function xAt(a: AxisTransform, p: number, start: boolean): number {
  const x = a.x(p);
  if (!Number.isNaN(x)) return x;
  const r = a.pieces;
  let best = Number.NaN;
  let bestD = Number.POSITIVE_INFINITY;
  for (let k = 0; k + 1 < r.length; k += 2) {
    const b = start ? r[k] : r[k + 1];
    const d = start ? b - p : p - b;
    if (d >= 0 && d < bestD) {
      bestD = d;
      best = b;
    }
  }
  return Number.isNaN(best) ? best : a.x(best);
}

/** A bridge placed: its band on the page and its two sides, `a` and `b` as the bridge names them. */
export class BridgePlacement {
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  readonly a = new BridgeSide();
  readonly b = new BridgeSide();

  /** Places the band between its sides: `a` above when `aAbove`, the band from the upper side's foot to the lower's top. */
  place(a: SectionLayout, b: SectionLayout, aAbove: boolean, dpr: number): this {
    this.a.set(a, aAbove, dpr);
    this.b.set(b, !aAbove, dpr);
    const top = aAbove ? this.a.edge : this.b.edge;
    const bottom = aAbove ? this.b.edge : this.a.edge;
    const x0 = Math.min(this.a.x0, this.b.x0);
    const x1 = Math.max(this.a.x0 + this.a.width, this.b.x0 + this.b.width);
    this.box = { x: x0, y: top, width: Math.max(0, x1 - x0), height: Math.max(0, bottom - top) };
    return this;
  }

  /** Whether a page point is in the band. */
  holds(x: number, y: number): boolean {
    const b = this.box;
    return y >= b.y && y < b.y + b.height && x >= b.x && x < b.x + b.width;
  }
}
