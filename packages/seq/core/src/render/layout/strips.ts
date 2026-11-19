/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Strips as a struct of arrays (_spec/README.md "Render", Section layouts): per strip its transform (its system the
// transform's), its block, its frame's place on the page (a translation: `dx`, `dy`; a polar strip's warp is its
// layout's), its top and height on the page, its first display unit (what its lays are kept by); per band of a strip
// (`strip * bands + band`) its top and height in the strip's frame and what its layers used. Buffers are kept frame to
// frame (`reset` zeroes the count).

import type { AxisTransform } from '../axis.js';
import { grown } from '../draw-lists.js';
import type { Box } from '../layer-type.js';
import { PlacePoint, type Placement } from './types.js';

/** How a strip's frame goes to the page: a translation (tracks, text) or the layout's polar warp (rings). */
export const WarpKind = { translate: 0, polar: 1 } as const;

export class Strips {
  n = 0;
  /** Bands a strip (the section's tracks as shown). */
  bands = 0;
  readonly axes: AxisTransform[] = [];
  /**
   * Each strip's second transform (`strip.y`, a grid's: down the strip's frame, its `x0` the plot's top and its `width`
   * the plot's height), absent in every other layout (_spec/README.md "Render", Section layouts).
   */
  readonly ys: (AxisTransform | undefined)[] = [];
  block = new Int32Array(4);
  /** 1 where the strip is a chip: a one-position block in a fixed slot of screen px (_spec/README.md "Render", Blocks). */
  chip = new Uint8Array(4);
  warp = new Uint8Array(4);
  dx = new Float64Array(4);
  dy = new Float64Array(4);
  /** The strip's top and height on the page. */
  top = new Float64Array(4);
  height = new Float64Array(4);
  /** Its first display unit: its lays are kept by it while it is in view. */
  first = new Float64Array(4);
  /** A grid over blocks: the strip's y block (its index in the section) and that block's first unit, else 0. */
  yBlock = new Int32Array(4);
  yFirst = new Float64Array(4);
  /** Per band of a strip: its top and height in the strip's frame (0: collapsed), and what its layers used. */
  bandY = new Float64Array(8);
  bandH = new Float64Array(8);
  used = new Float64Array(8);
  /**
   * Per band its well's padding above and below the band, in the strip's frame (a ring's radial, scaled with it): the
   * well is the band grown by them (_spec/README.md "Render", Tracks layout). 0 by default.
   */
  padT = new Float64Array(8);
  padB = new Float64Array(8);

  /** Clears the strips for a placement with `bands` bands a strip. */
  reset(bands: number): void {
    this.n = 0;
    this.bands = bands;
    this.axes.length = 0;
    this.ys.length = 0;
  }

  /** Adds a strip over `axis`, translated by `dx`, `dy`; its bands are set after. */
  add(axis: AxisTransform, block: number, dx: number, dy: number, first: number): number {
    const i = this.n++;
    if (i >= this.block.length) {
      const need = i + 1;
      this.block = grown(this.block, need);
      this.chip = grown(this.chip, need);
      this.warp = grown(this.warp, need);
      this.dx = grown(this.dx, need);
      this.dy = grown(this.dy, need);
      this.top = grown(this.top, need);
      this.height = grown(this.height, need);
      this.first = grown(this.first, need);
      this.yBlock = grown(this.yBlock, need);
      this.yFirst = grown(this.yFirst, need);
    }
    const nb = (i + 1) * this.bands;
    if (nb > this.bandY.length) {
      this.bandY = grown(this.bandY, nb);
      this.bandH = grown(this.bandH, nb);
      this.used = grown(this.used, nb);
      this.padT = grown(this.padT, nb);
      this.padB = grown(this.padB, nb);
    }
    this.axes[i] = axis;
    this.ys[i] = undefined;
    this.block[i] = block;
    this.chip[i] = 0;
    this.warp[i] = WarpKind.translate;
    this.dx[i] = dx;
    this.dy[i] = dy;
    this.top[i] = dy;
    this.height[i] = 0;
    this.first[i] = first;
    this.yBlock[i] = 0;
    this.yFirst[i] = 0;
    for (let j = i * this.bands; j < nb; j++) {
      this.bandY[j] = 0;
      this.bandH[j] = 0;
      this.used[j] = 0;
      this.padT[j] = 0;
      this.padB[j] = 0;
    }
    return i;
  }

  /** The strip whose data box holds frame x (strips side by side: a row of blocks, a ring's arcs), or −1. */
  atX(x: number): number {
    for (let i = 0; i < this.n; i++) {
      const a = this.axes[i];
      if (x >= a.x0 && x < a.x0 + a.width) return i;
    }
    return -1;
  }

  /** The strip with an edge nearest frame x (−1 with none). */
  nearestX(x: number): number {
    let best = Infinity;
    let i = -1;
    for (let k = 0; k < this.n; k++) {
      const a = this.axes[k];
      const d = Math.min(Math.abs(x - a.x0), Math.abs(x - a.x0 - a.width));
      if (d < best) {
        best = d;
        i = k;
      }
    }
    return i;
  }

  /** The section's display unit at frame x on strip `i`: its window's edge past its data box. */
  offsetIn(i: number, x: number): number {
    const a = this.axes[i];
    if (x < a.x0) return this.first[i] + a.u0;
    if (x >= a.x0 + a.width) return this.first[i] + a.u1;
    return this.first[i] + a.unitAt(x);
  }

  /** The first band of strip `i` with height (a collapsed band holds nothing), or −1. */
  firstBand(i: number): number {
    for (let j = 0; j < this.bands; j++) if (this.bandH[i * this.bands + j] > 0) return j;
    return -1;
  }

  /**
   * Places a section's display unit on the strip whose window holds it (strips side by side, each from its `first`):
   * its frame x, the middle of its first band with height (else the strip's top), into `out`; false when none holds it.
   */
  placeAcross(offset: number, out: PlacePoint): boolean {
    out.strip = -1;
    out.kind = 'empty';
    // a unit at a strip's end and the next strip's start is the next's (a chip at its own first unit is the chip)
    let at = -1;
    for (let i = 0; i < this.n; i++) {
      const a = this.axes[i];
      const u = offset - this.first[i];
      if (u < a.u0 || u > a.u1) continue;
      if (u < a.u1) {
        at = i;
        break;
      }
      if (at < 0) at = i;
    }
    if (at < 0) return false;
    const band = this.firstBand(at);
    const k = at * this.bands + (band < 0 ? 0 : band);
    out.strip = at;
    out.band = band;
    out.offset = offset;
    out.x = this.axes[at].xOfUnit(offset - this.first[at]);
    out.y = band < 0 ? this.top[at] : this.bandY[k] + this.bandH[k] / 2;
    out.kind = band < 0 ? 'gap' : 'data';
    return true;
  }

  /** The band's box in its strip's frame (the data box's x and width, the band's top and height). */
  bandBox(strip: number, band: number): Box {
    const a = this.axes[strip];
    const j = strip * this.bands + band;
    return { x: a.x0, y: this.bandY[j], width: a.width, height: this.bandH[j] };
  }

  /** The band's well in its strip's frame: its box grown by its padding (`padT`, `padB`). */
  wellBox(strip: number, band: number): Box {
    const a = this.axes[strip];
    const j = strip * this.bands + band;
    return {
      x: a.x0,
      y: this.bandY[j] - this.padT[j],
      width: a.width,
      height: this.bandH[j] + this.padT[j] + this.padB[j],
    };
  }
}

/** A strip-frame x and its display unit (`AxisTransform.inGapOf`, `unitInto`): one scratch, read at once. */
const StripX = new Float64Array(2);

/**
 * The placement of translated strips stacked down the page (tracks, text): `at` finds the strip by its top (a binary
 * search: strips are in reading order, top to bottom), the band by its top, the point in the strip's frame and the
 * display unit under it; `offsetAt` the nearest unit; `place` the strip holding a unit by its window.
 */
export class StripPlacement implements Placement {
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  contentHeight = 0;
  readonly strips = new Strips();
  /**
   * Its strips are lines of blocks (a text of blocks): each on its block's units, from its `first` section unit, and a
   * line's chips side by side on it.
   */
  lineBlocks = false;
  /** Its bands overlay one another (a grid's plot): `at` gives the last band holding the point, the topmost. */
  readonly overlay: boolean = false;

  /** The strip at page y, or the one above it (−1 above every strip). */
  private stripAt(y: number): number {
    const s = this.strips;
    let lo = 0;
    let hi = s.n;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (s.top[mid] <= y) lo = mid + 1;
      else hi = mid;
    }
    return lo - 1;
  }

  at(p: ArrayLike<number>, out: PlacePoint): boolean {
    const x = p[0];
    const y = p[1];
    const b = this.box;
    out.strip = -1;
    out.band = -1;
    out.kind = 'empty';
    out.offset = -1;
    if (y < b.y || y >= b.y + b.height) return false;
    const s = this.strips;
    let i = this.stripAt(y);
    // strips side by side on a line (a text's chips): the one under x
    if (this.lineBlocks) while (i > 0 && s.top[i - 1] === s.top[i] && x - s.dx[i] < s.axes[i].x0) i--;
    if (i < 0 || y >= s.top[i] + s.height[i]) {
      out.kind = 'gap';
      return true;
    }
    const fx = x - s.dx[i];
    const fy = y - s.dy[i];
    out.strip = i;
    out.x = fx;
    out.y = fy;
    for (let j = 0; j < s.bands; j++) {
      const k = i * s.bands + j;
      if (s.bandH[k] > 0 && fy >= s.bandY[k] && fy < s.bandY[k] + s.bandH[k]) {
        out.band = j;
        if (!this.overlay) break;
      }
    }
    const a = s.axes[i];
    // the strip-frame x and its unit through an array: `fx` is computed, and a number computed and passed to (or
    // returned from) a call that is not inlined is boxed (a hit allocates nothing, `Placement.at`)
    StripX[0] = fx;
    if (fx < a.x0) out.kind = 'gutter';
    else if (fx >= a.x0 + a.width) out.kind = 'gutter';
    else if (out.band < 0 || a.inGapOf(StripX)) out.kind = 'gap';
    else out.kind = 'data';
    if (out.kind === 'data') {
      a.unitInto(StripX);
      out.offset = this.lineBlocks ? s.first[i] + StripX[1] - a.u0 : StripX[1];
    }
    return true;
  }

  offsetAt(x: number, y: number): number {
    const s = this.strips;
    if (!s.n) return -1;
    let i = this.stripAt(y);
    if (i < 0) i = 0;
    if (this.lineBlocks) while (i > 0 && s.top[i - 1] === s.top[i] && x - s.dx[i] < s.axes[i].x0) i--;
    const a = s.axes[i];
    const fx = x - s.dx[i];
    // before a line's first cell or past its last (a gutter, the space after it): that cell, so a drag past a line's
    // end reaches its last unit and never the next line's first
    const cx = fx < a.x0 ? a.x0 : fx >= a.x0 + a.width ? a.x0 + a.width - 1e-6 : fx;
    const u = a.unitAt(cx);
    return this.lineBlocks ? s.first[i] + u - a.u0 : u;
  }

  place(offset: number, out: PlacePoint): boolean {
    const s = this.strips;
    out.strip = -1;
    out.kind = 'empty';
    // the strip whose window holds the unit: strips run in reading order (lines of blocks: by their first units)
    const by = this.lineBlocks;
    let lo = 0;
    let hi = s.n;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if ((by ? s.first[mid] : s.axes[mid].u0) <= offset) lo = mid + 1;
      else hi = mid;
    }
    const i = lo - 1;
    const a = i < 0 ? undefined : s.axes[i];
    const unit = a && by ? a.u0 + offset - s.first[i] : offset;
    if (!a || unit > a.u1) return false;
    // on the page: the unit's x, the middle of the strip's first band with height (a collapsed band holds nothing)
    const band = s.firstBand(i);
    out.strip = i;
    out.band = band;
    out.offset = offset;
    out.x = a.xOfUnit(unit) + s.dx[i];
    const k = i * s.bands + (band < 0 ? 0 : band);
    out.y = band < 0 ? s.top[i] : s.bandY[k] + s.bandH[k] / 2 + s.dy[i];
    out.kind = band < 0 ? 'gap' : 'data';
    return true;
  }
}

export { PlacePoint };
