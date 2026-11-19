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
  block = new Int32Array(4);
  warp = new Uint8Array(4);
  dx = new Float64Array(4);
  dy = new Float64Array(4);
  /** The strip's top and height on the page. */
  top = new Float64Array(4);
  height = new Float64Array(4);
  /** Its first display unit: its lays are kept by it while it is in view. */
  first = new Float64Array(4);
  /** Per band of a strip: its top and height in the strip's frame (0: collapsed), and what its layers used. */
  bandY = new Float64Array(8);
  bandH = new Float64Array(8);
  used = new Float64Array(8);

  /** Clears the strips for a placement with `bands` bands a strip. */
  reset(bands: number): void {
    this.n = 0;
    this.bands = bands;
    this.axes.length = 0;
  }

  /** Adds a strip over `axis`, translated by `dx`, `dy`; its bands are set after. */
  add(axis: AxisTransform, block: number, dx: number, dy: number, first: number): number {
    const i = this.n++;
    if (i >= this.block.length) {
      const need = i + 1;
      this.block = grown(this.block, need);
      this.warp = grown(this.warp, need);
      this.dx = grown(this.dx, need);
      this.dy = grown(this.dy, need);
      this.top = grown(this.top, need);
      this.height = grown(this.height, need);
      this.first = grown(this.first, need);
    }
    const nb = (i + 1) * this.bands;
    if (nb > this.bandY.length) {
      this.bandY = grown(this.bandY, nb);
      this.bandH = grown(this.bandH, nb);
      this.used = grown(this.used, nb);
    }
    this.axes[i] = axis;
    this.block[i] = block;
    this.warp[i] = WarpKind.translate;
    this.dx[i] = dx;
    this.dy[i] = dy;
    this.top[i] = dy;
    this.height[i] = 0;
    this.first[i] = first;
    for (let j = i * this.bands; j < nb; j++) {
      this.bandY[j] = 0;
      this.bandH[j] = 0;
      this.used[j] = 0;
    }
    return i;
  }

  /** The band's box in its strip's frame (the data box's x and width, the band's top and height). */
  bandBox(strip: number, band: number): Box {
    const a = this.axes[strip];
    const j = strip * this.bands + band;
    return { x: a.x0, y: this.bandY[j], width: a.width, height: this.bandH[j] };
  }
}

/**
 * The placement of translated strips stacked down the page (tracks, text): `at` finds the strip by its top (a binary
 * search: strips are in reading order, top to bottom), the band by its top, the point in the strip's frame and the
 * display unit under it; `offsetAt` the nearest unit; `place` the strip holding a unit by its window.
 */
export class StripPlacement implements Placement {
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  contentHeight = 0;
  readonly strips = new Strips();

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

  at(x: number, y: number, out: PlacePoint): boolean {
    const b = this.box;
    out.strip = -1;
    out.band = -1;
    out.kind = 'empty';
    out.offset = -1;
    if (y < b.y || y >= b.y + b.height) return false;
    const s = this.strips;
    const i = this.stripAt(y);
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
        break;
      }
    }
    const a = s.axes[i];
    if (fx < a.x0) out.kind = 'gutter';
    else if (fx >= a.x0 + a.width) out.kind = 'gutter';
    else if (out.band < 0 || a.inGap(fx)) out.kind = 'gap';
    else out.kind = 'data';
    if (out.kind === 'data') out.offset = a.unitAt(fx);
    return true;
  }

  offsetAt(x: number, y: number): number {
    const s = this.strips;
    if (!s.n) return -1;
    let i = this.stripAt(y);
    if (i < 0) i = 0;
    const a = s.axes[i];
    const fx = x - s.dx[i];
    const cx = fx < a.x0 ? a.x0 : fx >= a.x0 + a.width ? a.x0 + a.width - 1e-6 : fx;
    // past a line's last cell (its gutter or the space after it): its end
    if (fx >= a.x0 + a.width) return a.u1;
    return a.unitAt(cx);
  }

  place(offset: number, out: PlacePoint): boolean {
    const s = this.strips;
    out.strip = -1;
    out.kind = 'empty';
    // the strip whose window holds the unit: strips run in reading order
    let lo = 0;
    let hi = s.n;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (s.axes[mid].u0 <= offset) lo = mid + 1;
      else hi = mid;
    }
    const i = lo - 1;
    if (i < 0 || offset > s.axes[i].u1) return false;
    const a = s.axes[i];
    out.strip = i;
    out.band = 0;
    out.offset = offset;
    out.x = a.xOfUnit(offset) + s.dx[i];
    out.y = s.top[i];
    out.kind = 'data';
    return true;
  }
}

export { PlacePoint };
