/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Blocks in the tracks layout (_spec/README.md "Render", Blocks): a strip a block in view, each on its own system (its
// transform's), translated by nothing. `row` (the default): the blocks side by side along one axis, one viewport over
// them laid end to end, a gap between two (shrinking when zoomed out so positions keep three quarters of the width,
// the gaps in view sharing the rest, EC-102); `stack`: a block a row of the section's tracks, aligned at their first
// position, one zoom (EC-102, EC-115). Each block's label in a row above its tracks. Its lays are kept by the block
// and its first unit, so a pan keeps them. A one-position block is a chip: a fixed slot its label's width, never a
// position's share (zoom does not scale it), where it stands in a row, in a stack on one row the chips share.

import { piecesTotal, type ResolvedSection } from '@molstar/seqx-core/view';
import { AxisTransform } from '../axis.js';
import { fontSize, type RenderConfig } from '../config.js';
import { Anchor, type PassName, type PassWriter } from '../draw-lists.js';
import type { TextMetrics } from '../text.js';
import type { SectionLayout } from '../tracks.js';
import { PlacePoint, StripPlacement } from './strips.js';
import type { PlaceContext } from './types.js';

// the anchor a number read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;

/**
 * Fixed sizes (CSS px): the gap between two blocks in a row at most (a boundary that reads), the share of the width
 * the gaps in view may take (a quarter: positions keep three quarters, however many blocks), a chip's padding either
 * side of its label (the text chip's rule).
 */
export const BlockMetrics = { gap: 12, gapShare: 0.25, chipPad: 5 } as const;

/** Whether a block is a chip: one position (a ligand, an ion, a glycan standing as one). */
export function isChip(block: { readonly axis: { readonly ranges?: readonly number[] } }): boolean {
  return !!block.axis.ranges && piecesTotal(block.axis.ranges) === 1;
}

/** A chip's width, CSS px: its label (its `×N` too) in the label font, `chipPad` either side; zoom never scales it. */
export function chipWidth(text: TextMetrics, config: RenderConfig, label: string): number {
  return text.measure(config.labelFont, label) + 2 * BlockMetrics.chipPad;
}

/** A section of blocks placed: its strips by x (a row) or by y (a stack), each strip's units its block's. */
export class BlockPlacement extends StripPlacement {
  stacked = false;
  /** The label row above each block's tracks: its height and the labels' baseline in it. */
  labelRow = 0;
  /** Each strip's block's label (and the block's first unit: `strips.first`). */
  labels: string[] = [];

  /** The strip at a page point: by x in a row, by its rows' extent in a stack (the chips' row by x too); −1 when none. */
  private stripOf(x: number, y: number): number {
    const s = this.strips;
    if (!this.stacked) return s.atX(x);
    let row = -1;
    for (let i = 0; i < s.n; i++) {
      if (y < s.top[i] || y >= s.top[i] + s.height[i]) continue;
      if (!s.chip[i]) return i;
      if (x >= s.axes[i].x0 && x < s.axes[i].x0 + s.axes[i].width) return i;
      if (row < 0) row = i;
    }
    return row;
  }

  override at(p: ArrayLike<number>, out: PlacePoint): boolean {
    const x = p[0];
    const y = p[1];
    const b = this.box;
    out.strip = -1;
    out.band = -1;
    out.kind = 'empty';
    out.offset = -1;
    if (y < b.y || y >= b.y + b.height) return false;
    const s = this.strips;
    const i = this.stripOf(x, y);
    const first = s.n ? s.axes[0].x0 : 0;
    if (i < 0) {
      out.kind = x < first ? 'gutter' : 'gap';
      return true;
    }
    const a = s.axes[i];
    out.strip = i;
    out.x = x;
    out.y = y;
    for (let j = 0; j < s.bands; j++) {
      const k = i * s.bands + j;
      if (s.bandH[k] > 0 && y >= s.bandY[k] && y < s.bandY[k] + s.bandH[k]) {
        out.band = j;
        break;
      }
    }
    if (x < a.x0) out.kind = 'gutter';
    else if (x >= a.x0 + a.width) out.kind = this.stacked ? 'gap' : 'gutter';
    else if (out.band < 0) out.kind = 'gap';
    else {
      out.kind = 'data';
      out.offset = s.first[i] + a.unitAt(x);
    }
    return true;
  }

  override offsetAt(x: number, y: number): number {
    const s = this.strips;
    if (!s.n) return -1;
    let i = this.stripOf(x, y);
    if (i < 0 && !this.stacked) i = s.nearestX(x);
    else if (i < 0) {
      // a stack's nearest row by y
      let best = Infinity;
      for (let k = 0; k < s.n; k++) {
        const d = Math.min(Math.abs(y - s.top[k]), Math.abs(y - s.top[k] - s.height[k]));
        if (d < best) [best, i] = [d, k];
      }
    }
    return s.offsetIn(i, x);
  }

  override place(offset: number, out: PlacePoint): boolean {
    return this.strips.placeAcross(offset, out);
  }
}

/** Places a section's blocks in a row or a stack over the section's window (`ctx.axis`'s). */
export function placeBlocks(section: ResolvedSection, ctx: PlaceContext, prev?: BlockPlacement): BlockPlacement {
  const p = prev ?? new BlockPlacement();
  const blocks = section.blocks!;
  const { config, axis } = ctx;
  const strips = p.strips;
  const nb = ctx.bands.length;
  strips.reset(nb);
  p.labels.length = 0;
  p.stacked = section.params.arrange === 'stack';
  p.labelRow = fontSize(config.labelFont) + 4;
  const top = ctx.box.y;
  const u0 = axis.u0;
  const u1 = axis.u1;
  const x0 = axis.x0;
  const width = axis.width;
  // a row's gaps inside the window take their px from the width (as a text's groups do), and so do its chips in view:
  // the row's x is not linear over a chip, as over a gap (EC-102)
  let inner = 0;
  let chipPx = 0;
  let chipUnits = 0;
  for (let b = 0; b < blocks.length; b++) {
    const o = blocks[b].offset ?? 0;
    if (b > 0 && o > u0 && o < u1) inner++;
    if (p.stacked || !isChip(blocks[b]) || o + 1 <= u0 || o >= u1) continue;
    chipPx += chipWidth(ctx.text, config, blocks[b].label);
    chipUnits += Math.min(u1, o + 1) - Math.max(u0, o);
  }
  // the gaps in view share a quarter of the width, at most 12 px each
  const gap = p.stacked ? 0 : Math.min(BlockMetrics.gap, (BlockMetrics.gapShare * width) / Math.max(1, inner));
  const scale = Math.max(width - inner * gap - chipPx, 1e-9) / Math.max(u1 - u0 - chipUnits, 1e-9);
  // the rows' height: the label row, then the tracks
  let rows = p.labelRow;
  for (let j = 0; j < nb; j++) rows += (j ? config.trackGap : 0) + ctx.heights[j];
  let y = top;
  let gaps = 0;
  // a row: the chips before a block, their px and their units in the window
  let chipsPx = 0;
  let chipsUnits = 0;
  // a stack: the other blocks a row each, then the chips' one row, side by side from the first position's x
  let stackRows = 0;
  let chipRow = -1;
  let chipX = x0;
  for (let round = 0; round < (p.stacked ? 2 : 1); round++)
    for (let b = 0; b < blocks.length; b++) {
      const block = blocks[b];
      const chip = isChip(block);
      if (p.stacked && chip !== (round === 1)) continue;
      const o = p.stacked ? 0 : (block.offset ?? 0);
      const len = block.axis.ranges ? piecesTotal(block.axis.ranges) : 0;
      if (b > 0 && !p.stacked && o > u0 && o < u1) gaps++;
      let rowTop = y;
      if (p.stacked && (!chip || chipRow < 0)) {
        if (stackRows++) y += config.sectionGap;
        rowTop = y;
        y += rows;
        if (chip) chipRow = rowTop;
      } else if (p.stacked) rowTop = chipRow;
      const vs = Math.max(u0, o);
      const ve = Math.min(u1, o + len);
      // a stack's chips are on their row whatever the window
      if (!block.axis.ranges || (ve <= vs && !(chip && p.stacked))) continue;
      let a: AxisTransform;
      if (chip) {
        // a fixed slot its label's width over its one position whole: zoom never scales it
        const w = chipWidth(ctx.text, config, block.label);
        const x = p.stacked ? chipX : x0 + (vs - u0 - chipsUnits) * scale + gaps * gap + chipsPx;
        a = AxisTransform.on(block.axis, block.axis.ranges, x, w, ctx.dpr, { start: 0, end: 1 });
        chipX += w + BlockMetrics.gap;
        chipsPx += w;
        chipsUnits += ve - vs;
      } else
        a = AxisTransform.on(
          block.axis,
          block.axis.ranges,
          x0 + (vs - u0 - chipsUnits) * scale + gaps * gap + chipsPx,
          (ve - vs) * scale,
          ctx.dpr,
          { start: vs - o, end: ve - o },
        );
      const i = strips.add(a, b, 0, 0, o);
      if (chip) strips.chip[i] = 1;
      p.labels[i] = block.label;
      strips.top[i] = rowTop;
      strips.height[i] = rows;
      let by = rowTop + p.labelRow;
      for (let j = 0; j < nb; j++) {
        if (j) by += config.trackGap;
        strips.bandY[i * nb + j] = by;
        strips.bandH[i * nb + j] = ctx.heights[j];
        by += ctx.heights[j];
      }
    }
  const height = p.stacked ? y - top : rows;
  p.contentHeight = height;
  p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height };
  for (let i = 0; i < strips.n; i++) ctx.layStrip(p, i);
  return p;
}

/** A block's label above its tracks (the labels pass), faded at its strip's width; a stack's rows filled as tracks. */
export function blockChrome(
  s: SectionLayout,
  p: BlockPlacement,
  pass: PassName,
  out: PassWriter,
  strip: number,
  band: number,
): void {
  const strips = p.strips;
  if (pass === 'data' && band >= 0 && strip > 0 && p.stacked) {
    const k = strip * strips.bands + band;
    out.fill(
      s.tracks[band]?.box.x ?? 0,
      strips.bandY[k],
      s.tracks[band]?.box.width ?? 0,
      strips.bandH[k],
      s.tokens.colors['seq-track'],
    );
    return;
  }
  if (pass !== 'labels' || band >= 0) return;
  const g = out.args;
  const size = fontSize(s.config.labelFont);
  for (let i = 0; i < strips.n; i++) {
    // a chip's label is in its letters' pill
    if (strips.chip[i] && hasPill(s, i)) continue;
    const a = strips.axes[i];
    g[0] = a.x0;
    g[1] = strips.top[i] + size;
    g[2] = 0;
    g[3] = Math.max(0, a.width);
    g[4] = s.tokens.colors['msui-text-muted'];
    out.textArgs(p.labels[i], AnchorStart, s.config.labelFont);
  }
}

/** Whether a chip's strip draws its label in a pill: a letters layer laid on one of its bands. */
function hasPill(s: SectionLayout, strip: number): boolean {
  const nb = s.placement.strips.bands;
  for (let j = 0; j < nb; j++) {
    const lays = s.lays[strip * nb + j];
    for (let k = 0; lays && k < lays.length; k++) if (lays[k].type.name === 'letters') return true;
  }
  return false;
}
