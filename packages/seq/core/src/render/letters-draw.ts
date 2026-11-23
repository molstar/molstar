/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The letters layer's passes (_spec/README.md "Render", Layer types: letters): the data pass (bars, placeholders, gap
// marks and carets, cells or ticks, the letters or three-letter codes on them, a chip's pill) and the hover pass, each
// in a function of its own called through a table (`drawPass`: nothing allocated once warm).

import { Black, contrastRgba, White } from '@molstar/seqx-core/theme/index';
import type { AxisTransform } from './axis.js';
import { Anchor, NoColor, type PassWriter, passDraws } from './draw-lists.js';
import type { Emphasis } from './layer-type.js';
import type { LettersLay } from './letters.js';
import { drawFailed, drawSpans, hoverCell } from './marks.js';
import { inkOn, type SeqTokens } from './tokens.js';

const AnchorMiddle: number = Anchor.middle;

/** Three-letter codes by one-letter code (IUPAC), the stop as `Ter`. */
const ThreeLetter: Readonly<Record<string, string>> = {
  A: 'Ala',
  R: 'Arg',
  N: 'Asn',
  D: 'Asp',
  C: 'Cys',
  Q: 'Gln',
  E: 'Glu',
  G: 'Gly',
  H: 'His',
  I: 'Ile',
  L: 'Leu',
  K: 'Lys',
  M: 'Met',
  F: 'Phe',
  P: 'Pro',
  S: 'Ser',
  T: 'Thr',
  W: 'Trp',
  Y: 'Tyr',
  V: 'Val',
  U: 'Sec',
  O: 'Pyl',
  B: 'Asx',
  Z: 'Glx',
  J: 'Xle',
  X: 'Xaa',
  '*': 'Ter',
};
// dense (no holes: a read of one allocates nothing), `Xaa` for a code with none
const ThreeByCode: string[] = Array.from({ length: 128 }, () => 'Xaa');
for (const [k, v] of Object.entries(ThreeLetter)) {
  ThreeByCode[k.charCodeAt(0)] = v;
  ThreeByCode[k.toLowerCase().charCodeAt(0)] = v;
}

/** The data pass: bars, placeholders, gap marks, cells or ticks, the letters on them, carets. */
function drawCells(lay: LettersLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  drawSpans(out, lay.bars, box.y, box.height, tokens.colors['seq-sequence']);
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  // a path's gap marks and carets, and three-letter codes, each in a function of its own: this one stays small enough to
  // have its calls inlined (a call that is not boxes the numbers it is passed); every number through `args`
  if (lay.gaps.length || lay.carets.length) drawPathMarks(lay, out);
  if (lay.isChip && lay.n) {
    drawChip(lay, out);
    return;
  }
  // a tick a position at `dots`, else the cell
  const tick = Math.max(2, Math.round(box.height * 0.3));
  const dots = lay.step === 'dots';
  const y = dots ? box.y + (box.height - tick) / 2 : box.y;
  const h = dots ? tick : box.height;
  const px = 1 / lay.axis.dpr;
  // a letter cell's hairline (`seq-cell-edge`: dark only, none in light), so dim dark cells read apart from the page
  const edge = lay.step === 'letters' ? tokens.colors['seq-cell-edge'] : NoColor;
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    g[0] = lay.x0[j];
    g[1] = y;
    g[2] = Math.max(lay.x1[j] - lay.x0[j], px);
    g[3] = h;
    g[4] = lay.fill[j];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
    if (edge !== NoColor) {
      g[5] = edge;
      g[6] = px;
      cellEdge(out);
    }
  }
  if (lay.conflicts.length) drawConflicts(lay, out);
  if (lay.step !== 'letters' || !(lay.fontSize > 0) || !lay.n) return;
  // a molecule's cells carry no letter: their labels are the molecules' letters' to draw
  if (lay.labels) return;
  if (lay.three) {
    drawThree(lay, out);
    return;
  }
  // a letter centred in its whole cell overhangs the view at an edge: clipped to the band's box
  const clipped = overhangs(lay);
  if (clipped) clipBox(lay, out);
  g[1] = lay.glyphY;
  out.beginGlyphsArgs(lay.font, NoColor);
  for (let j = 0; j < lay.n; j++) {
    if (!lay.show[j]) continue;
    g[0] = lay.gx[j];
    g[4] = lay.ink[j];
    out.letterArgs(lay.code[j]);
  }
  out.endGlyphs();
  if (clipped) out.unclip();
}

/**
 * A letter cell's hairline (`seq-cell-edge`) over the cell `out.args` holds (x, y, w, h; its color in `[5]`, a device
 * pixel in `[6]`: read from the array, a number passed to a call that is not inlined is boxed): drawn inside the cell,
 * a stroke a device pixel wide on the rect inset by half of one, so it covers whole device pixels of its own cell and a
 * neighbour's fill never paints over it (a stroke on the edge itself straddles two pixels). Shared with the rows layer.
 */
export function cellEdge(out: PassWriter) {
  const g = out.args;
  const px = g[6];
  if (g[2] <= 2 * px || g[3] <= 2 * px) return;
  g[0] += px / 2;
  g[1] += px / 2;
  g[2] -= px;
  g[3] -= px;
  g[4] = NoColor;
  out.rectArgs();
}

/** A one-to-one path's gap marks (a thin mark across the middle of an empty cell) and carets (a letter with no position). */
function drawPathMarks(lay: LettersLay, out: PassWriter) {
  const { box, tokens } = lay;
  const px = 1 / lay.axis.dpr;
  const g = out.args;
  for (let k = 0; k < 2; k++) {
    const gaps = k === 0;
    const list = gaps ? lay.gaps : lay.carets;
    const step = gaps ? 3 : 1;
    for (let i = 0; i + step - 1 < list.length; i += step) {
      g[0] = gaps ? list[i] : list[i] - px;
      g[1] = gaps ? box.y + box.height / 2 - 0.5 : box.y;
      g[2] = gaps ? Math.max(list[i + 1] - list[i], 1) : 2 * px;
      g[3] = gaps ? 1 : box.height;
      g[4] = tokens.colors[gaps ? 'seq-missing' : 'msui-text'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
}

/** A conflict's mark: a bar along the bottom of the cell, a fifth of its height (at least 2 px), in the conflict ink. */
function drawConflicts(lay: LettersLay, out: PassWriter) {
  const { box, tokens } = lay;
  const h = Math.max(2, Math.round(box.height * 0.2));
  const g = out.args;
  const c = lay.conflicts;
  for (let i = 0; i + 1 < c.length; i += 2) {
    g[0] = c[i];
    g[1] = box.y + box.height - h;
    g[2] = Math.max(c[i + 1] - c[i], 1);
    g[3] = h;
    g[4] = tokens.colors['seq-conflict'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

/**
 * A chip's pill (_spec/README.md "Render", Blocks): its one cell's fill as a pill a pixel inside the cell, its label
 * centred on it in the ink that holds 4.5:1 against the fill, faded at the pill's width (the label is elided to its
 * slot: it fits); the letter cells' hairline (`seq-cell-edge`, dark only) round it, as round a cell.
 */
function drawChip(lay: LettersLay, out: PassWriter) {
  const { box } = lay;
  const g = out.args;
  const h = Math.max(box.height - 2, 1);
  const edge = lay.tokens.colors['seq-cell-edge'];
  // a chip of several positions (a `whole` molecule): one pill over its cells
  const x0 = lay.x0[0];
  const x1 = lay.x1[lay.n - 1];
  const lo = x0 < x1 ? x0 : lay.x0[lay.n - 1];
  const hi = x0 < x1 ? x1 : lay.x1[0];
  g[0] = lo + 1;
  g[1] = box.y + 1;
  g[2] = hi - lo - 2;
  g[3] = h;
  g[4] = lay.fill[0];
  g[5] = edge;
  g[6] = edge === NoColor ? 0 : 1 / lay.axis.dpr;
  g[7] = h / 2;
  g[8] = 3;
  out.rectArgs();
  // no label fits the slot (one of `labels.minChars` letters or fewer, or not even a cut of so many): the pill alone
  if (!lay.chip) return;
  g[0] = (lo + hi) / 2;
  g[1] = lay.chipY;
  g[2] = lay.chipWidth;
  g[3] = hi - lo;
  g[4] = lay.chipInk;
  out.textArgs(lay.chip, AnchorMiddle, lay.chipFont);
}

/** A chip's label ink against its fill: the tokens' ink when it holds 4.5:1, else black or white, the higher. */
export function chipInk(fill: number, tokens: SeqTokens): number {
  const ink = inkOn(fill, tokens);
  if (contrastRgba(ink, fill) >= 4.5) return ink;
  return contrastRgba(Black, fill) >= contrastRgba(White, fill) ? Black : White;
}

/** Whether a shown cell's whole extent reaches past the data box (its letter is cut by the band's box there). */
export function overhangs(lay: {
  n: number;
  show: Uint8Array;
  gx: Float32Array;
  gw: Float32Array;
  axis: AxisTransform;
}): boolean {
  // (a cell snapped to the device pixel grid may sit a pixel past an unsnapped edge: that is not an overhang)
  const px = 1 / lay.axis.dpr;
  const left = lay.axis.x0 - px;
  const right = lay.axis.x0 + lay.axis.width + px;
  for (let j = 0; j < lay.n; j++) {
    if (!lay.show[j]) continue;
    const h = lay.gw[j] / 2;
    if (lay.gx[j] - h < left || lay.gx[j] + h > right) return true;
  }
  return false;
}

/** A clip to the data box of the band's row (numbers by `args`: a fractional number passed to a call is boxed). */
export function clipBox(
  lay: { box: { readonly y: number; readonly height: number }; axis: AxisTransform },
  out: PassWriter,
) {
  const g = out.args;
  g[0] = lay.axis.x0;
  g[1] = lay.box.y;
  g[2] = lay.axis.width;
  g[3] = lay.box.height;
  out.clipArgs();
}

/** Three-letter codes, centred in their cells. */
function drawThree(lay: LettersLay, out: PassWriter) {
  const clipped = overhangs(lay);
  if (clipped) clipBox(lay, out);
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    if (!lay.show[j]) continue;
    g[0] = lay.gx[j];
    g[1] = lay.glyphY;
    g[2] = 0;
    g[3] = lay.gw[j];
    g[4] = lay.ink[j];
    out.textArgs(ThreeByCode[lay.code[j] & 127], AnchorMiddle, lay.font);
  }
  if (clipped) out.unclip();
}

/**
 * The hover pass: a letter through a path outlines each part of its element (a residue's codon, split by an intron
 * too); a chip, its pill; else the hovered position's cell. Each a squircle (_spec/README.md "Render", Shapes: every
 * emphasis a squircle), concentric: a part's ends round where the element really ends (its outermost parts' outer
 * sides), square where an intron cuts it; square below `squircleMin`; a chip's round as its pill.
 */
function drawHover(lay: LettersLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key) return;
  if (lay.isChip && lay.n > 0) {
    drawChipHover(lay, out, emphasis);
    return;
  }
  const part = lay.through === 'own' ? undefined : h.extent?.parts[0];
  if (!part || part.system !== lay.dataSystem || !part.ranges.length) {
    hoverCell(lay, out, emphasis);
    return;
  }
  const p = part.ranges[0];
  const d = lay.hoverWidth / 2 + lay.hoverOffset;
  const g = out.args;
  const { box } = lay;
  // the element's outer sides: its parts' leftmost and rightmost x (the sides between parts are cut)
  let lo = Infinity;
  let hi = -Infinity;
  for (let j = 0; j < lay.n; j++) {
    if (lay.pos[j] !== p) continue;
    if (lay.x0[j] < lo) lo = lay.x0[j];
    if (lay.x1[j] > hi) hi = lay.x1[j];
  }
  const c = lay.config;
  const bh = box.height;
  for (let j = 0; j < lay.n; j++) {
    if (lay.pos[j] !== p) continue;
    const w = lay.x1[j] - lay.x0[j];
    const r =
      w < c.squircleMin || bh < c.squircleMin ? 0 : Math.min(bh * c.squircleRatio, c.squircleMax, w / 2, bh / 2);
    const round = (lay.x0[j] <= lo ? 1 : 0) | (lay.x1[j] >= hi ? 2 : 0);
    g[0] = lay.x0[j] - d;
    g[1] = box.y - d;
    g[2] = w + 2 * d;
    g[3] = bh + 2 * d;
    g[4] = NoColor;
    g[5] = lay.tokens.colors['seq-hover'];
    g[6] = lay.hoverWidth;
    g[7] = r > 0 && round !== 0 ? r + d : 0;
    g[8] = round;
    out.rectArgs();
  }
}

/** A chip's hover: outside its pill by the offset, concentric with it (the pill's radius plus the gap to it). */
function drawChipHover(lay: LettersLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis!.hover!;
  if (h.system !== lay.axis.system || h.position === undefined || !lay.axis.cell(h.position)) return;
  const { box } = lay;
  // the pill is a pixel inside the cells (`drawChip`); the outline `d` outside them
  const d = lay.hoverWidth / 2 + lay.hoverOffset;
  const ph = Math.max(box.height - 2, 1);
  const g = out.args;
  const x0 = lay.x0[0];
  const x1 = lay.x1[lay.n - 1];
  const lo = x0 < x1 ? x0 : lay.x0[lay.n - 1];
  const hi = x0 < x1 ? x1 : lay.x1[0];
  g[0] = lo - d;
  g[1] = box.y - d;
  g[2] = hi - lo + 2 * d;
  g[3] = box.height + 2 * d;
  g[4] = NoColor;
  g[5] = lay.tokens.colors['seq-hover'];
  g[6] = lay.hoverWidth;
  g[7] = ph / 2 + 1 + d;
  g[8] = 3;
  out.rectArgs();
}

export const LettersDraws = passDraws<LettersLay, Emphasis | undefined>({ data: drawCells, hover: drawHover });
