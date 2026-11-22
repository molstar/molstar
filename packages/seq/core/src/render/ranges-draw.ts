/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Anchor, grown, NoColor, passDraws, type PassWriter } from './draw-lists.js';
import { drawTicks } from './items.js';
import type { Emphasis } from './layer-type.js';
import { drawConnectors, hoverConnectors, selectConnectors } from './ranges-connector.js';
import { arrowArgs, drawFailed, drawSpans } from './marks.js';
import { labelText, labelX } from './range-labels.js';
import { Gene, type RangesLay, Solid } from './ranges-lay.js';

// A ranges layer's draw (_spec/README.md "Render", Layer types; Shapes): its marks as squircles, lines or arrows, its
// labels inside and beside with their leaders, the hovered item's outline and a selected boundary's hairline.

// the anchor a number read once: a draw reads no property of an imported object (that site's feedback, cleared, made
// a warm draw deoptimize and allocate without end once the 150,000-piece scenario had run first)
const AnchorStart: number = Anchor.start;

/** The dash of a dashed and a dotted line (on, then as long off), the wave's half period and amplitude, CSS px. */
const Dash = 6;
const Dot = 2;
export const WaveStep = 3;
const WaveAmplitude = 2;
const LineWidth = 2;

/** A gene line's thin box (a UTR's): this share of a thick one's height, about its lane's middle. */
const ThinShare = 0.5;

/** A highlight's alpha: 0.35 of 255, the channel's default. */
const HighlightAlpha = 89;

/** A wave's scratch points (x, y pairs), kept and grown at layout. */
let wave = new Float64Array(64);
/** Grows the wave's scratch to hold `n` numbers (at layout: a draw makes none). */
export function growWave(n: number): void {
  wave = grown(wave, n);
}
const segment = new Float64Array(4);

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: placeholders, the density strip, the marks. */
function drawMarks(lay: RangesLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  const gap = lay.gap;
  // in a function of its own: this one stays small enough to have its calls inlined (a call that is not boxes the
  // numbers it is passed)
  if (lay.tn) drawTicks(lay, out);
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  const a = out.args;
  if (lay.density) {
    a[0] = box.x;
    a[1] = lay.densityY;
    a[2] = box.width;
    a[3] = lay.densityH;
    out.imageArgs(lay.density);
  }
  if (lay.mode === 'density') return;
  if (lay.cn) drawConnectors(lay, out);
  if (lay.style !== 0 || lay.arrow) {
    drawStyled(lay, out);
    return;
  }
  // items and merged runs alike are discrete marks: squircles where big enough, round only where a run really ends,
  // so the base mark has the shape its emphasis outlines (_spec/README.md "Render", Shapes)
  const full = lay.laneHeight - gap;
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    // the size channel shrinks the mark about its lane's middle; its lane and hit box stay (EC-137)
    const h = item >= 0 ? full * lay.sizes[item] : full;
    a[0] = lay.mx0[j];
    a[1] = lay.laneY + lay.mlane[j] * lay.laneHeight + gap / 2 + (full - h) / 2;
    a[2] = lay.mx1[j] - lay.mx0[j];
    a[3] = h;
    a[4] = item >= 0 ? lay.fills[item] : item <= -2 ? lay.fills[-2 - item] : tokens.colors['seq-range'];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = lay.mradius[j];
    a[8] = lay.mround[j];
    out.rectArgs();
  }
}

/**
 * The marks of a layer with a line style or arrowheads: a block with a pointed end where it holds an end the item really
 * has, or the line (its end shortened to where its triangle begins). The dash and the wave run from the content's start
 * (`manchor`), so a pan moves them with the data. No closure, nothing made per mark: the numbers go through `args`.
 */
function drawStyled(lay: RangesLay, out: PassWriter) {
  const { tokens } = lay;
  const gap = lay.gap;
  const a = out.args;
  const style = lay.style;
  const forward = lay.axis.dir === 1;
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    const fill = item >= 0 ? lay.fills[item] : item <= -2 ? lay.fills[-2 - item] : tokens.colors['seq-range'];
    const full = lay.laneHeight - gap;
    const h = item >= 0 ? full * lay.sizes[item] : full;
    const y = lay.laneY + lay.mlane[j] * lay.laneHeight + gap / 2 + (full - h) / 2;
    const cy = y + h / 2;
    const x0 = lay.mx0[j];
    const x1 = lay.mx1[j];
    const tips = lay.mtip[j];
    const head = Math.min(h * 0.9, x1 - x0);
    if (style === Gene) {
      // a box: full between the thick part's ends, thin outside (no arrowheads: the hairline's chevrons say the strand)
      const t = lay.mthin[j] ? h * ThinShare : h;
      a[0] = x0;
      a[1] = y + (h - t) / 2;
      a[2] = x1 - x0;
      a[3] = t;
      a[4] = fill;
      a[5] = NoColor;
      a[6] = 0;
      a[7] = 0;
      a[8] = 0;
      out.rectArgs();
      continue;
    }
    if (style === 0 || style === 7) {
      // a block (a highlight is a translucent one); a tip makes it an arrow, both tips two overlapping ones
      const color = style === 0 ? fill : ((fill & 0xffffff00) | HighlightAlpha) >>> 0;
      a[0] = x0;
      a[1] = y;
      a[2] = x1 - x0;
      a[3] = h;
      a[4] = color;
      a[5] = NoColor;
      a[6] = 0;
      if (!tips) {
        a[7] = lay.mradius[j];
        a[8] = lay.mround[j];
        out.rectArgs();
        continue;
      }
      if (tips & 2) {
        a[7] = head;
        a[8] = 1;
        arrowArgs(out);
        a[0] = x0;
        a[2] = x1 - x0;
        a[4] = color;
        a[5] = NoColor;
        a[6] = 0;
      }
      if (tips & 1) {
        a[7] = head;
        a[8] = -1;
        arrowArgs(out);
      }
      continue;
    }
    // a line: from where the content starts, over what is in view, ending where a triangle begins
    const l0 = tips & 1 ? x0 + head : x0;
    const l1 = tips & 2 ? x1 - head : x1;
    // the start's triangle, on the side the content starts from: the dash and the wave begin at the line, not under it
    const startTip = forward ? tips & 1 : tips & 2;
    if (style === Solid) out.fill(l0, cy - LineWidth / 2, l1 - l0, LineWidth, fill);
    else if (style === 2 || style === 3) {
      const dash = style === 2 ? Dash : Dot;
      const period = 2 * dash;
      const anchor = lay.manchor[j];
      const m = Math.max(0, (startTip ? Math.ceil : Math.floor)((forward ? l0 - anchor : anchor - l1) / period));
      segment[0] = forward ? anchor + m * period : anchor - m * period;
      segment[1] = cy;
      segment[2] = forward ? l1 : l0;
      segment[3] = cy;
      a[4] = NoColor;
      a[5] = fill;
      a[6] = LineWidth;
      a[7] = dash;
      out.pathArgs(segment, 2, false);
    } else if (style === 4) {
      // the wave from the content's start, a corner every half period, phase kept when it starts before the view
      const period = 2 * WaveStep;
      const anchor = lay.manchor[j];
      const m = Math.max(0, (startTip ? Math.ceil : Math.floor)((forward ? l0 - anchor : anchor - l1) / period));
      const xs = forward ? anchor + m * period : anchor - m * period;
      const xe = forward ? l1 : l0;
      const cap = wave.length / 2;
      const count = Math.min(cap, Math.floor(Math.abs(xe - xs) / WaveStep) + 2);
      for (let i = 0; i < count; i++) {
        const x = forward ? xs + i * WaveStep : xs - i * WaveStep;
        wave[2 * i] = forward ? Math.min(x, xe) : Math.max(x, xe);
        wave[2 * i + 1] = cy + (i & 1 ? WaveAmplitude : -WaveAmplitude);
      }
      a[4] = NoColor;
      a[5] = fill;
      a[6] = 1;
      a[7] = 0;
      out.pathArgs(wave, count, false);
    } else if (style === 5) {
      out.fill(l0, cy - LineWidth, l1 - l0, 1, fill);
      out.fill(l0, cy + LineWidth - 1, l1 - l0, 1, fill);
    } else {
      // a box: the lane's outline, inside the mark
      a[0] = x0 + 0.5;
      a[1] = y + 0.5;
      a[2] = x1 - x0 - 1;
      a[3] = h - 1;
      a[4] = NoColor;
      a[5] = fill;
      a[6] = 1;
      a[7] = 0;
      a[8] = 0;
      out.rectArgs();
    }
    if (tips) {
      // the triangle: a head with no body, over the line's end
      a[1] = y;
      a[3] = h;
      a[4] = fill;
      a[5] = NoColor;
      a[6] = 0;
      if (tips & 2) {
        a[0] = x1 - head;
        a[2] = head;
        a[7] = head;
        a[8] = 1;
        arrowArgs(out);
      }
      if (tips & 1) {
        a[0] = x0;
        a[2] = head;
        a[4] = fill;
        a[5] = NoColor;
        a[6] = 0;
        a[7] = head;
        a[8] = -1;
        arrowArgs(out);
      }
    }
  }
}

/** The width of a scrolled band's rail (a slim line at its right edge), CSS px. */
const RailPx = 3;

/**
 * A band whose lanes scroll inside it (_spec/README.md "Render", Track heights): a 1 px edge in the line ink where lanes
 * are cut (above the first lane in view, below the last), and a slim rail at its right edge whose window is the lanes in
 * view; none while every lane fits.
 */
function drawScroll(lay: RangesLay, out: PassWriter) {
  if (!lay.scrolls || !(lay.scrollMax > 0)) return;
  const { box, tokens } = lay;
  const g = out.args;
  const edge = tokens.colors['msui-line'];
  // unrolled: a number passed to a call that is not inlined is boxed (a warm frame allocates nothing)
  g[4] = edge;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  if (lay.top > 0) {
    g[0] = box.x;
    g[1] = box.y;
    g[2] = box.width;
    g[3] = 1;
    out.rectArgs();
  }
  if (lay.top < lay.scrollMax) {
    g[0] = box.x;
    g[1] = box.y + box.height - 1;
    g[2] = box.width;
    g[3] = 1;
    out.rectArgs();
  }
  const room = box.height - 4;
  const h = Math.max(8, (room * lay.room) / (lay.scrollMax + lay.room));
  g[0] = box.x + box.width - RailPx - 1;
  g[1] = box.y + 2;
  g[2] = RailPx;
  g[3] = room;
  g[4] = edge;
  g[6] = RailPx / 2;
  out.rectArgs();
  g[1] = box.y + 2 + ((room - h) * lay.top) / lay.scrollMax;
  g[3] = h;
  g[4] = tokens.colors['seq-window'];
  out.rectArgs();
}

/** The labels pass: the labels inside their marks (whole or cut), and those beside in the label lane with their leaders. */
function drawLabels(lay: RangesLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const gap = lay.gap;
  const g = out.args;
  drawScroll(lay, out);
  if (lay.style === Gene) {
    // a gene model's label: plain text in the text ink, left of its model (placed whole, `ranges-gene.ts`)
    const gl = lay.gene;
    const ink = lay.tokens.colors['msui-text'];
    const bg = lay.tokens.colors['msui-bg'];
    for (let k = 0; k < gl.ln; k++) {
      if (gl.lfloat[k]) {
        // floating over its model's hairline: the page's color under the text
        g[0] = gl.lx[k] - 2;
        g[1] = gl.ly[k] - lay.fontSize * 0.85;
        g[2] = gl.lw[k] + 4;
        g[3] = lay.fontSize;
        g[4] = bg;
        g[5] = NoColor;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        out.rectArgs();
      }
      g[0] = gl.lx[k];
      g[1] = gl.ly[k];
      g[2] = gl.lw[k];
      g[3] = gl.lw[k];
      g[4] = ink;
      out.textArgs(lay.labels[gl.litem[k]] as string, AnchorStart, lay.font);
    }
    return;
  }
  for (let j = 0; j < lay.m; j++) {
    if (!(lay.mlabel[j] > 0)) continue;
    const y = lay.laneY + lay.mlane[j] * lay.laneHeight + gap / 2;
    const h = lay.laneHeight - gap;
    const baseline = y + h / 2 + lay.fontSize * 0.35;
    if (lay.style >= 1 && lay.style <= 5) {
      // a line runs through its label's place: the page's color under the text, so the label reads (the line stops at it)
      g[0] = labelX(lay, j) - 2;
      g[1] = baseline - lay.fontSize * 0.85;
      g[2] = lay.mlabel[j] + 4;
      g[3] = lay.fontSize;
      g[4] = lay.tokens.colors['msui-bg'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
    g[0] = labelX(lay, j);
    g[1] = baseline;
    g[2] = lay.mlabel[j];
    g[3] = lay.mlabel[j];
    g[4] = lay.mink[j];
    out.textArgs(labelText(lay, j), AnchorStart, lay.font);
  }
  if (!lay.ln) return;
  const mid = lay.labelY + lay.laneHeight / 2;
  const baseline = mid + lay.fontSize * 0.35;
  const ink = lay.tokens.colors['msui-text'];
  const muted = lay.tokens.colors['msui-text-muted'];
  for (let k = 0; k < lay.ln; k++) {
    const x = lay.lx[k];
    const want = lay.lwant[k];
    if (x - want >= 1) {
      // moved: a short leader from its mark's start (its own lane's top) up to the label's line and on to the label
      const xy = Leader;
      xy[0] = want + 0.5;
      xy[1] = lay.laneY + Math.max(0, lay.lanes[lay.litem[k]]) * lay.laneHeight + gap / 2;
      xy[2] = want + 0.5;
      xy[3] = mid;
      xy[4] = x - 2;
      xy[5] = mid;
      g[4] = NoColor;
      g[5] = muted;
      g[6] = 1;
      g[7] = 0;
      out.pathArgs(xy, 3, false);
    }
    g[0] = x;
    g[1] = baseline;
    g[2] = lay.lw[k];
    g[3] = lay.lw[k];
    g[4] = ink;
    out.textArgs(lay.labels[lay.litem[k]] as string, AnchorStart, lay.font);
  }
}

/** A leader's three points, written per label (kept). */
const Leader = new Float64Array(6);

/**
 * The hover pass: an outline around each of the hovered item's own marks (its hit marks: inside a merged run too),
 * concentric.
 */
function drawHover(lay: RangesLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const { tokens } = lay;
  const gap = lay.gap;
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || !h.element) return;
  // outside the mark by the offset, concentric: it never touches an outline drawn on the mark
  const half = lay.hoverWidth / 2 + lay.hoverOffset;
  const a = out.args;
  if (lay.cn) hoverConnectors(lay, out, h.element.key);
  for (let j = 0; j < lay.hn; j++) {
    const item = lay.hitem[j];
    if (item < 0 || lay.keys[item] !== h.element.key) continue;
    const r = lay.hradius[j];
    // concentric: the outline's radius is the mark's plus its offset
    a[0] = lay.hx0[j] - half;
    a[1] = lay.laneY + lay.hlane[j] * lay.laneHeight + gap / 2 - half;
    a[2] = lay.hx1[j] - lay.hx0[j] + 2 * half;
    a[3] = lay.laneHeight - gap + 2 * half;
    a[4] = NoColor;
    a[5] = tokens.colors['seq-hover'];
    a[6] = lay.hoverWidth;
    a[7] = r > 0 ? r + half : 0;
    a[8] = lay.hround[j];
    out.rectArgs();
  }
}

/**
 * The selection pass: a selected item of zero length (a boundary, an insertion point) is a hairline at its place over
 * its lane in the accent (EC-127), which the selection's Location cannot hold; the rest is the view's.
 */
function drawSelected(lay: RangesLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const els = emphasis?.elements;
  if (!els || !els.keys.length) return;
  const { tokens } = lay;
  const dpr = lay.axis.dpr;
  const a = out.args;
  if (lay.cn) selectConnectors(lay, out, els.keys);
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    if (item < 0 || lay.starts[item] !== lay.ends[item]) continue;
    let selected = false;
    for (let k = 0; k < els.keys.length && !selected; k++)
      selected = els.keys[k].owner === lay.owner && els.keys[k].key === lay.keys[item];
    if (!selected) continue;
    const x = Math.round(((lay.mx0[j] + lay.mx1[j]) / 2) * dpr) / dpr;
    a[0] = x - 1;
    a[1] = lay.laneY + lay.mlane[j] * lay.laneHeight;
    a[2] = 2;
    a[3] = lay.laneHeight;
    a[4] = tokens.colors['msui-accent'];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
  }
}

export const RangesDraws = passDraws<RangesLay, Emphasis | undefined>({
  data: drawMarks,
  labels: drawLabels,
  selection: drawSelected,
  hover: drawHover,
});
