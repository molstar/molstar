/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { Anchor, NoColor, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { LayoutContext } from './layer-type.js';
import { resolveColor } from './tokens.js';
import { elide } from './text.js';
import { domainHi, domainLo, Fill, ys } from './values-marks.js';
import type { ValuesLay } from './values.js';

// The values layer's reference lines (_spec/README.md "Render", Layer types, `values`): up to four hairlines across the
// well, each at a value of the layer's y scale (a threshold, a zero, a mean), through the same `ys` and headroom as a
// value, so a dot at the threshold and its line coincide. Drawn after the data and under the hover, never hit. A line
// outside a closed domain is not drawn (said once); an open end of the domain includes the lines. A label draws at the
// well's right end in the `labels` pass (a patch of the well's ground under the text): above the line, below it near
// the well's top, only from `labelMinWellPx` of well, the later of two overlapping ones dropped. Resolved once a layout
// into preallocated arrays: a draw reads numbers only (nothing allocated once warm).

/** The most lines a layer takes (the description's bound). */
export const MaxReferenceLines = 4;

const Anchor_Start: number = Anchor.start;

/** A layout holding its reference lines' state (the chunk landed). */
export type Held = ValuesLay & { ref: ReferenceState };

/** A line as the spec writes it: its style is 0 solid, 1 dashed, 2 dotted. */
export interface ReferenceLine {
  readonly at: number;
  readonly label: string;
  readonly color: string | undefined;
  readonly style: 0 | 1 | 2;
}

const NoLines: readonly ReferenceLine[] = [];
const Styles = ['solid', 'dashed', 'dotted'] as const;

/** The lines a layer's `reference` param holds (a record without a number `at` is left out; at most four). */
export function readReference(value: unknown): readonly ReferenceLine[] {
  if (!Array.isArray(value) || !value.length) return NoLines;
  const out: ReferenceLine[] = [];
  for (const r of value) {
    if (out.length >= MaxReferenceLines || !r || typeof r !== 'object') continue;
    const { at, label, color, style } = r as Record<string, unknown>;
    if (typeof at !== 'number' || !Number.isFinite(at)) continue;
    const s = Styles.indexOf(style as (typeof Styles)[number]);
    out.push({
      at,
      label: typeof label === 'string' ? label : '',
      color: typeof color === 'string' ? color : undefined,
      style: (s < 0 ? 1 : s) as 0 | 1 | 2,
    });
  }
  return out.length ? out : NoLines;
}

/** A layout's resolved lines, kept on the layer's layout (its arrays reused). */
export class ReferenceState {
  lines: readonly ReferenceLine[] = NoLines;
  readonly shown = new Uint8Array(MaxReferenceLines);
  readonly y = new Float64Array(MaxReferenceLines);
  /** Each line's top edge, snapped to a device pixel, and the hairline's thickness (the draw reads numbers only). */
  readonly top = new Float64Array(MaxReferenceLines);
  width = 1;
  readonly color = new Uint32Array(MaxReferenceLines);
  /** The dash of a line (0: solid). */
  readonly dash = new Float64Array(MaxReferenceLines);
  /** Each label: drawn, as drawn (cut), its patch (x, y, w, h) and its baseline. */
  readonly labelShown = new Uint8Array(MaxReferenceLines);
  readonly text: string[] = ['', '', '', ''];
  readonly patch = new Float64Array(4 * MaxReferenceLines);
  readonly baseline = new Float64Array(MaxReferenceLines);
  font = '10px sans';
  ink = 0;
  ground = 0;
}

const extent = new Float64Array(3);

/**
 * What the lines add to an open domain, the seeds of its extent: the lowest and the highest line, and the lowest
 * positive (a log scale's low end; its lines at or under zero are not on it).
 */
export function referenceExtent(lay: Pick<Held, 'ref' | 'yLog'>): Float64Array {
  extent[0] = Number.POSITIVE_INFINITY;
  extent[1] = Number.NEGATIVE_INFINITY;
  extent[2] = Number.POSITIVE_INFINITY;
  for (const l of lay.ref.lines) {
    if (lay.yLog && l.at <= 0) continue;
    if (l.at < extent[0]) extent[0] = l.at;
    if (l.at > extent[1]) extent[1] = l.at;
    if (l.at > 0 && l.at < extent[2]) extent[2] = l.at;
  }
  return extent;
}

/** The ends the writer closed: bit 0 the low, bit 1 the high (the layer's `domain`, or the scale node's own). */
export function closedEnds(
  domain: readonly [number, number] | undefined,
  scale: readonly (number | null | undefined)[] | undefined,
): number {
  if (domain) return 3;
  return (typeof scale?.[0] === 'number' ? 1 : 0) | (typeof scale?.[1] === 'number' ? 2 : 0);
}

const Pad = 3;

/**
 * Resolves the lines against the final domain: each one's y, color and dash, which are not drawn (outside the domain,
 * or at or under zero on a log scale), then the labels' places. Says once the first line outside a closed end.
 */
export function layoutReference(lay: Held, ctx: LayoutContext, closed: number): void {
  const r = lay.ref;
  const lines = r.lines;
  r.shown.fill(0);
  r.labelShown.fill(0);
  if (!lines.length) return;
  const cfg = ctx.config.values.reference;
  const { box, tokens } = lay;
  const lo = ys(lay, domainLo(lay));
  const hi = ys(lay, domainHi(lay));
  const span = (hi - lo || 1) / Fill;
  const eps = 1e-9 * (hi - lo || 1);
  const y1 = box.y + box.height;
  const muted = tokens.colors['msui-text-muted'];
  r.width = Math.max(cfg.widthPx, 1 / lay.axis.dpr);
  let said = false;
  for (let i = 0; i < lines.length; i++) {
    const l = lines[i];
    const v = ys(lay, l.at);
    // a log scale's line at or under zero is not on the scale (the layer's `seq.render.scale.log` says why)
    if (!(v > Number.NEGATIVE_INFINITY)) continue;
    const below = v < lo - eps;
    if (below || v > hi + eps) {
      if (!said && (below ? closed & 1 : closed & 2)) {
        said = true;
        ctx.report(
          issueOf(
            RenderIssues,
            'seq.render.reference.outside',
            `The reference line ${i} at ${l.at} is outside the domain [${domainLo(lay)}, ${domainHi(lay)}]`,
            { address: ctx.key, hint: RenderIssues.codes['seq.render.reference.outside'].hint },
          ),
        );
      }
      continue;
    }
    r.shown[i] = 1;
    r.y[i] = y1 - ((v - lo) / span) * box.height;
    r.top[i] = Math.round((r.y[i] - r.width / 2) * lay.axis.dpr) / lay.axis.dpr;
    r.color[i] = (l.color !== undefined ? resolveColor(l.color, tokens) : undefined) ?? tokens.colors['msui-text'];
    r.dash[i] = l.style === 0 ? 0 : l.style === 1 ? cfg.dashPx : cfg.dotPx;
  }
  layoutLabels(lay, ctx, muted);
}

function layoutLabels(lay: Held, ctx: LayoutContext, muted: number) {
  const r = lay.ref;
  const cfg = ctx.config.values.reference;
  const { box } = lay;
  if (box.height < cfg.labelMinWellPx) return;
  const h = cfg.labelPx + 2;
  r.font = `${cfg.labelPx}px sans`;
  r.ink = muted;
  r.ground = lay.tokens.colors['seq-track'];
  const maxText = box.width / 3 - 2 * Pad;
  for (let i = 0; i < r.lines.length; i++) {
    const label = r.lines[i].label;
    if (!r.shown[i] || !label || maxText <= 0) continue;
    const text = elide(ctx.text, r.font, label, maxText);
    if (!text) continue;
    const w = ctx.text.measure(r.font, text) + 2 * Pad;
    const y = r.y[i];
    // above the line, below it when the well's top is within a label height
    // clear of the hairline's whole thickness: above its top edge, below its bottom edge
    const py = y - box.y >= h ? r.top[i] - h : r.top[i] + r.width;
    if (py + h > box.y + box.height) continue;
    const px = box.x + box.width - w;
    let clash = false;
    for (let k = 0; k < i && !clash; k++) {
      if (!r.labelShown[k]) continue;
      clash =
        px < r.patch[4 * k] + r.patch[4 * k + 2] &&
        r.patch[4 * k] < px + w &&
        py < r.patch[4 * k + 1] + h &&
        r.patch[4 * k + 1] < py + h;
    }
    if (clash) continue;
    r.labelShown[i] = 1;
    r.text[i] = text;
    r.patch[4 * i] = px;
    r.patch[4 * i + 1] = py;
    r.patch[4 * i + 2] = w;
    r.patch[4 * i + 3] = h;
    r.baseline[i] = py + (h + 0.7 * cfg.labelPx) / 2;
  }
}

const pts = new Float64Array(4);

/** The lines, after the data: a solid one a rect, a dashed or dotted one a dashed path. */
export function drawReference(lay: Held, out: PassWriter): void {
  const r = lay.ref;
  const n = r.lines.length;
  if (!n) return;
  const { box } = lay;
  const a = out.args;
  for (let i = 0; i < n; i++) {
    if (!r.shown[i]) continue;
    if (r.dash[i] === 0) {
      a[0] = box.x;
      a[1] = r.top[i];
      a[2] = box.width;
      a[3] = r.width;
      a[4] = r.color[i];
      a[5] = NoColor;
      a[6] = 0;
      a[7] = 0;
      a[8] = 0;
      out.rectArgs();
    } else {
      pts[0] = box.x;
      pts[1] = r.top[i] + r.width / 2;
      pts[2] = box.x + box.width;
      pts[3] = pts[1];
      a[4] = NoColor;
      a[5] = r.color[i];
      a[6] = r.width;
      a[7] = r.dash[i];
      out.pathArgs(pts, 2, false);
    }
  }
}

/** The labels: a patch of the well's ground, then the text over it. */
export function drawReferenceLabels(lay: Held, out: PassWriter): void {
  const r = lay.ref;
  const n = r.lines.length;
  const a = out.args;
  for (let i = 0; i < n; i++) {
    if (!r.labelShown[i]) continue;
    a[0] = r.patch[4 * i];
    a[1] = r.patch[4 * i + 1];
    a[2] = r.patch[4 * i + 2];
    a[3] = r.patch[4 * i + 3];
    a[4] = r.ground;
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
    a[0] = r.patch[4 * i] + Pad;
    a[1] = r.baseline[i];
    a[2] = r.patch[4 * i + 2] - 2 * Pad;
    a[3] = Number.POSITIVE_INFINITY;
    a[4] = r.ink;
    out.textArgs(r.text[i], Anchor_Start, r.font);
  }
}

/** The layer's accessible words: "reference line at 7.3 (genome-wide)" a drawn line, in order. */
export function referenceWords(lay: Pick<Held, 'ref'>): string | undefined {
  const r = lay.ref;
  const words: string[] = [];
  for (let i = 0; i < r.lines.length; i++) {
    if (!r.shown[i]) continue;
    const l = r.lines[i];
    words.push(`reference line at ${Number(l.at.toPrecision(6))}${l.label ? ` (${l.label})` : ''}`);
  }
  return words.length ? words.join(', ') : undefined;
}
