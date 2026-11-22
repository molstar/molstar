/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid's categorical axes' chrome (_spec/README.md "Render", The grid layout: categorical axes; EC-35): no ruler and
// no first and last positions; the members' labels, thinned past `rows.labelMinPx` a label, down the label column for
// y and along the top at 45° for x, the plot's box grown by the x labels' room (fixed: never moved by a zoom). A click
// on a label is the member (`target: 'row'`, its rows system and position). A draw allocates nothing: a label's letters'
// advances are measured at place and kept by font and text.

import type { AxisTransform } from '../axis.js';
import type { RenderConfig } from '../config.js';
import { Anchor, type PassWriter } from '../draw-lists.js';
import type { Hit } from '../layer-type.js';
import type { TextMetrics } from '../text.js';
import type { SectionLayout } from '../tracks-types.js';

const AnchorEnd: number = Anchor.end;

/** The x labels' slant: up and to the right at 45°. */
const Slant = Math.SQRT1_2;
const Angle = -Math.PI / 4;

/** The x labels' gap above the plot, px. */
const Gap = 4;

/** A label's letters' advances (CSS px), by font and text, kept by the metrics (made at place, read by a draw). */
const Advances = new WeakMap<TextMetrics, Map<string, Float32Array>>();

function advancesOf(text: TextMetrics, font: string, label: string): Float32Array {
  let byKey = Advances.get(text);
  if (!byKey) Advances.set(text, (byKey = new Map()));
  const key = `${font}|${label}`;
  let a = byKey.get(key);
  if (!a) {
    a = new Float32Array(label.length);
    for (let i = 0; i < label.length; i++) a[i] = text.measure(font, label[i]);
    byKey.set(key, a);
  }
  return a;
}

/** The members' labels a categorical axis shows, and where (made at place). */
export class MemberLabels {
  /** x's labels' letters' advances by member (empty: x is no categorical axis). */
  advances: Float32Array[] = [];
  /** The x labels' room above the plot (0: none), and their longest run along the slant, px. */
  room = 0;
  run = 0;
  /** What the advances were measured for. */
  labels: readonly string[] | undefined;
  font = '';
  text: TextMetrics | undefined;
  textVersion = -1;
  cap = -1;
  /** The hit a label gives (one, reused). */
  readonly hit: { -readonly [K in keyof Hit]: Hit[K] } = { layer: '', section: '', system: '' };
}

/**
 * The x labels' room above the plot when x is categorical: the longest label's run (capped at the label column's width)
 * at 45°, its gap; their advances measured (each label once a font).
 */
export function placeMemberLabels(
  m: MemberLabels,
  axis: AxisTransform | undefined,
  text: TextMetrics,
  config: RenderConfig,
): number {
  const labels = axis?.members?.labels;
  // measured once a set of labels, font, metrics and cap: a warm place measures nothing
  if (
    labels === m.labels &&
    config.itemFont === m.font &&
    text === m.text &&
    text.version === m.textVersion &&
    config.labelWidth === m.cap
  )
    return m.room;
  m.labels = labels;
  m.font = config.itemFont;
  m.text = text;
  m.textVersion = text.version;
  m.cap = config.labelWidth;
  m.advances.length = 0;
  m.room = 0;
  m.run = 0;
  if (!labels) return 0;
  let run = 0;
  for (let i = 0; i < labels.length; i++) {
    const a = advancesOf(text, config.itemFont, labels[i]);
    m.advances.push(a);
    let w = 0;
    for (let k = 0; k < a.length; k++) w += a[k];
    if (w > run) run = w;
  }
  m.run = Math.min(run, config.labelWidth);
  m.room = Math.ceil(m.run * Slant + Gap + 2);
  return m.room;
}

/** Every `step`-th display offset is labelled: a label needs `rows.labelMinPx`. */
function stepOf(axis: AxisTransform, config: RenderConfig): number {
  return Math.max(1, Math.ceil(config.rows.labelMinPx / Math.max(axis.scale, 1e-9)));
}

/** The members' labels: y's down the label column at its right edge, x's at 45° above the plot. */
export function drawMemberLabels(
  s: SectionLayout,
  m: MemberLabels,
  ax: AxisTransform,
  ay: AxisTransform | undefined,
  plotTop: number,
  out: PassWriter,
) {
  const { config, tokens } = s;
  const g = out.args;
  const muted = tokens.colors['msui-text-muted'];
  const ym = ay?.members;
  if (ay && ym) {
    const step = stepOf(ay, config);
    const right = ax.x0 - config.labelGap / 2;
    const first = Math.max(0, Math.floor(ay.u0));
    const end = Math.min(ym.order.length, Math.ceil(ay.u1));
    for (let d = first - (first % step); d < end; d += step) {
      if (d < first) continue;
      const y = (ay.xOfUnit(d) + ay.xOfUnit(d + 1)) / 2;
      g[0] = right;
      g[1] = y + 4;
      g[2] = 0;
      g[3] = config.labelWidth;
      g[4] = muted;
      out.textArgs(ym.labels[ym.order[d]], AnchorEnd, config.itemFont);
    }
  }
  const xm = ax.members;
  if (!xm || !m.advances.length) return;
  const step = stepOf(ax, config);
  const first = Math.max(0, Math.floor(ax.u0));
  const end = Math.min(xm.order.length, Math.ceil(ax.u1));
  const base = plotTop - Gap;
  g[1] = base;
  out.beginGlyphsArgs(config.itemFont, muted);
  for (let d = first - (first % step); d < end; d += step) {
    if (d < first) continue;
    const cx = (ax.xOfUnit(d) + ax.xOfUnit(d + 1)) / 2;
    const label = xm.labels[xm.order[d]];
    const adv = m.advances[xm.order[d]];
    let t = 0;
    for (let k = 0; k < label.length && k < adv.length; k++) {
      const c = t + adv[k] / 2;
      if (t + adv[k] > m.run) break;
      g[0] = cx + c * Slant;
      g[1] = base - c * Slant;
      g[2] = Angle;
      g[3] = 1;
      g[4] = muted;
      out.placedLetterArgs(label.charCodeAt(k));
      t += adv[k];
    }
  }
  out.endGlyphs();
}

/**
 * A label under a point: on y's column (left of the plot, within it down) or on x's slanted labels (above the plot,
 * the column its run leaves from), its member as a `row` hit; undefined elsewhere.
 */
export function hitMemberLabel(
  s: SectionLayout,
  m: MemberLabels,
  ax: AxisTransform,
  ay: AxisTransform | undefined,
  plotTop: number,
  plotHeight: number,
  x: number,
  y: number,
): Hit | undefined {
  const h = m.hit;
  if (ay?.members && x < ax.x0 && x >= s.config.inset && y >= plotTop && y < plotTop + plotHeight) {
    const p = ay.positionAt(y);
    if (p < 0) return undefined;
    return rowHit(h, s, ay.system, p);
  }
  if (ax.members && m.room > 0 && y < plotTop && y >= plotTop - m.room) {
    // the slanted run through the point leaves the plot's top at x − (top − y)
    const p = ax.positionAt(x - (plotTop - Gap - y));
    if (p < 0) return undefined;
    return rowHit(h, s, ax.system, p);
  }
  return undefined;
}

function rowHit(h: MemberLabels['hit'], s: SectionLayout, system: string, position: number): Hit {
  h.layer = s.section.address;
  h.section = s.id;
  h.system = system;
  h.position = position;
  h.target = 'row';
  return h;
}
