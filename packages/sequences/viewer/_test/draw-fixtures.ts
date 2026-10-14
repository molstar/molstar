/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for RN1's draws (seqspec:viewer#extension-points; seqspec:sources#levels-and-the-zoom-gate): a context with
// the kernel renderers on the fixture measurer, a layer drawn into a backend through the conformance kit's composer
// (`LayerFrame`), documents reaching the zoomed-out paths (dense values and a row field under a continuous palette,
// 1,000 distinct values, many emphasized positions), and the readers tests use on a recording: rects by device column,
// texts with their x.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { Pass } from '../src/draw/backend.js';
import { drawGroupLabels } from '../src/draw/numbers.js';
import { DrawView } from '../src/draw/view.js';
import {
  EntryKind,
  KernelRenderers,
  type LayerRenderer,
  type RenderProps,
  type SeqViewerContext,
} from '../src/index.js';
import { LayerFrame, Op, RecordingBackend } from '../src/testing/index.js';
import { opened } from './layout-fixtures.js';

const PROTEIN = 'MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH';
const residues = (n: number, seed = 0) =>
  Array.from({ length: n }, (_, i) => PROTEIN[(i * 7 + seed * 13) % PROTEIN.length]).join('');

/** A context showing `doc` with the kernel renderers (every viewport `width` × `height` at `dpr`). */
export async function drawn(
  doc: SeqViewSpecDoc,
  o: { width?: number; height?: number; dpr?: number; render?: Partial<RenderProps> } = {},
) {
  const t = await opened(doc, { renderers: [...KernelRenderers] as LayerRenderer[], render: o.render });
  for (const ref of t.ctx.viewports.sections.value)
    t.ctx.viewports.get(ref)?.setSize(o.width ?? 1000, o.height ?? 400, o.dpr ?? 1);
  await t.ctx.graph.update();
  return t;
}

/** The layer of a track (by its label) and its section's ref. */
export function layerOf(ctx: SeqViewerContext, label: string): { k: number; section: string } {
  const m = ctx.graph.assembled.value?.view;
  if (!m) throw new Error('no view');
  const r = m.resolved;
  for (let k = 0; k < r.layers.count; k++) {
    const t = r.layers.track[k];
    if (m.set.strings.get(r.tracks.label[t]) !== label) continue;
    return { k, section: m.set.strings.get(r.sections.ref[r.tracks.section[t]]) };
  }
  throw new Error(`no track "${label}"`);
}

/** A track's layer drawn into a fresh recording (or `into`) through the kit's composer; its frame. */
export function record(ctx: SeqViewerContext, label: string, into = new RecordingBackend(), tiny = false) {
  const { k, section } = layerOf(ctx, label);
  const entry = ctx.registry.get('renderers', ctx.graph.assembled.value?.view?.prepared[k].renderer ?? '');
  if (!entry) throw new Error(`no renderer for "${label}"`);
  const f = new LayerFrame(ctx, entry, section, k);
  f.tiny = tiny;
  f.draw(into);
  return { rec: into, frame: f };
}

/** A recording's rects as [x, y, w, h, color, alpha] rows. */
export function rects(r: RecordingBackend): number[][] {
  const b = r.rectItems;
  return Array.from({ length: b.n }, (_, i) => [b.x[i], b.y[i], b.w[i], b.h[i], b.color[i], b.alpha[i]]);
}
/** A recording's texts as { text, x, y, color } with their batch's font and alignment (`TextAlign`). */
export function texts(r: RecordingBackend) {
  const b = r.textItems;
  const out: { text: string; x: number; y: number; color: number; font: number; align: number }[] = [];
  for (let k = 0; k < r.opCount; k++) {
    if (r.ops[4 * k] !== Op.text) continue;
    const [start, count, a] = [r.ops[4 * k + 1], r.ops[4 * k + 2], r.ops[4 * k + 3]];
    for (let i = start; i < start + count; i++)
      out.push({ text: b.text[i], x: b.x[i], y: b.y[i], color: b.color[i], font: r.args[a], align: r.args[a + 1] });
  }
  return out;
}
/** Per device column, the color of the rects covering it (later ones over earlier, as drawn); -2: none. */
export function columns(r: RecordingBackend, x0: number, x1: number, y?: number): number[] {
  const out = new Array<number>(x1 - x0).fill(-2);
  for (const [x, ry, w, h, c] of rects(r))
    for (let p = Math.max(x0, x); p < Math.min(x1, x + w); p++)
      if (y === undefined || (y >= ry && y < ry + h)) out[p - x0] = c;
  return out;
}

/**
 * Zoomed-out cells: `Shade` a `color` layer over `n` dense values (1,000 distinct, Viridis: the mean),
 * `Scores` one over ranges colored by a row field under Viridis (the row field's values), `Codes` one over ranges
 * colored by category (classes), `Levels` the dense values under a categorical palette (by value, past 255 colors),
 * `Marks` emphasis on every 7th position, and letters and a ruler; `blocks` chains.
 */
export function meanDoc(n = 100_000, blocks = 1) {
  const b = createSeqViewSpec({ title: 'Means' });
  const view = b.view({ ref: 'v', label: 'Means', height: 'box' });
  const shades = [];
  const scores = [];
  for (let c = 0; c < blocks; c++) {
    const a = b.data.sequence({ ref: `s${c}`, label: `S${c}`, value: residues(n, c) });
    view.block({ sequence: a });
    const value = Array.from({ length: n }, (_, i) => Math.floor((i * 1000) / n) / 1000);
    shades.push(b.data.annotation({ ref: `shade${c}`, geometry: 'values', target: a, value: { value } }));
    const rows = 400;
    scores.push(
      b.data.annotation({
        ref: `scores${c}`,
        geometry: 'ranges',
        target: a,
        value: {
          start: Array.from({ length: rows }, (_, r) => Math.floor((r * n) / rows)),
          stop: Array.from({ length: rows }, (_, r) => Math.floor(((r + 1) * n) / rows)),
          score: Array.from({ length: rows }, (_, r) => r / rows),
          category: Array.from({ length: rows }, (_, r) => `c${r % 5}`),
        },
      } as never),
    );
  }
  const viridis = {
    theme: 'field_palette',
    params: { field: 'score', palette: { kind: 'continuous', colors: 'Viridis' } },
  };
  const sites = b.data.annotation({
    ref: 'sites',
    geometry: 'ranges',
    target: 's0',
    value: { start: Array.from({ length: Math.floor(n / 7) }, (_, i) => i * 7) },
  } as never);
  view.track({ label: '' }).ruler();
  view.track({ label: 'Letters' }).letters();
  view.track({ label: 'Marks' }).emphasis({ data: sites });
  const shade = {
    theme: 'field_palette',
    params: { field: 'value', palette: { kind: 'continuous', colors: 'Viridis' } },
  };
  view.track({ label: 'Shade' }).color({ data: shades as never, color: shade } as never);
  view.track({ label: 'Scores' }).color({ data: scores as never, color: viridis } as never);
  view
    .track({ label: 'Codes' })
    .color({ data: scores as never, color: { theme: 'field_palette', params: { field: 'category' } } } as never);
  // A categorical palette over the dense values: by value, past 255 colors (the raw path).
  const levels = {
    theme: 'field_palette',
    params: { field: 'value', palette: { kind: 'categorical', colors: 'Category10' } },
  };
  view.track({ label: 'Levels' }).color({ data: shades as never, color: levels } as never);
  return b.getDocument();
}

/** A wrapped section's group labels (`aboveEnd`) drawn for each of its `labels` entries, as the composer will. */
export function recordLabels(ctx: SeqViewerContext, section: string, into = new RecordingBackend()) {
  const l = ctx.layout.get(section)?.current;
  const ax = ctx.layout.axis(section);
  if (!l || !ax) throw new Error('no layout');
  const v = new DrawView();
  into.reset();
  into.beginPass(Pass.main, Math.ceil(l.width * ax.dpr), Math.ceil(l.height * ax.dpr), ax.dpr);
  for (let e = 0; e < l.entryCount; e++) {
    if (l.entryKind[e] !== EntryKind.labels) continue;
    v.from(ax);
    v.backend = into;
    v.skin = ctx.skin.value;
    v.text = ctx.text;
    v.top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - ax.scrollY;
    v.height = l.entryHeight[e];
    const line = l.entryRow[e];
    v.line = line;
    v.onlyBlock = l.lines.block[line];
    v.lineStart = l.lines.start[line];
    v.lineStop = l.lines.stop[line];
    drawGroupLabels(v);
  }
  into.endPass();
  return into;
}

/** A section's layout entry box of a track (CSS px, by its label): its top and height (the first entry). */
export function entryOf(ctx: SeqViewerContext, label: string) {
  const { k, section } = layerOf(ctx, label);
  const l = ctx.layout.get(section)?.current;
  const ax = ctx.layout.axis(section);
  const t = ctx.graph.assembled.value?.view?.resolved.layers.track[k];
  if (!l || !ax) throw new Error('no layout');
  for (let e = 0; e < l.entryCount; e++)
    if (l.entryKind[e] === EntryKind.track && l.entryTrack[e] === t) {
      const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - ax.scrollY;
      return { top, height: l.entryHeight[e], dpr: ax.dpr, left: l.dataLeft, right: l.dataRight };
    }
  throw new Error(`no entry for "${label}"`);
}
