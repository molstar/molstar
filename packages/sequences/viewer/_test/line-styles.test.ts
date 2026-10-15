/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Line styles on ranges (seqspec:format#themes; seqspec:viewer#the-render-backend, seqspec:viewer#the-text-view): each
// item in its style wherever a range draws, in tracks (lanes, brackets), the text's underline lanes and its tracks
// drawn in full, and the SVG export: dashed and dotted a dashed stroke anchored to the content, wavy a sine polyline,
// double two thin lines, box an outline, highlight a translucent fill under the rest (in tracks the entry's rows, in
// text behind the line's letters too). A value that is no line style is drawn solid and said once.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { Dirty, type PathBatch, type RenderProps } from '../src/draw/index.js';
import { renderSvg } from '../src/svg/index.js';
import { Op, type RecordingBackend } from '../src/testing/index.js';
import { composed, recording } from './composer-fixtures.js';
import { drawn, record } from './draw-fixtures.js';
import { opened, preparedRenderers } from './layout-fixtures.js';
import { seqLoci } from './overlay-docs.js';

const STYLES = ['solid', 'dashed', 'dotted', 'wavy', 'double', 'box', 'highlight'];
/** A highlight's alpha: 72/256 of opaque. */
const HL = (255 * 72) >> 8;

/** A protein of 120 and one item per style (a column `line_style`), on a ranges track, the layer as `layer` gives. */
function doc(o: { layout?: 'text'; height?: number; layer?: Record<string, unknown>; styles?: string[] } = {}) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', alphabet: 'protein', value: 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(4) });
  const styles = o.styles ?? STYLES;
  // One group of ten each (in text a line holds whole groups: one piece an item).
  const start = styles.map((_, i) => 1 + 10 * i);
  const r = b.data.annotation({
    ref: 'r',
    geometry: 'ranges',
    target: a,
    value: { start, stop: start.map((x) => x + 8), line_style: styles },
  } as never);
  const v = b.view({ ref: 'v', label: 'V', ...(o.layout ? { layout: o.layout } : {}) });
  if (o.layout) v.track({ ref: 'v.letters', label: 'Letters' }).letters();
  else v.block({ sequence: a });
  const h = o.height ? { height: o.height } : {};
  v.track({ ref: 'v.r', label: 'Styled', ...h }).ranges({
    data: r,
    line_style: { theme: 'field_line_style' },
    ...o.layer,
  } as never);
  return b.getDocument() as SeqViewSpecDoc;
}

/** A recording's stroked subpaths: dash, gap, verbs, width (half px). */
function strokes(p: PathBatch) {
  const out: { dash: number; gap: number; verbs: number; width: number; phase: number }[] = [];
  for (let q = 0; q < p.n; q++) {
    if (p.mode[q] !== 1) continue;
    const v1 = q + 1 < p.n ? p.start[q + 1] : p.nv;
    out.push({ dash: p.dash[q], gap: p.gap[q], verbs: v1 - p.start[q], width: p.width[q], phase: p.phase[q] });
  }
  return out;
}
/** A recording's rects as {w, h, alpha, stroke}. */
function rectList(r: RecordingBackend) {
  const b = r.rectItems;
  return Array.from({ length: b.n }, (_, i) => ({
    y: b.y[i],
    w: b.w[i],
    h: b.h[i],
    alpha: b.alpha[i],
    stroke: b.stroke[i],
  }));
}
/** The styles' marks in a recording: dashed, dotted and wavy strokes, double's thin pairs, boxes, highlights. */
function kinds(r: RecordingBackend) {
  const s = strokes(r.pathItems);
  const rs = rectList(r);
  return {
    dashed: s.filter((x) => x.dash > 0 && x.dash === (3 * x.gap) / 2).length,
    dotted: s.filter((x) => x.dash > 0 && x.dash === x.gap).length,
    wavy: s.filter((x) => x.dash === 0 && x.verbs > 8).length,
    boxes: rs.filter((x) => x.stroke > 0).length,
    highlights: rs.filter((x) => x.alpha === HL).length,
  };
}

describe('line styles in tracks', () => {
  it('each item in its style: the solid bar, dashed, dotted, wavy, two lines, an outline, a highlight under all', async () => {
    const t = await drawn(doc(), { width: 1000, height: 400 });
    const { rec } = record(t.ctx, 'Styled');
    expect(kinds(rec)).toEqual({ dashed: 1, dotted: 1, wavy: 1, boxes: 1, highlights: 1 });
    const rs = rectList(rec);
    // The highlight is drawn first (under the rest) over the entry's rows; the solid bar is the tallest opaque fill.
    expect(rs[0].alpha).toBe(HL);
    expect(rs[0].h).toBeGreaterThan(Math.max(...rs.slice(1).map((x) => x.h)));
    const bar = Math.max(...rs.filter((x) => x.alpha === 255 && x.stroke === 0).map((x) => x.h));
    // Double: two 1 px lines, the bar's band apart.
    const thin = rs.filter((x) => x.alpha === 255 && x.stroke === 0 && x.h === 1);
    expect(thin).toHaveLength(2);
    expect(thin[1].y - thin[0].y).toBe(bar - 1);
    // Dashed and dotted: a line at most 3 px wide, the pattern in its widths.
    const dash = strokes(rec.pathItems).filter((x) => x.dash > 0);
    expect(dash.map((x) => [x.width / 2, x.dash, x.gap])).toEqual([
      [3, 9, 6],
      [3, 3, 3],
    ]);
    t.ctx.dispose();
  });

  it('a constant style for the whole layer; none (null) draws as before', async () => {
    const one = await drawn(doc({ layer: { line_style: { theme: 'uniform', params: { line_style: 'dashed' } } } }));
    expect(kinds(record(one.ctx, 'Styled').rec).dashed).toBe(STYLES.length);
    one.ctx.dispose();
    const none = await drawn(doc({ layer: { line_style: null } }));
    const rec = record(none.ctx, 'Styled').rec;
    expect(rec.pathItems.n).toBe(0);
    expect(rectList(rec).every((x) => x.alpha === 255 && x.stroke === 0)).toBe(true);
    expect(rec.rectItems.n).toBe(STYLES.length);
    none.ctx.dispose();
  });

  it('brackets: the styles on the bracket line, its ticks solid', async () => {
    const t = await drawn(doc({ layer: { display: 'brackets' } }));
    expect(kinds(record(t.ctx, 'Styled').rec)).toEqual({ dashed: 1, dotted: 1, wavy: 1, boxes: 1, highlights: 1 });
    t.ctx.dispose();
  });

  it('a dash pattern is anchored to the content: an item cut by the view keeps the phase it has in full', async () => {
    const t = await drawn(doc({ layer: { line_style: { theme: 'uniform', params: { line_style: 'dashed' } } } }), {
      width: 1000,
    });
    const vp = t.ctx.viewports.get(t.ctx.viewports.sections.value[0]);
    /** The first dashed subpath's start x less its phase, modulo its period, with the view at scrollX `x`. */
    const origin = async (x: number) => {
      vp?.setView(20, x, 'api');
      await t.ctx.graph.update();
      const p = record(t.ctx, 'Styled').rec.pathItems;
      const period = p.dash[0] + p.gap[0];
      return { at: (((p.coords[p.cstart[0]] - p.phase[0]) % period) + period) % period, period };
    };
    const whole = await origin(0);
    // Positions 0 to 2 scrolled off: item [1, 9) starts at position 3's edge, 47 px into its pattern.
    const cut = await origin(47);
    expect(whole.period).toBe(15);
    expect(cut.at).toBe((((whole.at - 47) % 15) + 15) % 15);
    t.ctx.dispose();
  });

  it('a value that is no line style is drawn solid, and the layer says so once', async () => {
    const t = await drawn(doc({ styles: ['dashed', 'zigzag', 'zigzag', 'Wavy'] }));
    const said = t.ctx.issues.current.value.filter((i) => i.code === 'svs.render.line_style');
    expect(said).toHaveLength(1);
    expect(said[0].message).toContain('2 row(s)');
    expect(said[0].message).toContain('"zigzag"');
    const k = kinds(record(t.ctx, 'Styled').rec);
    expect([k.dashed, k.wavy]).toEqual([1, 1]);
    t.ctx.dispose();
  });
});

describe('line styles in the text', () => {
  it('underline lanes: each in its style, a wave over the lane', async () => {
    const t = await drawn(doc({ layout: 'text' }), { width: 600, height: 1200 });
    const { rec } = record(t.ctx, 'Styled');
    expect(kinds(rec)).toEqual({ dashed: 1, dotted: 1, wavy: 1, boxes: 1, highlights: 1 });
    // Dashed on the underline's 3 px.
    expect(strokes(rec.pathItems).find((x) => x.dash > 0)?.width).toBe(6);
    t.ctx.dispose();
  });

  it("highlight and box behind the line's letters, drawn before them", async () => {
    const t = await composed(doc({ layout: 'text' }), { width: 600, height: 1200 });
    const rec = t.frame(recording(), Dirty.all);
    const b = rec.rectItems;
    // The letters' band (12 px letters: 20 px a text row).
    const band = t.entries(t.section).find((e) => e[1] === 'v.letters')?.[4];
    expect(band).toBe(20);
    // On the letters' band (its whole height): the highlight and the box.
    const on: number[] = [];
    for (let i = 0; i < b.n; i++) if (b.h[i] === band && (b.alpha[i] === HL || b.stroke[i] > 0)) on.push(i);
    expect(on.map((i) => (b.stroke[i] > 0 ? 'box' : 'highlight')).sort()).toEqual(['box', 'highlight']);
    // Their rects op comes before the line's letters (a text op).
    let rectOp = -1;
    let textAfter = false;
    for (let k = 0; k < rec.opCount; k++) {
      const [op, start, count] = [rec.ops[4 * k], rec.ops[4 * k + 1], rec.ops[4 * k + 2]];
      if (op === Op.rects && on[0] >= start && on[0] < start + count) rectOp = k;
      if (rectOp >= 0 && op === Op.text) textAfter = true;
    }
    expect(rectOp).toBeGreaterThanOrEqual(0);
    expect(textAfter).toBe(true);
    t.ctx.dispose();
  });

  it('tracks drawn in full: their track forms, the highlight over the track’s rows', async () => {
    const render: Partial<RenderProps> = { wrappedTracks: 'full', wrappedGroupLabels: 'none' };
    const t = await drawn(doc({ layout: 'text', height: 1 }), { width: 600, height: 1200, render });
    const { rec } = record(t.ctx, 'Styled');
    expect(kinds(rec)).toMatchObject({ dashed: 1, dotted: 1, wavy: 1, boxes: 1 });
    expect(kinds(rec).highlights).toBeGreaterThanOrEqual(1);
    t.ctx.dispose();
  });
});

describe('line styles in the SVG export', () => {
  it('writes the dash patterns and the wave; the highlight translucent', async () => {
    for (const layout of [undefined, 'text'] as const) {
      const xml = await renderSvg(doc({ layout }), { width: 800 });
      expect(xml).toMatch(/stroke-dasharray="9 6"/);
      expect(xml).toMatch(/stroke-dasharray="3 3"/);
      expect(xml).toMatch(/fill="none" stroke-width="1"[^>]*d="M[^"]*(L[^"]*){8,}/);
      expect(xml).not.toMatch(/NaN|Infinity|undefined/);
    }
  });
});

describe('bookmark.add with a line_style', () => {
  it('writes the `line_style` column, shown by field_line_style; one not a line style is refused', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ' });
    const notes = b.data.annotation({ ref: 's.bookmarks', geometry: 'ranges', value: [] } as never);
    const sec = b.view({ ref: 'v', label: 'V' }).section({ ref: 's' });
    sec.block({ sequence: a });
    sec.track({ label: 'Notes' }).ranges({ data: notes, line_style: { theme: 'field_line_style' } } as never);
    const { ctx } = await opened(b.getDocument() as SeqViewSpecDoc, { renderers: preparedRenderers() });
    const set = () => ctx.graph.assembled.value?.set as never;
    const rows = () =>
      (ctx.document.shown.value?.root.children ?? [])
        .find((c) => c.kind === 'data')
        ?.children?.find((c) => c.ref === 's.bookmarks')?.params?.value;
    expect(
      ctx.actions.can('bookmark.add', { section: 's', loci: seqLoci(set(), 'a', 2, 6), line_style: 'zigzag' }),
    ).toBe(false);
    const add = { section: 's', loci: seqLoci(set(), 'a', 2, 6), color: '#e76f51', line_style: 'wavy' };
    expect(ctx.actions.run('bookmark.add', add)).toBe(true);
    for (let i = 0; i < 5; i++) await new Promise((r) => setTimeout(r, 0));
    expect(rows()).toEqual([{ target: 'a/unit', start: 2, stop: 6, id: 'b1', color: '#e76f51', line_style: 'wavy' }]);
    expect(ctx.issues.current.value.filter((i) => i.severity !== 'info')).toEqual([]);
    ctx.dispose();
  });
});
