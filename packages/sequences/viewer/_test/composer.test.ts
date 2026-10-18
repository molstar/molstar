/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer on the recording backend (CS1a; seqspec:viewer#the-render-backend): every layer through its registry
// entry (a host's replacement, a namespaced type and its fallback), `tiny` or bars below `minHeight`, a throwing
// renderer omitted and reported once with nothing kept, caches per prepared state and skin; headers (labels cut with an
// ellipsis, the group's and the split track's triangles only where they expand; row labels through `rowModel`; the data
// columns' heads band with no pinned entry, and in each `dataColumnHeads` mode), the wrapped gutter's numbers by the
// group labeler, drawn once per copied frame, the break hairlines from the runs; copies on pan; a stream's version
// redrawing only its layer's entries, in its section only; and presentation: a 30-frame resize presents each frame
// whole, after all its batches, and clears nothing it does not present.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { case13 } from '../../core/_test/format/cases.js';
import { Dirty, PathMode, SectionComposer } from '../src/draw/index.js';
import type { DrawView } from '../src/draw/view.js';
import { DefaultSkin, EntryKind, KernelRenderers, type LayerRenderer, Lod } from '../src/index.js';
import type { StreamView } from '../src/streams.js';
import { Op, type RecordingBackend, StarsFormat, StarsRenderer } from '../src/testing/index.js';
import { COPY_BACK, composed, opArgs, opNames, recording, textsOf } from './composer-fixtures.js';
import { asText, starsDoc } from './conformance-fixtures.js';
import { rects, texts } from './draw-fixtures.js';
import { allTypes, rangesDoc } from './prepare-fixtures.js';

/** A frame's ops with their arguments, rects and texts. */
const snapshot = (r: RecordingBackend) => [opNames(r), Array.from(r.args.subarray(0, r.argCount)), rects(r), texts(r)];
const kernel = (name: string) => (KernelRenderers as LayerRenderer[]).find((r) => r.name === name) as LayerRenderer;
const issues = (t: Awaited<ReturnType<typeof composed>>) =>
  t.ctx.issues.current.value.filter((i) => i.code === 'seqv.renderer.error').map((i) => i.message);
/** A renderer that records its calls (the layer drawn and the view's height), else as `base`. */
function spy(base: LayerRenderer, over: Partial<LayerRenderer> = {}) {
  const calls: string[] = [];
  const r: LayerRenderer = {
    ...base,
    draw: (v: DrawView) => calls.push(`draw ${v.layer} ${v.height}`),
    tiny: (v: DrawView) => calls.push(`tiny ${v.layer}`),
    ...over,
  };
  return { r, calls };
}

describe('layers through their registry entries', () => {
  it("a host's replacement, a namespaced type and, without it, its fallback", async () => {
    const values = spy(kernel('values'));
    const t = await composed(allTypes(), { replace: { values: values.r } });
    t.frame();
    expect(values.calls.length).toBeGreaterThan(0);
    // `org.example:stars` draws through its entry; a context without it, through its fallback's (`points`).
    const stars = spy(StarsRenderer as LayerRenderer);
    const s = await composed(starsDoc(), { registry: [{ svs: StarsFormat, renderers: [stars.r] }] });
    s.frame();
    const points = spy(kernel('points'));
    const f = await composed(starsDoc(), { registry: [{ svs: StarsFormat }], replace: { points: points.r } });
    f.frame();
    expect([stars.calls.length, points.calls.length]).toEqual([1, 1]);
    expect(f.ctx.graph.assembled.value?.view?.prepared.map((p) => p.renderer)).toContain('points');
  });

  it('below its minimum height: `tiny`, else its extents as bars', async () => {
    const tall = spy(kernel('values'), { minHeight: 1000 });
    const t = await composed(allTypes(), { replace: { values: tall.r } });
    t.frame();
    expect(tall.calls.every((c) => c.startsWith('tiny'))).toBe(true);
    const bare = spy(kernel('values'), { minHeight: 1000, tiny: undefined });
    const u = await composed(allTypes(), { replace: { values: bare.r } });
    const r = u.frame();
    expect(bare.calls).toEqual([]);
    // Bars in the faint ink, a quarter down the entry, half its height.
    expect(Array.from(r.rectItems.color.subarray(0, r.rectItems.n))).toContain(-6);
  });

  it('a renderer that throws is omitted, reported once, the frame completes, and the view keeps nothing', async () => {
    const boom = spy(kernel('values'), {
      draw() {
        throw new Error('boom');
      },
      cache: { create: () => ({}) },
    });
    const t = await composed(allTypes(), { replace: { values: boom.r } });
    const r = t.frame();
    expect(issues(t)).toEqual(['layer "v.values-2" threw while drawing and is omitted: boom']);
    expect(opNames(r).at(-1)).toBe('endPass');
    const v = (t.composer as unknown as { v: DrawView }).v;
    expect([v.cache, v.stream]).toEqual([null, null]);
    t.composer.redraw();
    t.frame();
    expect(issues(t).length).toBe(1);
  });

  it('caches per prepared state and skin: made once, rebuilt once on a skin change', async () => {
    let made = 0;
    let disposed = 0;
    const cached = spy(kernel('values'), { cache: { create: () => ++made, dispose: () => disposed++ } });
    const t = await composed(allTypes(), { replace: { values: cached.r } });
    t.frame();
    t.vp?.setView(4, 30, 'api');
    t.frame();
    expect([made, disposed]).toEqual([1, 0]);
    t.ctx.setSkin({ ...DefaultSkin, dark: true });
    await t.ctx.graph.update();
    t.frame();
    t.vp?.panBy(20, 0, 'api');
    t.frame();
    expect([made, disposed, t.composer.stats.cacheBuilds]).toEqual([2, 1, 2]);
  });
});

describe('headers', () => {
  it('track labels cut with an ellipsis to the header; titles bold', async () => {
    const t = await composed(allTypes(), { render: { headerWidth: 60 } });
    const texts = textsOf(t.frame());
    expect(texts).toContain('Letters');
    expect(texts.some((s) => s.endsWith('…'))).toBe(true);
    expect(texts.every((s) => !s.startsWith('Secondary') || s.length < 9)).toBe(true);
  });

  it("a split track's triangle only where it expands: a scrolled entry whose split has values", async () => {
    const t = await composed(case13(), { width: 1000, height: 900, section: 'features.main' });
    const r = t.frame();
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error('no layout');
    // One triangle per group header and per split label entry (row -2) in view, none on a pinned entry.
    let expected = 0;
    const m = t.ctx.graph.assembled.value?.view?.resolved;
    for (let e = l.pinnedCount; e < l.entryCount; e++)
      if (l.entryKind[e] === EntryKind.group || (l.entryRow[e] === -2 && m && m.tracks.split[l.entryTrack[e]] >= 0))
        if (l.scrolledTop + l.entryTop[e] < l.height) expected++;
    expect(expected).toBeGreaterThan(0);
    expect(r.count(Op.paths)).toBe(expected);
  });

  it("a group's header: regular weight, its label where a track's starts, a thin stroked chevron after its text", async () => {
    const t = await composed(case13(), { width: 1000, height: 900, section: 'features.main' });
    const r = t.frame();
    const m = t.ctx.graph.assembled.value?.view;
    const l = t.ctx.layout.get(t.section)?.current;
    if (!m || !l) throw new Error('no view');
    const e = Array.from({ length: l.entryCount }, (_, k) => k).find(
      (k) =>
        l.entryKind[k] === EntryKind.group && l.entryGroup[k] >= 0 && m.resolved.groups.parent[l.entryGroup[k]] < 0,
    ) as number;
    const name = m.set.strings.get(m.resolved.groups.label[l.entryGroup[e]]);
    const y = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2;
    const label = texts(r).find((x) => x.text === name && Math.abs(x.y - y) < 2);
    // A top-level label starts at 10 px (no inset for a triangle), in the weight and size of a track's label.
    const track = texts(r).find((x) => x.text === 'Sequence' || x.text === 'Chain');
    expect([label?.x, label?.font]).toEqual([10, track?.font]);
    const font = label?.font as number;
    expect(t.ctx.text.fontString(font)).toMatch(/\b400\b/);
    // One stroked path of three points whose box is centered after the text (CHEVRON_GAP + CHEVRON_HALF past it).
    const batch = r.pathItems;
    const p = Array.from({ length: batch.n }, (_, k) => k).find(
      (k) => Math.abs(batch.coords[batch.cstart[k] + 1] - y) < 6,
    );
    expect(p).toBeDefined();
    const c = batch.coords.subarray(batch.cstart[p as number], batch.cstart[p as number] + 6);
    const center = (Math.min(c[0], c[2], c[4]) + Math.max(c[0], c[2], c[4])) / 2;
    expect([batch.mode[p as number], center]).toEqual([PathMode.stroke, 10 + t.ctx.text.label(name, font).width + 8]);
  });

  it("a rows layer's labels through its `rowModel`, and the data columns' heads with no pinned entry", async () => {
    const rows = kernel('rows');
    const named = { ...rows, rowModel: { label: (_s: unknown, row: number) => `row ${row}` } } as LayerRenderer;
    const b = createSeqViewSpec();
    const sto = ['# STOCKHOLM 1.0', 'a1 MKV-LA', 'a2 MKVQLA', 'a3 MRV-LA', '//', ''].join('\n');
    const aln = b.data.alignment({ ref: 'aln', source: { value: sto, format: 'stockholm' } });
    const meta = b.data.annotation({
      ref: 'meta',
      geometry: 'rows',
      target: aln,
      value: { row: [0, 1, 2], score: [1, 2, 3] },
    } as never);
    const v = b.view({ ref: 'f', label: 'F', height: 'box' });
    v.block({ alignment: aln });
    v.track({ label: 'Rows' }).rows({ data_columns: [{ data: meta, field: 'score', label: 'Score' }] } as never);
    const t = await composed(b.getDocument(), { replace: { rows: named }, height: 300 });
    const l = t.ctx.layout.get(t.section)?.current;
    const band = (l?.scrolledTop ?? 0) - (l?.navHeight ?? 0);
    expect([l?.pinnedCount, band]).toEqual([0, t.ctx.layout.props(t.section).rowHeight]);
    const texts = textsOf(t.frame());
    expect(texts).toEqual(expect.arrayContaining(['row 0', 'row 1', 'row 2', 'Score']));
  });

  it("the data columns' heads directly above the rows with several pinned entries: on the last one's line", async () => {
    const b = createSeqViewSpec();
    const sto = ['# STOCKHOLM 1.0', 'a1 MKV-LA', 'a2 MKVQLA', 'a3 MRV-LA', '//', ''].join('\n');
    const aln = b.data.alignment({ ref: 'aln', source: { value: sto, format: 'stockholm' } });
    const meta = b.data.annotation({
      ref: 'meta',
      geometry: 'rows',
      target: aln,
      value: { row: [0, 1, 2], score: [1, 2, 3] },
    } as never);
    const v = b.view({ ref: 'f', label: 'F', height: 'box' });
    v.block({ alignment: aln });
    v.track({ label: '' }).ruler();
    v.track({ label: 'Second' }).ruler();
    v.track({ label: 'Third' }).ruler();
    v.track({ label: 'Rows' }).rows({ data_columns: [{ data: meta, field: 'score', label: 'Score' }] } as never);
    const t = await composed(b.getDocument(), { height: 300 });
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error('no layout');
    expect(l.pinnedCount).toBe(3);
    const r = t.frame();
    const i = textsOf(r).indexOf('Score');
    expect(i).toBeGreaterThanOrEqual(0);
    // The last pinned entry's center (dpr 1, the label's half pixel), not the middle of the three.
    const last = l.pinnedCount - 1;
    expect(l.entryTop[last] + l.entryHeight[last]).toBe(l.scrolledTop);
    expect(r.textItems.y[i]).toBe(Math.round(l.entryTop[last] + l.entryHeight[last] / 2 + 0.5));
    expect(r.textItems.y[i]).toBeGreaterThan(Math.round((l.navHeight + l.scrolledTop) / 2 + 0.5));
  });

  /** A rows layer with a data column, under three pinned rulers or none. */
  const headsDoc = (pinned: boolean) => {
    const b = createSeqViewSpec();
    const sto = ['# STOCKHOLM 1.0', 'a1 MKV-LA', 'a2 MKVQLA', 'a3 MRV-LA', '//', ''].join('\n');
    const aln = b.data.alignment({ ref: 'aln', source: { value: sto, format: 'stockholm' } });
    const meta = b.data.annotation({
      ref: 'meta',
      geometry: 'rows',
      target: aln,
      value: { row: [0, 1, 2], score: [1, 2, 3] },
    } as never);
    const v = b.view({ ref: 'f', label: 'F', height: 'box' });
    v.block({ alignment: aln });
    if (pinned) for (const label of ['', 'Second', 'Third']) v.track({ label }).ruler();
    v.track({ label: 'Rows' }).rows({ data_columns: [{ data: meta, field: 'score', label: 'Score' }] } as never);
    return b.getDocument();
  };
  /** The heads' text y, the rows' top and the last pinned entry's bottom (the navigator's without one), a row's height. */
  const headsAt = async (pinned: boolean, mode: 'shared' | 'ownRow') => {
    const t = await composed(headsDoc(pinned), { height: 300, render: { dataColumnHeads: mode } });
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error('no layout');
    const r = t.frame();
    const i = textsOf(r).indexOf('Score');
    const last = l.pinnedCount - 1;
    const bottom = last >= 0 ? l.entryTop[last] + l.entryHeight[last] : l.navHeight;
    return { l, y: i >= 0 ? r.textItems.y[i] : -1, bottom, row: t.ctx.layout.props(t.section).rowHeight };
  };

  it('`ownRow`: the heads on a row of their own between the last pinned entry and the rows, never its line', async () => {
    const own = await headsAt(true, 'ownRow');
    expect(own.l.pinnedCount).toBe(3);
    expect(own.l.scrolledTop).toBe(own.bottom + own.row);
    // Centered in the reserved row (dpr 1, the label's half pixel), below the last pinned entry's bottom.
    expect(own.y).toBe(Math.round(own.bottom + own.row / 2 + 0.5));
    expect(own.y).toBeGreaterThan(own.bottom);
    // `shared` over the same document: no row reserved, the heads on the last pinned line, above its bottom.
    const shared = await headsAt(true, 'shared');
    expect(shared.l.scrolledTop).toBe(shared.bottom);
    expect(shared.y).toBeLessThan(shared.bottom);
    expect(own.l.scrolledTop - shared.l.scrolledTop).toBe(own.row);
  });

  it('with nothing pinned both modes reserve the one row under the navigator and draw the heads there', async () => {
    for (const mode of ['shared', 'ownRow'] as const) {
      const h = await headsAt(false, mode);
      expect([h.l.pinnedCount, h.l.scrolledTop - h.l.navHeight]).toEqual([0, h.row]);
      expect(h.y).toBe(Math.round(h.l.navHeight + h.row / 2 + 0.5));
    }
  });
});

describe('gutters and breaks', () => {
  /** A protein numbered from 1001 and one without a numbering, wrapped, its gutter numbering the lines. */
  function wrappedDoc() {
    const b = createSeqViewSpec();
    const a = b.data
      .sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(30) })
      .numbering({ name: 'auth', method: 'offset', first: 1001 });
    const c = b.data.sequence({ ref: 'c', value: 'MKVLATGHEK'.repeat(12) });
    const v = b.view({ ref: 'v', label: 'V', layout: 'text', numbering: 'auth', height: 'box' });
    v.block({ sequence: a });
    v.block({ sequence: c });
    v.track({ label: 'L' }).letters();
    return b.getDocument();
  }

  it("wrapped: each line's first position by the group labeler, the numbering's or ordinal", async () => {
    const t = await composed(wrappedDoc(), { width: 400, height: 2000, render: { wrappedGroupLabels: 'none' } });
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error('no layout');
    const numbers = textsOf(t.frame()).filter((s) => /^\d+$/.test(s));
    const want: string[] = [];
    for (let k = 0; k < l.lines.count; k++) want.push(String(l.lines.start[k] + (l.lines.block[k] === 0 ? 1001 : 1)));
    expect(numbers).toEqual(want);
  });

  it("a copied wrapped scroll draws each line's number once, and the frame equals a full redraw's", async () => {
    // At 15 px a letter: three copies over these scrolls.
    const t = await composed(wrappedDoc(), {
      width: 400,
      height: 300,
      render: { wrappedGroupLabels: 'none', wrappedFontPx: 15 },
    });
    const be = recording(COPY_BACK);
    t.frame(be);
    for (const dy of [37, -20, 55]) {
      t.vp?.panBy(0, dy, 'api');
      const r = t.frame(be);
      expect(r.count(Op.copy)).toBe(3);
      const numbers = textsOf(r).filter((s) => /^\d+$/.test(s));
      expect(new Set(numbers).size).toBe(numbers.length);
    }
  });

  it('break hairlines where a filter skips positions, from the shown runs', async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: 's', value: 'MKVLAAGIVGLLLAQRSTVW'.repeat(3) });
    const site = b.data.annotation({
      ref: 'site',
      geometry: 'ranges',
      target: s,
      value: { start: [5, 30], stop: [15, 40] },
    });
    b.view({ ref: 'f', label: 'F', positions: site }).block({ sequence: s }).track({ label: 'L' }).letters();
    const t = await composed(b.getDocument(), { width: 600 });
    const r = t.frame();
    const ax = t.ctx.layout.axis(t.section);
    const l = t.ctx.layout.get(t.section)?.current;
    if (!ax || !l) throw new Error('no axis');
    const tall = [];
    for (let i = 0; i < r.rectItems.n; i++)
      if (r.rectItems.color[i] === -6 && r.rectItems.h[i] >= Math.floor((l.height - l.navHeight) * ax.dpr) - 1)
        tall.push(r.rectItems.x[i]);
    // One break, between the filter's two runs: at shown 10.
    expect(tall).toEqual([Math.round(ax.x(0, 10 as never) * ax.dpr) - 0]);
  });
});

describe('frames', () => {
  it('a pan at the same zoom copies the data, the pinned band and the header, then draws the strip', async () => {
    const t = await composed(allTypes({ length: 600 }), { width: 900 });
    t.vp?.setView(6, 100, 'api');
    const be = recording(COPY_BACK);
    t.frame(be);
    t.vp?.panBy(40, 0, 'api');
    const r = t.frame(be);
    const names = opNames(r);
    expect(names.filter((n) => n === 'copy').length).toBe(3);
    expect(names.indexOf('present')).toBe(names.length - 2);
    // Nothing moved: nothing drawn.
    expect(t.frame(be).opCount).toBe(0);
    expect(t.composer.stats.copies).toBe(1);
  });

  it("at letter zoom a pan's strip is widened to a cell edge: no letter cut where strip and copy meet", async () => {
    const t = await composed(allTypes({ length: 600 }), { width: 900, dpr: 2 });
    t.vp?.setView(16, 100, 'api');
    const be = recording(COPY_BACK);
    t.frame(be);
    const v = (t.composer as unknown as { v: DrawView }).v;
    // The strip (the third clear after the copies): its inner edge (device px) on a cell's edge, wider than the move.
    const strip = (r: RecordingBackend, dx: number) => {
      const names = opNames(r);
      let k = names.lastIndexOf('copy');
      for (let n = 0; n < 3; n++) k = names.indexOf('clear', k + 1);
      const [x, , w] = opArgs(r, k, 4);
      const inner = dx > 0 ? x : x + w;
      const cell = Math.round((inner - v.xd(0, 0)) / 32);
      return [w >= Math.abs(dx) * 2, inner === v.xd(0, cell), v.lod];
    };
    t.vp?.panBy(5, 0, 'api');
    expect(strip(t.frame(be), 5)).toEqual([true, true, Lod.letters]);
    t.vp?.panBy(-7, 0, 'api');
    expect(strip(t.frame(be), -7)).toEqual([true, true, Lod.letters]);
    expect(t.composer.stats.copies).toBe(2);
    // Moves that leave the strip's inner edge in either half of a cell, each way (LC-39: a cell edge rounded to the
    // nearer one, not beyond the move, passed the two above).
    for (const dx of [3, 6, 9, 13, 2, -3, -6, -9, -13, -2]) {
      t.vp?.panBy(dx, 0, 'api');
      expect(strip(t.frame(be), dx), `${dx}`).toEqual([true, true, Lod.letters]);
    }
  });

  it("a stream's arrival redraws only its layer's entries, in its own section only", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(6) });
    const vals = b.data.annotation({ ref: 'vals', geometry: 'values', target: a, value: { value: [0.5, 1, 0.2] } });
    const g = b.view({ ref: 'g', label: 'G' });
    const one = g.section({ ref: 'g.one' }).block({ sequence: a });
    one.track({ label: 'L' }).letters();
    one.track({ label: 'V' }).values({ ref: 'g.vals', data: vals });
    g.section({ ref: 'g.two' }).block({ sequence: a }).track({ label: 'L' }).letters();
    const t = await composed(b.getDocument(), { section: 'g.one' });
    const other = new SectionComposer(t.ctx, 'g.two');
    const be = recording(COPY_BACK);
    const quiet = recording(COPY_BACK);
    // A stream on the values layer, as `ctx.streams.get` gives one for a layer held by region.
    const fake = { version: 0 } as StreamView;
    const get = t.ctx.streams.get.bind(t.ctx.streams);
    (t.ctx.streams as { get: typeof get }).get = (ref: string) => (ref === 'g.vals' ? fake : get(ref));
    t.frame(be);
    other.frame(quiet, Dirty.tracks);
    (fake as { version: number }).version++;
    const r = t.frame(be);
    quiet.reset();
    other.frame(quiet, Dirty.tracks);
    // One entry's data box cleared and drawn again, no copy; the other section draws nothing.
    const l = t.ctx.layout.get(t.section)?.current;
    expect([r.count(Op.clear), r.count(Op.copy), quiet.opCount]).toEqual([1, 0, 0]);
    expect(opArgs(r, opNames(r).indexOf('clear'), 4)[2]).toBe(
      Math.round(l?.dataRight ?? 0) - Math.round(l?.dataLeft ?? 0),
    );
  });

  it('a 30-frame resize presents each frame whole, after all its batches, and clears nothing it does not present', async () => {
    // The plan's budget of 1 ms does not apply: no frame is sliced (measured in CS1a).
    const t = await composed(asText(allTypes()), { width: 900, height: 600 });
    const be = recording(COPY_BACK);
    for (let f = 0; f < 30; f++) {
      const w = 900 - f * 13;
      for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(w, 600 - f * 7, 1);
      const r = t.frame(be);
      const names = opNames(r);
      expect(names.slice(-2), `frame ${f}`).toEqual(['present', 'endPass']);
      expect(opArgs(r, names.length - 2, 4), `frame ${f}`).toEqual([0, 0, Math.ceil(w), Math.ceil(600 - f * 7)]);
      expect(names.filter((n) => n === 'clear').length, `frame ${f}`).toBe(1);
      // What is presented equals a fresh composer's frame of the same view: every entry's batches, all of them.
      const fresh = new SectionComposer(t.ctx, t.section);
      const ref = recording(COPY_BACK);
      fresh.frame(ref, Dirty.tracks);
      expect(snapshot(r), `frame ${f}`).toEqual(snapshot(ref));
    }
  });
});

describe("the navigator's passes", () => {
  /** The passes a frame began, by `Pass`, each with its op count. */
  const passes = (r: RecordingBackend) => {
    const out: [number, number][] = [];
    const names = opNames(r);
    names.forEach((n, k) => {
      if (n === 'beginPass') out.push([opArgs(r, k, 1)[0], 0]);
      else if (out.length) out[out.length - 1][1]++;
    });
    return out;
  };

  it('its row on the nav pass (once per layout and size), its window on the nav overlay', async () => {
    const t = await composed(allTypes({ length: 600 }));
    t.vp?.setView(4, 100, 'api');
    const r = t.frame(recording(), Dirty.all);
    const p = passes(r);
    expect(p.map(([pass]) => pass)).toEqual([0, 2, 1, 3]);
    // The nav pass and its overlay hold batches past their clear and end; the overlay is CS1b's.
    expect([p[1][1] > 2, p[3][1] > 2]).toEqual([true, true]);
    // The axis is drawn as at the fit: a pan redraws the window only.
    t.vp?.panBy(30, 0, 'api');
    expect(passes(t.frame(recording(), Dirty.all)).map(([pass]) => pass)).toEqual([0, 1, 3]);
  });

  it("the rulerRail on the overlay, along the ruler's bottom border", async () => {
    const rail = { name: 'lab:rail', geometry: 'linear', navigator: 'rulerRail' } as const;
    const t = await composed(allTypes({ length: 600 }), { registry: [{ presets: [rail] }] });
    t.ctx.ui.setStyle(t.section, rail.name);
    await t.ctx.graph.update();
    t.vp?.setView(4, 100, 'api');
    const r = t.frame(recording(), Dirty.all);
    const p = passes(r);
    expect(p.map(([pass]) => pass)).toEqual([0, 1]);
    expect(p[1][1]).toBeGreaterThan(2);
  });
});

describe('other paths', () => {
  it('a wrapped section draws group labels above each line through `drawGroupLabels`', async () => {
    const t = await composed(asText(rangesDoc([[2, 9]], 120)), { width: 300, height: 1500 });
    const texts = textsOf(t.frame());
    expect(texts).toContain('10');
  });
});
