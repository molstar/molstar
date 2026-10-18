/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `renderSvg` and the SVG side of the tee (SV1a): a case rendered in Node equals the same case on the recording backend,
// element for element (rectangles and texts, by position and color); a document or a context rendered whole, sections
// stacked; the overlay and navigator on request; the conformance kit's SVG check catching a backend that drops a call;
// a live context left as it was (zoom, scroll, row scale, size, measurer, families, skin, view; no viewport committed;
// resized and put back in one task; paused while the export shows another view, skin or measurer); and a genome-sized document (2 × 200,000 positions) in bounded time and size, zoomed out and at letters.
// A section on rings (CG1b): a square figure of arcs, every item the recording saw written, marks in the overlay.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { case1, case13, case14 } from '../../../core/_test/format/cases.js';
import { itWithin } from '../../../core/_test/model/timing.js';
import { ChannelMarks } from '../../src/channels.js';
import { Dirty, SectionComposer } from '../../src/draw/index.js';
import { drawResources } from '../../src/draw/resources.js';
import { VerticalNav } from '../../src/draw/vnav.js';
import { Hover } from '../../src/hover.js';
import { DarkSkin, DefaultSkin, KernelRenderers, type LayerRenderer, type SeqViewerContext } from '../../src/index.js';
import { DefaultPresets } from '../../src/presets.js';
import { cssColor, hex } from '../../src/svg/format.js';
import { DefaultFontMetrics, metricsMeasurer, renderSvg, SvgBackend } from '../../src/svg/index.js';
import { svgDiff, TeeBackend } from '../../src/testing/index.js';
import type { Viewport } from '../../src/viewport.js';
import { composed, opNames, recording } from '../composer-fixtures.js';
import { syntheticFamily } from '../families.js';
import { opened } from '../layout-fixtures.js';
import { LAZY_COLUMNS, LAZY_URI, lazyStoryLoaders } from '../lazy-story-fixtures.js';
import { allTypes, scaleDoc } from '../prepare-fixtures.js';
import { all, parseXml, rectsOf } from './xml.js';

const kernel = [...KernelRenderers] as LayerRenderer[];
const hund = (v: number) => Math.round(v * 100) / 100;

describe.each([1, 2])('a case on the recording backend and as SVG at DPR %i', (dpr) => {
  it('every rectangle and text element, by position and color; counts of the rest', async () => {
    const t = await composed(allTypes(), { width: 900, height: 420, dpr });
    const res = drawResources(t.ctx);
    const rec = recording();
    const svg = new SvgBackend(res, { overlay: true, nav: true, trace: true });
    const composer = new SectionComposer(t.ctx, t.section);
    composer.frame(new TeeBackend(rec, svg), Dirty.all);
    expect(svgDiff(rec, svg)).toBeNull();

    const color = (c: number) => (c < -1 ? cssColor(res.token(c))[0] : hex(c));
    // The recorded rectangles, outlines as the four edges the writer makes of them.
    const want: string[] = [];
    const r = rec.rectItems;
    for (let i = 0; i < r.n; i++) {
      if (r.color[i] === -1) continue;
      const [x, y, w, h, s] = [r.x[i], r.y[i], r.w[i], r.h[i], r.stroke[i]];
      const cut =
        s === 0 || 2 * s >= w || 2 * s >= h
          ? [[x, y, w, h]]
          : [
              [x, y, w, s],
              [x, y + h - s, w, s],
              [x, y + s, s, h - 2 * s],
              [x + w - s, y + s, s, h - 2 * s],
            ];
      for (const [a, b, c, d] of cut)
        want.push(`${color(r.color[i])} ${hund(a / dpr)} ${hund(b / dpr)} ${hund(c / dpr)} ${hund(d / dpr)}`);
    }
    const doc = parseXml(await svg.toString());
    // The passes, not the defs (a layer's rectangles are recorded in its own recording).
    const body = { ...doc, children: doc.children.filter((c) => c.name !== 'defs') };
    const got: string[] = [];
    for (const p of all(body, 'path'))
      for (const [x, y, w, h] of rectsOf(p.attrs.d ?? ''))
        if (p.attrs['shape-rendering']) got.push(`${p.attrs.fill} ${x} ${y} ${w} ${h}`);
    expect(want.length).toBeGreaterThan(200);
    // Edges are rounded one by one, so a width may differ from the rounded width by a hundredth: compare the edges.
    const edges = (list: string[]) =>
      list
        .map((e) => e.split(' '))
        .map(([c, x, y, w, h]) => `${c} ${x} ${y} ${hund(+x + +w)} ${hund(+y + +h)}`)
        .sort();
    expect(got.length).toBe(want.length);
    const [eg, ew] = [edges(got), edges(want)];
    expect(eg.filter((e, i) => e !== ew[i]).length).toBe(0);

    // Texts: the string and its anchor.
    const tx = rec.textItems;
    const wantText: string[] = [];
    for (let i = 0; i < tx.n; i++)
      if (tx.color[i] !== -1) wantText.push(`${tx.text[i]} ${hund(tx.x[i] / dpr)} ${hund(tx.y[i] / dpr)}`);
    const gotText = all(body, 'text').map((e) => {
      const at = e.attrs.transform
        ? /translate\(([-\d.]+) ([-\d.]+)\)/.exec(e.attrs.transform)?.slice(1)
        : [e.attrs.x, e.attrs.y];
      return `${e.text} ${at?.join(' ')}`;
    });
    expect(wantText.length).toBeGreaterThan(20);
    expect(gotText.sort()).toEqual(wantText.sort());
    // Everything else by count (the kit's check, per kind).
    expect(svg.stats).toMatchObject({ sectors: 0, nonFinite: 0 });
    expect(svg.stats.paths).toBeGreaterThan(0);
    expect(svg.stats.glyphs).toBeGreaterThan(0);
  });
});

describe("the kit's SVG check", () => {
  it('is null for the same calls, and names the first call the SVG side missed or the item it did not write', async () => {
    const t = await composed(allTypes(), { width: 600, height: 300 });
    const res = drawResources(t.ctx);
    const rec = recording();
    const svg = new SvgBackend(res, { trace: true, overlay: true, nav: true });
    new SectionComposer(t.ctx, t.section).frame(new TeeBackend(rec, svg), Dirty.all);
    expect(svgDiff(rec, svg)).toBeNull();
    const missed = svg.trace?.findIndex((l) => l.startsWith('text')) ?? -1;
    svg.trace?.splice(missed, 1);
    expect(svgDiff(rec, svg)).toMatch(/^svg call \d+: .*, recorded text /);
    svg.reset();
    expect(svgDiff(rec, svg)).toMatch(/^svg call 0: \(none\), recorded beginPass/);
    const wrote = new SvgBackend(res, { trace: true, overlay: true, nav: true });
    rec.replay(wrote);
    expect(svgDiff(rec, wrote)).toBeNull();
    wrote.stats.rects--;
    expect(svgDiff(rec, wrote)).toMatch(/^svg wrote \d+ rects, recorded \d+/);
  });
});

describe('renderSvg', () => {
  it("a document: sections stacked in one well-formed SVG, sized by the width and each section's content", async () => {
    const xml = await renderSvg(allTypes(), { width: 800, gap: 10, background: true, title: 'All types' });
    const svg = parseXml(xml);
    const sections = all(svg, 'g').filter((g) => g.attrs.class === 'seqv-section');
    expect(sections.length).toBeGreaterThanOrEqual(1);
    expect(svg.attrs.width).toBe('800');
    expect(Number(svg.attrs.height)).toBeGreaterThan(150);
    expect(svg.attrs.viewBox).toBe(`0 0 800 ${svg.attrs.height}`);
    expect(all(svg, 'text').map((e) => e.text)).toContain('Letters');
    // The main pass only: no overlay or navigator unless asked.
    expect(all(svg, 'g').map((g) => g.attrs.class)).not.toContain('seqv-overlay');
    expect(xml).not.toMatch(/NaN|Infinity|undefined/);
  });

  it('a second section sits below the first, at the gap', async () => {
    const t = await opened(case13(), { renderers: kernel });
    const refs = t.ctx.viewports.sections.value;
    const svg = parseXml(await renderSvg(t.ctx, { width: 700, gap: 12, sections: [refs[0], refs[0]] }));
    const [a, b] = all(svg, 'g').filter((g) => g.attrs.class === 'seqv-section');
    const h = Number(svg.attrs.height);
    expect(a.attrs.transform).toBe('translate(0 0)');
    const y = Number(/translate\(0 ([\d.]+)\)/.exec(b.attrs.transform)?.[1]);
    expect(y).toBeGreaterThan(50);
    expect(h).toBeCloseTo(2 * y - 12, 0);
    t.ctx.dispose();
  });

  it("the overlay and the navigator on request; the scale sets the file's size and the layout's DPR", async () => {
    const t = await opened(allTypes(), { renderers: kernel });
    const vp = t.ctx.viewports.get(t.ctx.viewports.sections.value[0]);
    const before = [vp?.live.current.width, vp?.live.current.height, vp?.live.current.dpr];
    const plain = parseXml(await renderSvg(t.ctx, { width: 600 }));
    // A context passed in keeps its size and DPR.
    expect([vp?.live.current.width, vp?.live.current.height, vp?.live.current.dpr]).toEqual(before);
    const full = parseXml(await renderSvg(t.ctx, { width: 600, overlay: true, nav: true, scale: 2 }));
    expect(
      all(plain, 'g')
        .map((g) => g.attrs.class)
        .filter((c) => c?.startsWith('seqv-') && c !== 'seqv-section'),
    ).toEqual(['seqv-main']);
    expect(all(full, 'g').map((g) => g.attrs.class)).toEqual(
      expect.arrayContaining(['seqv-main', 'seqv-overlay', 'seqv-nav']),
    );
    expect(full.attrs.width).toBe('1200');
    expect(full.attrs.viewBox.startsWith('0 0 600 ')).toBe(true);
    // The navigator is the same section's, drawn over its top rows.
    expect(all(full, 'path').length).toBeGreaterThan(all(plain, 'path').length);
    t.ctx.dispose();
  });

  it("a section on rings: a square figure of arcs and the ruler's numbers, every item written, marks in the overlay", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: 'p', value: 'ACGTACGTTA'.repeat(24), alphabet: 'dna', cyclic: true });
    const f = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: s,
      value: { start: [200, 40], stop: [30, 112] },
    });
    const v = b.view({ ref: 'v', label: 'V', style: 'circular' }).block({ sequence: s });
    v.track({ label: '' }).ruler();
    v.track({ label: 'S' }).letters();
    v.track({ label: 'F' }).ranges({ data: f });
    const doc = b.getDocument();
    const xml = await renderSvg(doc, { width: 480, background: true });
    const svg = parseXml(xml);
    // As tall as it is wide (within [320, 560] px), the ring centered in it.
    expect([svg.attrs.width, svg.attrs.height]).toEqual(['480', '480']);
    expect(all(svg, 'path').filter((p) => /A/.test(p.attrs.d ?? '')).length).toBeGreaterThan(0);
    // Every number round the ring, 240 beside 10 at the top (a ring has no edge to cut a label).
    const numbers = [...xml.matchAll(/<text[^>]*>(\d+)</g)].map((m) => Number(m[1])).sort((x, y) => x - y);
    expect(numbers).toEqual(Array.from({ length: 24 }, (_, i) => 10 * (i + 1)));
    expect(xml).not.toMatch(/NaN|Infinity|undefined/);
    // The SVG side of a tee wrote every item of the ring's frame, overlay included.
    const t = await composed(doc, {
      width: 480,
      height: 480,
      registry: [{ presets: DefaultPresets.filter((p) => p.name === 'circular') }],
    });
    t.ctx.selection.set({ systems: ['seqstar:p/unit'], runs: Int32Array.of(0, 225, 240, 0, 0, 12) }, 'api');
    const tee = new TeeBackend(recording(), new SvgBackend(drawResources(t.ctx), { overlay: true, trace: true }));
    t.composer.frame(tee, Dirty.tracks | Dirty.overlay);
    expect(svgDiff(tee.a as never, tee.b as SvgBackend)).toBeNull();
    const lit = parseXml(await renderSvg(t.ctx, { width: 480, overlay: true }));
    const overlay = all(lit, 'g').find((g) => g.attrs.class === 'seqv-overlay');
    expect(overlay && all(overlay, 'path').length).toBeGreaterThan(0);
    t.ctx.dispose();
  });

  it('rejects a section the view does not have', async () => {
    const t = await opened(allTypes(), { renderers: kernel });
    await expect(renderSvg(t.ctx, { width: 600, sections: ['nope'] })).rejects.toThrow('no section "nope"');
    t.ctx.dispose();
  });
});

describe("a document's streams (VB-01, VB-16)", () => {
  /** The lazy story's document, `rows` rows tall: its alignment streamed by tiles and labels. */
  const lazy = (rows: number) => {
    const b = createSeqViewSpec({ title: 'lazy' });
    const aln = b.data.alignment({
      ref: 'lazy',
      row_count: rows,
      column_count: LAZY_COLUMNS,
      source: { uri: LAZY_URI, protocol: 'seqstar-tiles/1' },
    } as never);
    const v = b.view({ ref: 'rows', label: 'Lazy', height: 'box' });
    v.block({ alignment: aln });
    v.track({ label: 'Rows' }).rows();
    return b.getDocument();
  };
  const labels = (svg: string) => (svg.match(/>SEQ\d+</g) ?? []).length;

  it('settle at the size they are drawn: every row of an 18,000 px alignment labeled and drawn', async () => {
    const svg = await renderSvg(lazy(1000), { width: 600, loaders: lazyStoryLoaders() });
    expect(Number(/height="([\d.]+)"/.exec(svg)?.[1])).toBeGreaterThan(15_000);
    expect(labels(svg)).toBe(1000);
    expect(svg.split('<image').length - 1).toBeGreaterThan(0);
  });

  it('a `timeout` draws what stands, and says the streams had not settled', async () => {
    let unsettled = 0;
    const never = lazyStoryLoaders(() => new Promise<Uint8Array>(() => {}));
    const svg = await renderSvg(lazy(200), {
      width: 600,
      loaders: never,
      timeout: 50,
      unsettled: () => unsettled++,
    });
    expect(unsettled).toBe(1);
    expect(svg.startsWith('<svg')).toBe(true);
  });

  it("`timeout` bounds a static source's load too; `timeout: 0` on a document with nothing to wait for is settled (LB-34)", async () => {
    // A FASTA source whose read never ends: the export draws what stands (nothing loaded) at the deadline.
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: 'f', host_file: 's.fasta', format: 'fasta' });
    const seq = b.data.sequence({ ref: 's', alphabet: 'protein', source_ref: src, record: 's' });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: seq }).track({ label: 'S' }).letters();
    let unsettled = 0;
    const svg = await renderSvg(b.getDocument(), {
      width: 600,
      loaders: { load: () => new Promise<Uint8Array>(() => {}) },
      timeout: 100,
      unsettled: () => unsettled++,
    });
    expect([unsettled, svg.startsWith('<svg')]).toEqual([1, true]);
    unsettled = 0;
    await renderSvg(case14(), { width: 600, timeout: 0, unsettled: () => unsettled++ });
    expect(unsettled).toBe(0);
  });
});

describe('the vertical overview in an export (`VerticalNav`)', () => {
  it.each([1, 2])(
    'with `nav`, at its box: the rectangles and images the recording backend gets, at DPR %i',
    async (dpr) => {
      const t = await opened(syntheticFamily(120, 8), { renderers: kernel });
      const ref = t.ctx.viewports.sections.value[0];
      const plain = parseXml(await renderSvg(t.ctx, { width: 801, scale: dpr, overlay: true }));
      const svg = parseXml(await renderSvg(t.ctx, { width: 801, scale: dpr, overlay: true, nav: true }));
      const classes = (e: ReturnType<typeof parseXml>) => all(e, 'g').map((g) => g.attrs.class);
      expect(classes(plain)).not.toContain('seqv-vnav');
      expect(classes(svg)).toEqual(expect.arrayContaining(['seqv-vnav', 'seqv-vnav-overlay']));
      // The same frame on the recording backend: the section at the export's size, the overview at its box.
      const vp = t.ctx.viewports.get(ref) as Viewport;
      vp.setSize(801, 1000, dpr);
      vp.setSize(801, Math.max(1, Math.ceil(t.ctx.layout.contentHeight(ref)?.value ?? 1000)), dpr);
      const { graph, selection, highlight, systems, layout } = t.ctx;
      const marks = new ChannelMarks({ graph, hover: new Hover(), selection, highlight, systems, layout });
      const vnav = new VerticalNav(t.ctx, ref, marks);
      const rec = recording();
      vnav.frame(rec, Dirty.all);
      const box = Array.from(vnav.box);
      expect(box[2]).toBeGreaterThan(0);
      // The offset group holding both passes, at the box's origin in CSS px.
      const at = all(svg, 'g').find((g) => g.children.some((c) => c.attrs?.class === 'seqv-vnav'));
      expect(at?.attrs.transform).toBe(`translate(${hund(box[0] / dpr)} ${hund(box[1] / dpr)})`);
      const want: string[] = [];
      const r = rec.rectItems;
      for (let i = 0; i < r.n; i++) {
        if (r.color[i] === -1) continue;
        // Outlines as the four edges the writer makes of them.
        const [x, y, w, h, s] = [r.x[i], r.y[i], r.w[i], r.h[i], r.stroke[i]];
        const cut =
          s === 0 || 2 * s >= w || 2 * s >= h
            ? [[x, y, w, h]]
            : [
                [x, y, w, s],
                [x, y + h - s, w, s],
                [x, y + s, s, h - 2 * s],
                [x + w - s, y + s, s, h - 2 * s],
              ];
        for (const [a, b, c, d] of cut)
          want.push([a, b, a + c, b + d].map((v, k) => hund((v + box[k % 2]) / dpr)).join(' '));
      }
      const got: string[] = [];
      for (const g of all(at as ReturnType<typeof parseXml>, 'g').filter((e) => e.attrs.class?.startsWith('seqv-vnav')))
        for (const p of all(g, 'path'))
          for (const [x, y, w, h] of rectsOf(p.attrs.d ?? ''))
            if (p.attrs['shape-rendering'])
              got.push([x, y, x + w, y + h].map((v, k) => hund(v + box[k % 2] / dpr)).join(' '));
      expect(want.length).toBeGreaterThan(8);
      expect(got.sort()).toEqual(want.sort());
      // The overview's image: one draw of it in each.
      const images = opNames(rec).filter((n) => n === 'image').length;
      expect(images).toBe(1);
      expect(
        all(at as ReturnType<typeof parseXml>, 'use').filter((u) => /^#i/.test(u.attrs['xlink:href'] ?? '')),
      ).toHaveLength(images);
      vnav.dispose();
      marks.dispose();
      t.ctx.dispose();
    },
  );
});

/** Every section's size, DPR, zoom, scroll and row scale. */
const views = (ctx: SeqViewerContext) =>
  ctx.viewports.sections.value.map((ref) => {
    const s = ctx.viewports.get(ref)?.live.current;
    return [ref, s?.width, s?.height, s?.dpr, s?.cell, s?.scrollX, s?.scrollY, s?.rowScale];
  });

describe('a context passed in', () => {
  it("is left as it was: zoom, scroll and row scale too, which the export's size would clamp", async () => {
    for (const doc of [case14(), allTypes()]) {
      // Zoomed in and scrolled to the far end and the bottom: a wider export clamps the scroll, a narrower one the zoom.
      const t = await opened(doc, { renderers: kernel });
      for (const ref of t.ctx.viewports.sections.value) {
        t.ctx.viewports.get(ref)?.setSize(1000, 120, 1);
        t.ctx.viewports.get(ref)?.zoomAt(3, 400);
        t.ctx.viewports.get(ref)?.scrollTo(1e9, 1e9);
        t.ctx.layout.get(ref);
      }
      const before = views(t.ctx);
      expect(before.some(([, , , , , x]) => (x as number) > 0)).toBe(true);
      for (const width of [1800, 300]) await renderSvg(t.ctx, { width, scale: 2, overlay: true, nav: true });
      expect(views(t.ctx)).toEqual(before);
      t.ctx.dispose();
    }
  });

  it('commits no viewport, so a host syncing from `committed` sees no export; a change pending before stays pending', async () => {
    const t = await opened(case14(), { renderers: kernel });
    const refs = t.ctx.viewports.sections.value;
    const commits: string[] = [];
    for (const ref of refs) {
      const vp = t.ctx.viewports.get(ref);
      vp?.setSize(1000, 120, 1);
      vp?.zoomAt(3, 400);
      vp?.scrollTo(1e9, 1e9);
      t.ctx.layout.get(ref);
      let replay = true;
      vp?.committed.subscribe((c) => {
        if (!replay) commits.push(`${ref} ${c.origin}`);
      });
      replay = false;
    }
    for (const width of [1800, 300]) await renderSvg(t.ctx, { width, overlay: true });
    expect(commits).toEqual([]);
    // A pointer's pan not settled yet: still pending after an export, and committed once, as it would have been.
    const vp = t.ctx.viewports.get(refs[0]);
    vp?.panBy(-40, 0, 'pointer');
    await renderSvg(t.ctx, { width: 1800 });
    expect(commits).toEqual([]);
    vp?.release();
    expect(commits).toEqual([`${refs[0]} pointer`]);
    // Under `quiet`, a view not put back commits once, at the end.
    const done = (vp as Viewport).quiet();
    vp?.scrollTo(0, 0, 'api');
    vp?.zoomAt(2, 0, 'api');
    expect(commits.length).toBe(1);
    done();
    expect(commits).toEqual([`${refs[0]} pointer`, `${refs[0]} api`]);
    t.ctx.dispose();
  });

  it("is resized and put back in one task: nothing that runs later sees the export's size", async () => {
    const t = await opened(case14(), { renderers: kernel });
    const vp = t.ctx.viewports.get(t.ctx.viewports.sections.value[0]);
    const seen: number[] = [];
    const later: Promise<void>[] = [];
    // What a frame or a ResizeObserver callback would read: a task queued on every change.
    vp?.live.version.subscribe(() => {
      later.push(new Promise((r) => setTimeout(() => r(void seen.push(vp.live.current.width)))));
    });
    await renderSvg(t.ctx, { width: 1500 });
    await Promise.all(later);
    expect(seen.length).toBeGreaterThan(0);
    expect(new Set(seen)).toEqual(new Set([1000]));
    t.ctx.dispose();
  });

  it("another view, a skin and a measurer: paused while they show, each put back, a host's skin override kept", async () => {
    const t = await opened(case1(), { renderers: kernel });
    const active: boolean[] = [];
    t.ctx.ui.active.subscribe((on) => void active.push(on));
    const port = t.ctx.text.port;
    const families = { ...t.ctx.text.families };
    const xml = await renderSvg(t.ctx, {
      width: 500,
      view: 'tracks',
      skin: { ...DarkSkin, fontMono: 'Menlo', fontUi: 'Inter' },
      ports: { measurer: metricsMeasurer(DefaultFontMetrics) },
    });
    expect(xml).toContain('font-family="Inter"');
    expect(active).toEqual([true, false, true]);
    expect([t.ctx.views.current.value, t.ctx.text.port, t.ctx.text.families]).toEqual(['text', port, families]);
    // No override was held: the skin is the reader's again, and the reader's next one shows.
    expect(t.ctx.skin.value).toBe(DefaultSkin);
    const read = { ...DefaultSkin, accent: '#123456' };
    t.ctx.publishSkin(read);
    expect(t.ctx.skin.value).toBe(read);
    // A host's override is held again after an export with a skin of its own; the reader's next one does not show.
    const host = { ...DefaultSkin, accent: '#654321' };
    t.ctx.setSkin(host);
    await renderSvg(t.ctx, { width: 500, skin: DarkSkin });
    expect(t.ctx.skin.value).toBe(host);
    t.ctx.publishSkin({ ...read });
    expect(t.ctx.skin.value).toBe(host);
    // The current view and no skin or measurer: nothing to pause for.
    await renderSvg(t.ctx, { width: 500 });
    expect(active).toEqual([true, false, true, false, true]);
    t.ctx.dispose();
  });
});

describe('a genome-sized document', () => {
  itWithin(
    '2 × 200,000 positions: zoomed out and at letters, bounded in time and size',
    async (within) => {
      const t = await opened(scaleDoc({ length: 200_000 }), { renderers: kernel });
      const ref = t.ctx.viewports.sections.value[0];
      for (const [name, zoom] of [['fit', null], ['letters', [100_000, 100_120]] as const] as const) {
        if (zoom) t.ctx.viewports.get(ref)?.zoomTo(0, zoom[0], zoom[1]);
        const xml = await within(3000, () => renderSvg(t.ctx, { width: 1200, overlay: true, nav: true }));
        // What the screen can show, not what the data holds: about 120 elements and 180 KB zoomed out, 230 and 53 KB at letters (measured).
        const elements = (xml.match(/<[a-z]/g) ?? []).length;
        expect(parseXml(xml).name).toBe('svg');
        expect(xml.length, `${name} ${xml.length} bytes`).toBeLessThan(600_000);
        expect(elements, `${name}: ${elements} elements`).toBeLessThan(5000);
      }
      t.ctx.dispose();
    },
    60_000,
  );
});
