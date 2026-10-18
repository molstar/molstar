/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tracks drawn in full in the text view (seqspec:viewer#the-text-view): a track with a `height` takes that many of the
// text's rows on every line, and its layers draw their track forms over each line: a values line (one path a line,
// broken across lines, one scale for every line), an area, glyphs, ranges as items with labels, secondary structure;
// without a height a track keeps its lane. Hits by drawn geometry; the axes and the labels (`wrappedTrackAxes`,
// `wrappedTrackLabels`), off by default, drawn as the SVG export draws them.

import { describe, expect, it } from 'vitest';
import { Dirty, Pass } from '../src/draw/index.js';
import { EntryKind, Hit, HoverKind, KernelRenderers, type LayerRenderer, type RenderProps } from '../src/index.js';
import { renderSvg } from '../src/svg/index.js';
import { recording } from './composer-fixtures.js';
import { drawn, layerOf, record, rects, texts } from './draw-fixtures.js';
import { opened } from './layout-fixtures.js';
import { lit, passRects } from './overlay-fixtures.js';
import { publicationDoc } from './prepare-fixtures.js';

const N = 120;
const NO_LABELS: Partial<RenderProps> = { wrappedGroupLabels: 'none', wrappedGutter: 0, wrappedTracks: 'full' };

describe('a track with a height in the text', () => {
  it("takes `height` of the text's rows on every line; without one, its lane", async () => {
    const t = await opened(publicationDoc(), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      size: { 'v.main': [600, 2000] },
      render: NO_LABELS,
    });
    const l = t.layout('v.main');
    const height = (ref: string) => t.entries('v.main').find((e) => e[1] === ref)?.[4];
    // 12 px letters: a text row is 20 px.
    expect([height('v.letters'), height('v.line'), height('v.area'), height('v.sites'), height('v.dom')]).toEqual([
      20, 40, 40, 20, 20,
    ]);
    expect(l.lines.count).toBeGreaterThan(1);
    const lane = await opened(publicationDoc({ lane: true }), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      size: { 'v.main': [600, 2000] },
      render: NO_LABELS,
    });
    expect(lane.entries('v.main').find((e) => e[1] === 'v.line')?.[4]).toBe(5);
  });

  it('is opt-in: under `wrappedTracks: lanes` (the default) the same document keeps its lanes whatever the height', async () => {
    const heights = async (render: Partial<RenderProps>) => {
      const t = await opened(publicationDoc(), {
        renderers: [...KernelRenderers] as LayerRenderer[],
        size: { 'v.main': [600, 2000] },
        render,
      });
      const first = (ref: string) => t.entries('v.main').find((e) => e[1] === ref)?.[4];
      return ['v.letters', 'v.line', 'v.area', 'v.sites', 'v.dom', 'v.ss'].map(first);
    };
    // Lanes: the letters, two strips (2 + 3), and underlines (2 + 5 a lane).
    expect(await heights({ wrappedGroupLabels: 'none', wrappedGutter: 0 })).toEqual([20, 5, 5, 7, 7, 7]);
    expect(await heights(NO_LABELS)).toEqual([20, 40, 40, 20, 20, 20]);
  });

  it('values: a line per line, broken across lines, on one scale for every line; an area', async () => {
    const t = await drawn(publicationDoc(), { width: 600, height: 2000, render: NO_LABELS });
    const lines = t.ctx.layout.get(layerOf(t.ctx, 'Hydrophobicity').section)?.current.lines.count ?? 0;
    const { rec } = record(t.ctx, 'Hydrophobicity');
    const p = rec.pathItems;
    // One subpath a line at least (a color change would add more; there is no theme here).
    expect(p.n).toBe(lines);
    // The same value is at the same height in its band on every line: positions 0 and 6·2π (sin equal) on two lines.
    const ys = new Map<number, number[]>();
    for (let q = 0; q < p.n; q++) {
      const c0 = p.cstart[q];
      const c1 = q + 1 < p.n ? p.cstart[q + 1] : p.nc;
      ys.set(q, [p.coords[c0 + 1], p.coords[c1 - 1]]);
    }
    const bands = t.entries(layerOf(t.ctx, 'Hydrophobicity').section).filter((e) => e[1] === 'v.line');
    const rel = (q: number) => (ys.get(q)?.[0] ?? 0) - (bands[q][3] + t.layout('v.main').scrolledTop);
    // Every line's path inside its band.
    for (let q = 0; q < p.n; q++) expect(rel(q)).toBeGreaterThanOrEqual(0);
    for (let q = 0; q < p.n; q++) expect(rel(q)).toBeLessThanOrEqual(40);
    const area = rects(record(t.ctx, 'Confidence').rec).filter((r) => r[5] === 77);
    expect(area.length).toBeGreaterThan(N / 2);
  });

  it("the line runs from each text line's first cell's left edge to its last's right, as the area does", async () => {
    const t = await drawn(publicationDoc(), { width: 600, height: 2000, render: NO_LABELS });
    const { rec } = record(t.ctx, 'Confidence');
    const p = rec.pathItems;
    const area = rects(rec).filter((r) => r[5] === 77);
    const l = t.layout('v.main');
    const bands = t.entries(layerOf(t.ctx, 'Confidence').section).filter((e) => e[1] === 'v.area');
    expect(p.n).toBe(l.lines.count);
    for (let q = 0; q < p.n; q++) {
      const c0 = p.cstart[q];
      const c1 = q + 1 < p.n ? p.cstart[q + 1] : p.nc;
      const top = l.scrolledTop + bands[q][3];
      const mine = area.filter((r) => r[1] >= top - 1 && r[1] + r[3] <= top + bands[q][4] + 1);
      expect(p.coords[c0], `line ${q} start`).toBe(Math.min(...mine.map((r) => r[0])));
      expect(p.coords[c1 - 2], `line ${q} end`).toBe(Math.max(...mine.map((r) => r[0] + r[2])));
    }
  });

  it('an area is continuous through the gaps between groups, as the line; bars stay a cell each', async () => {
    /** The x extents (device px) the encoding's rects cover on the first line's band, merged. */
    const covered = async (encoding: 'area' | 'bars') => {
      const doc = publicationDoc();
      type N = { kind: string; params?: Record<string, unknown>; children?: N[] };
      const walk = (n: N) => {
        if (n.kind === 'layer' && n.params?.encoding === 'area') n.params.encoding = encoding;
        n.children?.forEach(walk);
      };
      walk(doc.root as N);
      const t = await drawn(doc, { width: 600, height: 2000, render: NO_LABELS });
      const e = t.entries(layerOf(t.ctx, 'Confidence').section).find((x) => x[1] === 'v.area');
      const top = (e?.[3] ?? 0) + t.layout('v.main').scrolledTop;
      const r = rects(record(t.ctx, 'Confidence').rec)
        .filter((x) => x[1] + x[3] <= top + 41 && x[1] >= top - 1 && (encoding === 'bars' || x[5] === 77))
        .sort((a, b) => a[0] - b[0]);
      const runs: number[][] = [];
      for (const [x, , w] of r) {
        const last = runs[runs.length - 1];
        if (last && x <= last[1]) last[1] = Math.max(last[1], x + w);
        else runs.push([x, x + w]);
      }
      return runs;
    };
    // The first line holds several groups of 10: the area one run across them, the bars a run a cell.
    const area = await covered('area');
    expect(area.length).toBe(1);
    const bars = await covered('bars');
    expect(bars.length).toBeGreaterThan(20);
  });

  it('points as glyphs, ranges as items with their labels, secondary structure as shapes: not underlines', async () => {
    const t = await drawn(publicationDoc(), { width: 600, height: 2000, render: NO_LABELS });
    expect(record(t.ctx, 'Sites').rec.glyphItems.n).toBe(4);
    const dom = record(t.ctx, 'Domains').rec;
    expect(rects(dom).some((r) => r[3] > 6)).toBe(true);
    expect(texts(dom).some((d) => d.text.startsWith('Kinase') || d.text.startsWith('Linker'))).toBe(true);
    // A helix is a cylinder (a path), drawn in the track's band, not an underline.
    expect(record(t.ctx, 'Secondary').rec.pathItems.n).toBeGreaterThan(0);
  });

  it("a glyph is hit as its item; a hover's outline is drawn on the glyph's line", async () => {
    const t = await lit(publicationDoc(), { width: 600, height: 2000, render: NO_LABELS });
    const { rec, frame } = record(t.ctx, 'Sites');
    const g = rec.glyphItems;
    let hits = 0;
    for (let i = 0; i < g.n; i++) {
      const h = Hit.create();
      t.ctx.layout.hitTest(t.section, g.x[i] / frame.v.dpr, g.y[i] / frame.v.dpr, h);
      if (h.kind === HoverKind.item) hits++;
    }
    expect(hits).toBe(g.n);
    t.hover(g.x[0], g.y[0]);
    expect(passRects(t.overlay(), Pass.overlay).length + t.overlay().glyphItems.n).toBeGreaterThan(0);
  });

  it('axes and labels: off by default; on, the domain at each band and the label on the first line, in SVG too', async () => {
    const off = await lit(publicationDoc(), { width: 600, height: 2000, render: NO_LABELS });
    const all = (r: ReturnType<typeof recording>) => {
      const out: string[] = [];
      for (let i = 0; i < r.textItems.n; i++) out.push(r.textItems.text[i]);
      return out;
    };
    const plain = all(off.frame(recording(), Dirty.all));
    expect(plain).not.toContain('Hydrophobicity');
    expect(plain).not.toContain('100');
    const render = { ...NO_LABELS, wrappedTrackAxes: true, wrappedTrackLabels: true };
    const on = await lit(publicationDoc(), { width: 600, height: 2000, render });
    const drawnTexts = all(on.frame(recording(), Dirty.all));
    for (const label of ['Hydrophobicity', 'Confidence', 'Sites', 'Domains', 'Secondary'])
      expect(drawnTexts.filter((x) => x === label).length, label).toBe(1);
    const lines = on.layout(on.section).lines.count;
    // Each values track: its domain's ends on every line.
    expect(drawnTexts.filter((x) => x === '100').length).toBe(2 * lines);
    expect(drawnTexts.filter((x) => x === '0').length).toBe(2 * lines);
    // The axes' room: the lines start 30 px further in.
    expect(on.layout(on.section).dataLeft - off.layout(off.section).dataLeft).toBe(30);
    const xml = await renderSvg(publicationDoc(), { width: 600, render });
    expect(xml).toContain('Hydrophobicity');
    expect(xml).toContain('>100<');
    expect(xml).not.toMatch(/NaN|Infinity|undefined/);
    // Entries are the lines' bands: no labels band (group labels off).
    expect(on.layout(on.section).entryKind.subarray(0, on.layout(on.section).entryCount)).not.toContain(
      EntryKind.labels,
    );
  });
});

describe("a line's block in the text drawn in full", () => {
  /** The figure's spacing and labels, group labels above each line; `wrappedMaxWidth` 160: 2 groups (20 positions) a line. */
  const FIGURE: Partial<RenderProps> = {
    wrappedTracks: 'full',
    wrappedTrackLabels: true,
    wrappedTrackGap: 3,
    wrappedGap: 18,
    wrappedGroupLabels: 'aboveEnd',
    wrappedGutter: 0,
    wrappedMaxWidth: 160,
  };
  const open = (render: Partial<RenderProps>, doc = publicationDoc()) =>
    opened(doc, { renderers: [...KernelRenderers] as LayerRenderer[], size: { 'v.main': [600, 4000] }, render });
  /** Per line, its entries as [kind, track ref, top, height]. */
  const byLine = (t: Awaited<ReturnType<typeof open>>) => {
    const l = t.layout('v.main');
    const all = t.entries('v.main');
    const out: [number, string, number, number][][] = Array.from({ length: l.lines.count }, () => []);
    all.forEach((e, k) => out[l.entryRow[k]]?.push([e[0], e[1], e[3], e[4]]));
    return out;
  };

  it("a range's label on each line its item crosses, inside its part there, never past the line's end", async () => {
    const t = await drawn(publicationDoc(), { width: 600, height: 4000, render: FIGURE });
    const { rec, frame } = record(t.ctx, 'Domains');
    const l = t.layout('v.main');
    const v = frame.v;
    const names = texts(rec).filter((d) => /^(Kinase|Linker|Lin|Kin)/.test(d.text) || d.text.endsWith('…'));
    // "Linker region" over [60, 110): on lines [60, 80), [80, 100) and [100, 110), each its own label.
    expect(names.filter((d) => d.text.startsWith('Lin')).length).toBeGreaterThanOrEqual(2);
    const dom = t.trackIndex('v.dom');
    for (const d of names) {
      // The line whose band holds the label.
      let line = -1;
      for (let e = 0; e < l.entryCount; e++)
        if (l.entryTrack[e] === dom && l.entryKind[e] === EntryKind.track) {
          const y0 = (l.scrolledTop + l.entryTop[e]) * v.dpr;
          if (d.y >= y0 && d.y <= y0 + l.entryHeight[e] * v.dpr) line = l.entryRow[e];
        }
      expect(line).toBeGreaterThanOrEqual(0);
      v.onlyBlock = l.lines.block[line];
      v.lineStart = l.lines.start[line];
      v.lineStop = l.lines.stop[line];
      v.lineLeft = l.lines.left[line];
      const half = (t.ctx.text.label(d.text, d.font).width * v.dpr) / 2;
      expect(d.x - half, d.text).toBeGreaterThanOrEqual(v.xd(0, v.lineStart) - 1);
      expect(d.x + half, d.text).toBeLessThanOrEqual(v.xe(0, v.lineStop - 1) + 1);
    }
  });

  it("small gaps inside a line's block, a larger one between blocks; the group labels topmost", async () => {
    const t = await open(FIGURE);
    const l = t.layout('v.main');
    const lines = byLine(t);
    let inner = 0;
    for (const row of lines) {
      // Its group labels first.
      expect(row[0][0]).toBe(EntryKind.labels);
      for (let i = 1; i < row.length; i++) inner = Math.max(inner, row[i][2] - (row[i - 1][2] + row[i - 1][3]));
    }
    for (let k = 1; k < l.lines.count; k++) expect(l.lines.top[k] - l.lines.bottom[k - 1]).toBe(18);
    expect(inner).toBeLessThanOrEqual(3);
    expect(18).toBeGreaterThan(2 * inner);
  });

  it('a track of features with nothing on a line takes no band there; values keep theirs; a prop turns it off', async () => {
    const t = await open(FIGURE);
    const lines = byLine(t);
    // Line [40, 60): no domain (Kinase [2, 40), Linker region [60, 110)).
    const l = t.layout('v.main');
    const k = Array.from({ length: l.lines.count }, (_, i) => i).find((i) => l.lines.start[i] === 40) ?? -1;
    const height = (ref: string) => lines[k].find((e) => e[1] === ref)?.[3] ?? 0;
    expect([height('v.dom'), height('v.line'), height('v.area')]).toEqual([0, 40, 40]);
    const kept = await open({ ...FIGURE, wrappedCollapseEmpty: false });
    expect(byLine(kept)[k].find((e) => e[1] === 'v.dom')?.[3]).toBe(20);
  });

  it("a track's label in a row of its own above its band, on the first line it has one: never over its data", async () => {
    const t = await open(FIGURE);
    const all = t.entries('v.main');
    const rows = all.filter((e) => e[0] === EntryKind.trackLabel);
    expect(rows.map((e) => e[1])).toEqual(['v.line', 'v.area', 'v.sites', 'v.dom', 'v.ss']);
    for (const r of rows) {
      const band = all.find((e) => e[0] === EntryKind.track && e[1] === r[1] && e[4] > 0);
      expect(band && r[3] + r[4]).toBeLessThanOrEqual(band?.[3] ?? -1);
    }
  });

  it('a track ordered before the letters sits between the group labels and the letters', async () => {
    const doc = publicationDoc();
    type N = { kind: string; ref?: string; params?: Record<string, unknown>; children?: N[] };
    const root = doc.root as N;
    const data = root.children?.find((c) => c.kind === 'data') as N;
    data.children?.push({
      kind: 'annotation',
      ref: 'above',
      params: { geometry: 'ranges', target: 'a', value: { start: [3, 45, 85], stop: [15, 55, 95] } },
    });
    const section = root.children?.find((c) => c.kind === 'view')?.children?.find((c) => c.kind === 'section') as N;
    section.children?.unshift({
      kind: 'track',
      ref: 'v.above',
      params: {},
      children: [{ kind: 'layer', ref: 'v.above.layer', params: { type: 'ranges', data: 'above' } }],
    });
    const t = await open(FIGURE, doc);
    for (const row of byLine(t)) {
      const top = (pred: (e: [number, string, number, number]) => boolean) => row.find(pred)?.[2] ?? -1;
      const labels = top((e) => e[0] === EntryKind.labels);
      const above = top((e) => e[1] === 'v.above');
      const letters = top((e) => e[1] === 'v.letters');
      expect(labels).toBeLessThan(above);
      expect(above).toBeLessThan(letters);
    }
  });
});
