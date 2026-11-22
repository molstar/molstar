/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The examples' navigators drawn against the data: chr17's every band has its mark over its extent and its hit at the
// fit, one across a tile edge too, and the region's window is its place on the chromosome; the tracks' strip draws
// every domain over its extent; the protein's strip draws its gene follower's window; the Globins' overview has every
// row and draws them as the shared reduction's image over its column; the plasmid's strip wraps the origin. Each
// exported with `navigators: include`, titled, settled, with no issue.

import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { BuiltInNavigatorPictures } from '@molstar/seqx-core/render/navigator-picture';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore, DefaultStoreSettings, SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import {
  type BandsLay,
  bundledTextMetrics,
  DarkTokens,
  DefaultRenderConfig,
  LightTokens,
  Renderer,
  type RangesLay,
  RowsLay,
  type Viewport,
  openModel,
} from '@molstar/seqx-core/render';
import { NavigatorRenderer } from '@molstar/seqx-core/render/navigator';
import { describe, expect, it } from 'vitest';
import { generatedRegistry } from '../globins/data.js';
import { Chr17Bands, Chr17BandsLevel, RegionOnChr17 } from '../navigators.js';
import * as navigators from '../navigators.js';
import { nodeFixtureAssets } from '../node-assets.js';
import { registry as regionRegistry } from '../p53/gene-by-region.js';
import { ExampleViews } from '../views.js';

/** A view and its navigator settled over one store (the view's wants and the navigator's answered). */
async function settled(
  spec: ReturnType<typeof navigators.tracksStrip>,
  runtimes: readonly unknown[],
  input: { size: { width: number; height: number }; viewports?: Record<string, Viewport> },
  nav = 0,
  navSize = { width: input.size.width, height: 48 },
) {
  const store = new DataStore();
  const opened = await openModel(spec, store, {
    runtimes: [SeqRuntime, ...(runtimes as never[])],
    assets: nodeFixtureAssets().reader,
  });
  expect(opened.issues).toEqual([]);
  const { model } = opened;
  const view = opened.view!;
  const common = {
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    pictures: BuiltInNavigatorPictures,
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
  };
  const renderer = new Renderer(common);
  const n = new NavigatorRenderer({ ...common, navigator: view.navigators[nav] });
  const draw = () => {
    const f = renderer.frame({ ...input, settled: true });
    // settled as the viewer gives it (the store's): a still picture builds its level whole
    const g = n.frame({ size: navSize, view: f.layout, settled: true });
    return { f, g };
  };
  for (let i = 0; i < 8; i++) {
    const { f, g } = draw();
    store.want([...f.wants, ...g.wants]);
    await store.whenSettled();
  }
  const { f, g } = draw();
  return { view: f, nav: g, renderer: n, model };
}

/** Each item's marks (`mitem`) span its extent on the axis, within a px: none lost, none cut short. */
function marksOverExtents(
  lay: { m: number; mx0: Float32Array; mx1: Float32Array; mitem: Int32Array },
  extents: { starts: number[]; ends: number[] },
  x: (p: number) => number,
  items: number,
): string[] {
  const wrong: string[] = [];
  for (let i = 0; i < items; i++) {
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    for (let j = 0; j < lay.m; j++)
      if (lay.mitem[j] === i) {
        lo = Math.min(lo, lay.mx0[j]);
        hi = Math.max(hi, lay.mx1[j]);
      }
    const a = Math.min(x(extents.starts[i]), x(extents.ends[i]));
    const b = Math.max(x(extents.starts[i]), x(extents.ends[i]));
    if (!(Math.abs(lo - a) <= 1 && Math.abs(hi - b) <= 1)) wrong.push(`${i}: [${lo}, ${hi}) for [${a}, ${b})`);
  }
  return wrong;
}

describe('the examples’ navigators', () => {
  it('chr17’s ideogram: every band its mark and its hit; the region’s window its place on the chromosome', async () => {
    const r = await settled(
      navigators.geneRegionIdeogram(),
      regionRegistry({ delayMs: 0 }).map((e) => (e as { runtime: unknown }).runtime),
      { size: { width: 800, height: 300 } },
      0,
      { width: 800, height: 24 },
    );
    const pic = r.renderer.pictureSection!;
    const lay = pic.tracks[0].layers[0].lay as BandsLay;
    expect(lay.rows.n).toBe(Chr17Bands.start.length);
    // each band drawn over its extent, those across a tile edge of the coarse level whole
    const tile = DefaultStoreSettings.tileBins * Chr17BandsLevel;
    const across = Chr17Bands.start.filter(
      (s, i) => Math.floor(s / tile) !== Math.floor((Chr17Bands.end[i] - 1) / tile),
    );
    expect(across.length).toBeGreaterThan(0);
    expect(marksOverExtents(lay, lay.rows, (p) => pic.axis!.x(p), lay.rows.n)).toEqual([]);
    for (let i = 0; i < Chr17Bands.start.length; i++) {
      const x = (pic.axis!.x(Chr17Bands.start[i]) + pic.axis!.x(Chr17Bands.end[i] - 1)) / 2;
      expect(r.renderer.hit(x, 12)?.band, Chr17Bands.id[i]).toBe(Chr17Bands.id[i]);
    }
    const w = r.nav.windows[0];
    expect(w.shown).toBe(true);
    expect(w.p0).toBeCloseTo(RegionOnChr17, 0);
    // the region is a sliver of chr17: its window on the zoom strip, whose axis is the region alone
    const zoom = r.renderer.pictureLayout!.sections[1]!;
    expect(w.lo).toBeCloseTo(zoom.axis!.x(RegionOnChr17), 6);
    expect(r.nav.issues.filter((i) => i.severity !== 'info')).toEqual([]);
  });

  it('chr17’s ideogram has a zoom strip: the region’s bands, its mark and the connectors, drawn and exported, both themes', async () => {
    const r = await settled(
      navigators.geneRegionIdeogram(),
      regionRegistry({ delayMs: 0 }).map((e) => (e as { runtime: unknown }).runtime),
      { size: { width: 800, height: 300 } },
      0,
      { width: 800, height: 24 },
    );
    expect(r.renderer.zoomHeight).toBe(16);
    const ctx = r.renderer.pictureSection!;
    const zoom = r.renderer.pictureLayout!.sections[1]!;
    // the zoom strip under the ideogram, within the 16 px it adds; its axis the region (q's p13.1 band there)
    expect(ctx.y + ctx.height).toBeLessThanOrEqual(zoom.y);
    expect(zoom.y + zoom.height).toBeCloseTo(24 + 16, 6);
    expect(zoom.axis!.positionAt(zoom.axis!.x0)).toBeGreaterThanOrEqual(RegionOnChr17);
    const bands = zoom.tracks[0].layers[0].lay as BandsLay;
    expect(bands.rows.n).toBeGreaterThan(0);
    expect(r.renderer.hit(400, zoom.y + 2)?.band).toBe('p13.1');
    // both strips in the data pass; the mark and the connectors in the window pass
    const data = r.nav.lists.passes.data.rects;
    const ys = Array.from(data.y.subarray(0, data.n));
    expect(ys.some((y) => y < ctx.y + ctx.height)).toBe(true);
    expect(ys.some((y) => y >= zoom.y)).toBe(true);
    const win = r.nav.lists.passes.window;
    expect(win.paths.n).toBe(2);
    const mark = Array.from({ length: win.rects.n }, (_, i) => i).find(
      (i) => win.rects.y[i] === ctx.y && win.rects.h[i] === ctx.height,
    );
    expect(mark).toBeDefined();
    expect(win.rects.w[mark!]).toBeGreaterThanOrEqual(16 - 1e-6);
    for (const tokens of [LightTokens, DarkTokens]) {
      r.renderer.setTokens(tokens);
      const g = r.renderer.frame({ size: { width: 800, height: 24 }, view: r.view.layout });
      const ink = g.lists.passes.window.paths;
      // the connectors at the quiet alpha of the window's ink, one weight with every line of the part
      const quiet =
        ((tokens.colors['seq-window'] & 0xffffff00) |
          Math.round(DefaultRenderConfig.navigator.window.idleAlpha * 255)) >>>
        0;
      for (let i = 0; i < 2; i++) expect(ink.stroke[i]).toBe(quiet);
      const view = ExampleViews.find((v) => v.id === 'gene-region-ideogram')!;
      const c = createSeqContext({
        registry: [{ kind: 'format-runtime', runtime: SeqRuntime }, ...(view.registry?.({ delayMs: 0 }) ?? [])],
        ports: { assets: () => nodeFixtureAssets() },
      });
      const svg = (
        await c.renderSvg(view.build(), { size: { width: 720, height: 320 }, navigators: 'include', tokens })
      ).svg;
      c.retire();
      // the view's 320, the ideogram's 24, its zoom strip's 16
      expect(svg).toMatch(/^<svg xmlns="http:\/\/www.w3.org\/2000\/svg" width="720" height="360"/);
      const hex = `#${((tokens.colors['seq-window'] >>> 8) & 0xffffff).toString(16).padStart(6, '0')}`;
      const nav = svg.slice(svg.indexOf('data-navigator="chr17"'));
      expect(nav.match(new RegExp(`<path d="M [^"]+ L [^"]+" fill="none" stroke="${hex}"`, 'g'))?.length).toBe(2);
    }
  });

  it.each([
    ['tracks-strip', 'p53’s protein strip'],
    ['protein-gene-strip', 'the protein’s strip with its gene'],
    ['detail-strip', 'the plasmid’s strip'],
  ])('%s keeps one strip: its reachable span is wider than the floor (%s)', async (id) => {
    const view = ExampleViews.find((v) => v.id === id)!;
    const r = await settled(
      view.build(),
      (view.registry?.({ delayMs: 0 }) ?? []).map((e) => (e as { runtime: unknown }).runtime),
      { size: { width: 800, height: 600 } },
    );
    expect(r.renderer.zoomHeight).toBe(0);
    expect(r.renderer.pictureLayout!.sections.length).toBe(1);
    expect(r.nav.lists.passes.window.paths.n).toBe(0);
  });

  it('the tracks’ strip draws every domain small', async () => {
    const r = await settled(navigators.tracksStrip(), [], { size: { width: 800, height: 300 } });
    const pic = r.renderer.pictureSection!;
    const domains = pic.tracks.find((t) => t.track.address.endsWith('/domains'))!.layers[0].lay as RangesLay;
    expect(domains.rows.n).toBeGreaterThan(0);
    // each domain small over its extent on the picture's axis
    expect(domains.overflow).toBe(0);
    expect(marksOverExtents(domains, domains, (p) => pic.axis!.x(p), domains.items)).toEqual([]);
    const viewDomains = r.view.layout.sections[0].tracks.find((t) => t.track.address.endsWith('/domains'))!.layers[0]
      .lay as { rows: { n: number } };
    expect(domains.rows.n).toBe(viewDomains.rows.n);
  });

  it('the protein’s strip draws its gene follower’s window beside its own', async () => {
    const r = await settled(navigators.proteinGeneStrip(), [], {
      size: { width: 800, height: 300 },
      viewports: { '/views/main/protein': { start: 100, end: 200 } },
    });
    const [p, g] = r.nav.windows;
    expect(p.shown && g.shown).toBe(true);
    expect(Math.abs(g.x0 - p.x0)).toBeLessThan(800 / 393 + 1e-6);
  });

  it('the Globins’ overview has every row, its window the rows in view', { tags: ['medium'] }, async () => {
    const r = await settled(
      navigators.globinsNavigators(),
      generatedRegistry({ delayMs: 0 }).map((e) => (e as { runtime: unknown }).runtime),
      { size: { width: 800, height: 400 } },
      0,
      { width: 96, height: 400 },
    );
    const lay = r.renderer.pictureSection!.tracks[0].layers[0].lay as RowsLay;
    expect(lay.display.count).toBe(10_000);
    const w = r.nav.windows[0];
    expect(w.rows).toBe(10_000);
    expect(w.y1).toBeGreaterThan(w.y0);
    // the rows summarized: the shared reduction's image over the column, every row in it
    const images = r.nav.lists.passes.data.images;
    expect(images.n).toBeGreaterThan(0);
    let top = Number.POSITIVE_INFINITY;
    let bottom = Number.NEGATIVE_INFINITY;
    for (let i = 0; i < images.n; i++) {
      top = Math.min(top, images.y[i]);
      bottom = Math.max(bottom, images.y[i] + images.h[i]);
      expect(images.x[i]).toBeGreaterThanOrEqual(lay.box.x - 1e-3);
      expect(images.x[i] + images.w[i]).toBeLessThanOrEqual(lay.box.x + lay.box.width + 1e-3);
    }
    expect(top).toBeCloseTo(w.ly, 0);
    expect(bottom).toBeCloseTo(w.hy, 0);
  });

  it('the plasmid’s strip wraps the origin with the detail’s window', async () => {
    const r = await settled(navigators.plasmidStrip(), [RecordsRuntime], {
      size: { width: 800, height: 900 },
      viewports: { '/views/main/ring': { start: 2600, end: 2760 } },
    });
    expect(r.nav.windows[0].wrap).toBe(true);
  });

  const exportsWithNavigators = async (id: string, nav: string) => {
    const view = ExampleViews.find((v) => v.id === id)!;
    const ctx = createSeqContext({
      registry: [{ kind: 'format-runtime', runtime: SeqRuntime }, ...(view.registry?.({ delayMs: 0 }) ?? [])],
      ports: { assets: () => nodeFixtureAssets() },
    });
    const r = await ctx.renderSvg(view.build(), { size: { width: 720, height: 320 }, navigators: 'include' });
    expect(r.settled).toBe(true);
    expect(r.issues.filter((i) => i.severity !== 'info')).toEqual([]);
    expect(r.svg).toContain(`data-navigator="${nav}"`);
    expect(r.svg).toContain('<title>Window: ');
    ctx.retire();
  };

  it.each([
    ['gene-region-ideogram', 'chr17'],
    ['tracks-strip', 'strip'],
    ['protein-gene-strip', 'strip'],
    ['detail-strip', 'strip'],
  ])('%s exports with its navigators, titled', exportsWithNavigators);

  // 10,000 rows: the medium lane's, the others above its quick representatives
  it('scale-navigators exports with its navigators, titled', { tags: ['medium'] }, () =>
    exportsWithNavigators('scale-navigators', 'overview'),
  );
});
