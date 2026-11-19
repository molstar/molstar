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

import { createSeqContext } from '@molstar/seqx-core/context';
import { DataStore, DefaultStoreSettings, SeqRuntime } from '@molstar/seqx-core/data';
import {
  attachSources,
  type BandsLay,
  bundledTextMetrics,
  LightTokens,
  NavigatorRenderer,
  Renderer,
  type RangesLay,
  RowsLay,
  type Viewport,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
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
  const model = resolveSpec(spec, { runtimes: [SeqRuntime, ...(runtimes as never[])] });
  const [view] = resolveViews(model);
  const store = new DataStore();
  expect(await attachSources(model, store, { assets: nodeFixtureAssets().reader })).toEqual([]);
  const common = { model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() };
  const renderer = new Renderer(common);
  const n = new NavigatorRenderer({ ...common, navigator: view.navigators[nav] });
  const draw = () => {
    const f = renderer.frame({ ...input, settled: true });
    const g = n.frame({ size: navSize, view: f.layout });
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
    expect(w.lo).toBeCloseTo(pic.axis!.x(RegionOnChr17), 6);
    expect(r.nav.issues.filter((i) => i.severity !== 'info')).toEqual([]);
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

  it('the Globins’ overview has every row, its window the rows in view', async () => {
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
    const r = await settled(navigators.plasmidStrip(), [], {
      size: { width: 800, height: 900 },
      viewports: { '/views/main/ring': { start: 2600, end: 2760 } },
    });
    expect(r.nav.windows[0].wrap).toBe(true);
  });

  it.each([
    ['gene-region-ideogram', 'chr17'],
    ['tracks-strip', 'strip'],
    ['protein-gene-strip', 'strip'],
    ['detail-strip', 'strip'],
    ['scale-navigators', 'overview'],
  ])('%s exports with its navigators, titled', async (id, nav) => {
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
  });
});
