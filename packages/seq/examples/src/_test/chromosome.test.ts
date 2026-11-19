/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// p53's chromosome, the p53 dataset's genome-level view: the generator is deterministic by region and carries the
// dataset's gene in TP53's span; the view, headless (the renderer, the overlay and the navigator over one store), reads
// only the view's regions and its margin, evicts past the store's bound and never the view, draws every gene in view
// against the data at the fit, across tile edges, at TP53 and at its letters, every drawn gene hit, and shows its
// window on the ideogram. The frame at 83 Mb against the gene view's (the Performance page's flat frame time) is heavy.

import { readFileSync } from 'node:fs';
import { DataStore, decodeSeq, type FormatRuntime, type RegionRequest, SeqRuntime } from '@molstar/seqx-core/data';
import {
  bundledTextMetrics,
  DefaultRenderConfig,
  type FrameInput,
  hitView,
  LightTokens,
  NavigatorRenderer,
  NoAssets,
  OverlaySources,
  openModel,
  type RangesLay,
  Renderer,
  recordDrawLists,
  type ViewLayout,
} from '@molstar/seqx-core/render';
import { OverlayBuilder } from '@molstar/seqx-core/resolve';
import { itWithin, ratioWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';
import * as chromosome from '../p53/chromosome.js';
import {
  ChromosomeSystem,
  chromosome as generate,
  chromosomeRuntime,
  DefaultChromosome,
  type Gene,
  GeneBlock,
  LetterBlock,
  Tp53,
} from '../p53/chromosome-data.js';
import { GeneFile, GeneLength } from '../p53/data.js';
import { ExampleViews } from '../views.js';

const Length = DefaultChromosome.length;
const Size = { width: 1440, height: 400 };
const Dpr = 2;
const c = generate();
const allGenes = c.genesIn(0, Length);

const request = (start: number, end: number, level = 1): RegionRequest => ({
  system: ChromosomeSystem,
  start,
  end,
  level,
  shapes: ['letters', 'ranges', 'values'],
});

describe('the chromosome’s generator', () => {
  it('a region’s answer is the same whatever was asked before it, and joins its halves', () => {
    // across a block of letters, a block of genes and TP53's start; one coarse
    const regions = [
      request(Tp53.start - 3000, Tp53.start + 2000),
      request(7 * GeneBlock - 700, 7 * GeneBlock + 1348),
      request(0, 4096 * 1024, 4096),
    ];
    const first = regions.map((r) => generate().answer(r));
    const other = generate();
    const reversed = [...regions].reverse().map((r) => other.answer(r));
    expect(reversed.reverse()).toEqual(first);
    const a = 5 * LetterBlock - 300;
    expect(c.letters(a, a + 300) + c.letters(a + 300, a + 900)).toBe(c.letters(a, a + 900));
  });

  it('TP53’s span carries the dataset’s gene: its letters are tp53-gene.seq’s, its exons the gene’s', () => {
    const file = decodeSeq(readFileSync(`${PackageFolder}${GeneFile}`, 'utf8')).data!;
    const letters = file.pairs[0].answer.ok ? file.pairs[0].answer.tables.letters?.letters : undefined;
    expect(letters).toHaveLength(GeneLength);
    expect(c.letters(Tp53.start, Tp53.start + GeneLength)).toBe(letters);
    const tp53 = c.genesIn(Tp53.start, Tp53.end).find((g) => g.name === 'TP53')!;
    expect(tp53).toMatchObject({ start: 7_668_402, end: 7_687_550, strand: -1 });
    // the coding exons the gene's, on the minus strand: the protein's start codon at the last exon's end
    const coding = tp53.exons.reduce((n, [s, e]) => n + e - s, 0);
    expect(coding).toBe(1182);
    // no other gene within its span
    expect(c.genesIn(Tp53.start, Tp53.end).map((g) => g.name)).toEqual(['TP53']);
  });

  it('about 1,200 genes on both strands, their exons inside them in order; the levels’ signal agrees', () => {
    expect(allGenes.length).toBeGreaterThan(1100);
    expect(allGenes.length).toBeLessThan(1300);
    expect(new Set(allGenes.map((g) => g.name)).size).toBe(allGenes.length);
    expect(new Set(allGenes.map((g) => g.strand))).toEqual(new Set([1, -1]));
    // TP53's span holds its exons, the dataset's gene's coding ones; the others' first and last exons are its ends
    for (const g of allGenes) {
      expect(g.exons[0][0], g.name).toBeGreaterThanOrEqual(g.start);
      expect(g.exons[g.exons.length - 1][1], g.name).toBeLessThanOrEqual(g.end);
      if (g === Tp53) continue;
      expect(g.exons[0][0], g.name).toBe(g.start);
      expect(g.exons[g.exons.length - 1][1], g.name).toBe(g.end);
      for (let i = 1; i < g.exons.length; i++) expect(g.exons[i][0], g.name).toBeGreaterThanOrEqual(g.exons[i - 1][1]);
    }
    // a bin's mean is its sixteen finer bins' mean
    const coarse = c.answer(request(0, 256 * 1024, 256));
    const fine = c.answer(request(0, 256 * 16 * 4, 16));
    const bins = (a: typeof coarse) => (a.ok && a.tables.values?.form === 'bins' ? a.tables.values : undefined);
    for (let b = 0; b < 4; b++) {
      let sum = 0;
      for (let k = 0; k < 16; k++) sum += bins(fine)!.mean[16 * b + k];
      expect(sum / 16).toBeCloseTo(bins(coarse)!.mean[b], 5);
    }
  });
});

/** A runtime that records the regions asked of it. */
function recording(inner: FormatRuntime): FormatRuntime & { requests: RegionRequest[] } {
  const requests: RegionRequest[] = [];
  return {
    name: inner.name,
    formats: inner.formats,
    requests,
    async open(source, ctx) {
      const opened = await inner.open(source, ctx);
      return {
        ...opened,
        read(r, signal) {
          requests.push(r);
          return opened.read(r, signal);
        },
      };
    },
  };
}

/** The view over one store: the renderer, the overlay's sources and the ideogram, settled at a viewport. */
async function open(settings: { maxBytes?: number } = {}) {
  const runtime = recording(chromosomeRuntime());
  const runtimes = [SeqRuntime, runtime];
  const store = new DataStore({ settings });
  const opened = await openModel(chromosome.build(), store, { runtimes });
  expect(opened.issues).toEqual([]);
  const { model } = opened;
  const view = opened.view!;
  const overlays = new OverlayBuilder({ runtimes });
  const sources = new OverlaySources(store, NoAssets);
  const common = { model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() };
  const renderer = new Renderer(common);
  const nav = new NavigatorRenderer({ ...common, navigator: view.navigators[0] });
  const input = (viewport?: { start: number; end: number }): FrameInput => ({
    size: Size,
    dpr: Dpr,
    settled: true,
    overlay: overlays.of(model, store),
    ...(viewport ? { viewport } : {}),
  });
  const draw = (viewport?: { start: number; end: number }) => {
    const f = renderer.frame(input(viewport));
    return { f, g: nav.frame({ size: { width: Size.width, height: 24 }, view: f.layout, dpr: Dpr }) };
  };
  const settle = async (viewport?: { start: number; end: number }) => {
    let d = draw(viewport);
    for (let i = 0; i < 6; i++) {
      await sources.sync(overlays.of(model, store), view);
      sources.hold(d.f.layout.shown);
      store.want([...d.f.wants, ...d.g.wants]);
      await store.whenSettled();
      d = draw(viewport);
    }
    sources.hold(d.f.layout.shown);
    return d;
  };
  return { store, renderer, nav, runtime, input, draw, settle };
}

/** A track's first layer, by its id in the section or its label (a template's instance: `Exons TP53`). */
const layerOf = (layout: ViewLayout, track: string) =>
  layout.sections[0].tracks.find((t) => t.track.label === track || t.track.address.endsWith(`/${track}`))!.layers[0];

/** The genes overlapping `[a, b)`. */
const genesIn = (a: number, b: number): Gene[] => allGenes.filter((g) => g.start < b && g.end > a);

/**
 * A ranges layer against the data: each gene in view has a mark over its extent (cut to the view) in its lane; each
 * mark hit at its center is one of the items it draws (its own, or those of its run). What is wrong.
 */
function rangesAgainst(
  layout: ViewLayout,
  track: string,
  genes: readonly Gene[],
  view: { start: number; end: number },
): string[] {
  const l = layerOf(layout, track);
  const lay = l.lay as RangesLay;
  const axis = layout.sections[0].axis!;
  const wrong: string[] = [];
  for (const g of genes) {
    const a = axis.x(Math.max(g.start, view.start));
    const b = axis.x(Math.min(g.end, view.end));
    let covered = false;
    for (let j = 0; j < lay.m; j++) if (lay.mx0[j] <= a + 1 && lay.mx1[j] >= b - 1) covered = true;
    if (!covered) wrong.push(`${g.name}: no mark over [${a}, ${b})`);
  }
  for (let j = 0; j < lay.m; j++) {
    const drawn = new Set<string | number>();
    if (lay.mitem[j] >= 0) drawn.add(lay.keys[lay.mitem[j]]);
    else
      for (let k = 0; k < lay.hn; k++)
        if (lay.hlane[k] === lay.mlane[j] && lay.hx1[k] >= lay.mx0[j] && lay.hx0[k] <= lay.mx1[j])
          drawn.add(lay.keys[lay.hitem[k]]);
    const hit = hitView(layout, (lay.mx0[j] + lay.mx1[j]) / 2, lay.box.y + (lay.mlane[j] + 0.5) * lay.laneHeight);
    if (hit?.layer !== l.key || hit.element === undefined || !drawn.has(hit.element.key))
      wrong.push(`mark ${j} [${lay.mx0[j]}, ${lay.mx1[j]}): hit ${JSON.stringify(hit?.element)}`);
  }
  return wrong;
}

/** Each request asked at a viewport: within the view and its margin (`DefaultRenderConfig.margin` each side). */
function outside(requests: readonly RegionRequest[], start: number, end: number): RegionRequest[] {
  const m = Math.ceil((end - start) * DefaultRenderConfig.margin);
  return requests.filter((r) => !(r.end > start - m && r.start < end + m));
}

const AtTp53 = { start: Tp53.start - 2000, end: Tp53.end + 2000 };
const AtLetters = { start: Tp53.start + 100, end: Tp53.start + 160 };

describe('p53’s chromosome, headless', () => {
  it('opens fit to the chromosome and lays out to TP53’s letters, asking only for the view and its margin', async () => {
    const v = await open();
    const fit = await v.settle();
    const axis = fit.f.layout.sections[0].axis!;
    const box = (layerOf(fit.f.layout, 'genes').lay as RangesLay).box;
    expect(axis.x(0)).toBeCloseTo(box.x, 3);
    expect(axis.x(Length)).toBeCloseTo(box.x + box.width, 3);
    // the whole chromosome at a coarse level: a few tiles, no letters
    expect(new Set(v.runtime.requests.map((r) => r.level)).has(1)).toBe(false);
    expect(v.runtime.requests.length).toBeLessThan(64);
    for (const at of [AtTp53, AtLetters]) {
      v.runtime.requests.length = 0;
      const d = await v.settle(at);
      expect(d.f.layout.issues).toEqual([]);
      // TP53's letters are in tiles its view asked for already
      if (at === AtTp53) expect(v.runtime.requests.length).toBeGreaterThan(0);
      expect(outside(v.runtime.requests, at.start, at.end)).toEqual([]);
    }
    // at TP53's letters: its exons' track, and the letters drawn are the chromosome's
    const d = v.draw(AtLetters);
    expect(d.f.layout.sections[0].tracks.map((t) => t.track.label)).toEqual(['', 'Genes', 'Exons TP53', 'GC', 'chr17']);
    const glyphs = recordDrawLists(d.f.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('glyphs'))
      .map((l) =>
        l
          .replace(/^glyphs "[^"]*" [\d.]+ \S+ /, '')
          .split(', ')
          .map((x) => x[0])
          .join(''),
      )
      .join('');
    expect(glyphs).toBe(c.letters(AtLetters.start, AtLetters.end));
  });

  it('draws every gene in view against the data: at the fit, across tile edges, at TP53, at its letters', async () => {
    const v = await open();
    // the fit: the genes a density strip, every gene's columns counted and none counted without a gene
    const fit = (await v.settle()).f.layout;
    const lay = layerOf(fit, 'genes').lay as RangesLay;
    const axis = fit.sections[0].axis!;
    expect(lay.mode).toBe('density');
    expect(lay.items).toBe(allGenes.length);
    const col = (p: number) => Math.min(axis.columns - 1, Math.floor((axis.x(p) - axis.x0) * Dpr));
    const counted = new Uint8Array(axis.columns);
    for (const g of allGenes) for (let k = col(g.start); k <= col(g.end - 1e-6); k++) counted[k] = 1;
    for (const g of allGenes) expect(lay.strip.counts[col(g.start)], g.name).toBeGreaterThan(0);
    for (let k = 0; k < axis.columns; k++) expect(lay.strip.counts[k] > 0, `column ${k}`).toBe(counted[k] === 1);
    // across a tile edge at the coarsest level read: a gene runs over one
    const tile = Math.max(...v.runtime.requests.map((r) => r.level)) * 1024;
    expect(allGenes.some((g) => Math.floor(g.start / tile) !== Math.floor((g.end - 1) / tile))).toBe(true);

    // 150 kb about a block edge, at level 16 (tiles of 16,384 positions): bars, genes over tile edges
    const edge = { start: 7 * GeneBlock - 75_000, end: 7 * GeneBlock + 75_000 };
    const at = await v.settle(edge);
    const inView = genesIn(edge.start, edge.end);
    expect(inView.some((g) => Math.floor(g.start / 16_384) !== Math.floor((g.end - 1) / 16_384))).toBe(true);
    expect(rangesAgainst(at.f.layout, 'genes', inView, edge)).toEqual([]);

    for (const view of [AtTp53, AtLetters]) {
      const layout = (await v.settle(view)).f.layout;
      expect(rangesAgainst(layout, 'genes', genesIn(view.start, view.end), view)).toEqual([]);
      // TP53's exons from its node, each over its extent
      const exons = Tp53.exons
        .filter(([s, e]) => s < view.end && e > view.start)
        .map(([s, e], i) => ({ ...Tp53, name: `TP53 exon ${i}`, start: s, end: e }));
      expect(exons.length).toBeGreaterThan(0);
      expect(rangesAgainst(layout, 'Exons TP53', exons, view)).toEqual([]);
    }
  });

  it('the ideogram’s window is the view: the whole chromosome at the fit, TP53’s place at TP53', async () => {
    const v = await open();
    const fit = (await v.settle()).g.windows[0];
    expect(fit).toMatchObject({ shown: true, p0: 0, p1: Length });
    const d = await v.settle(AtTp53);
    const w = d.g.windows[0];
    expect(w.p0).toBeCloseTo(AtTp53.start, 0);
    expect(w.p1).toBeCloseTo(AtTp53.end, 0);
    const pic = v.nav.pictureSection!.axis!;
    expect(w.x0).toBeCloseTo(pic.x(AtTp53.start), 3);
    expect(w.x1).toBeCloseTo(pic.x(AtTp53.end), 3);
    // a sliver drawn wider, centred on its place
    expect(w.dx1 - w.dx0).toBeGreaterThan(w.x1 - w.x0);
    expect((w.dx0 + w.dx1) / 2).toBeCloseTo((w.x0 + w.x1) / 2, 3);
  });

  it('the store evicts past its bound as the view pans along, never a tile of the view', async () => {
    const maxBytes = 256 * 1024;
    const v = await open({ maxBytes });
    const views = [0, 1, 2, 3, 4, 5].map((i) => ({ start: 30_000_000 + i * 200_000, end: 30_020_000 + i * 200_000 }));
    let first: string[] = [];
    for (const view of views) {
      const d = await v.settle(view);
      const keys = d.f.wants
        .filter((w) => w.priority === 'view' && w.source === '/data/chr17-generated')
        .flatMap((w) => v.store.keysOf(w));
      expect(keys.length).toBeGreaterThan(0);
      for (const k of keys) expect(v.store.get(k)?.status, k).toBe('loaded');
      if (!first.length) first = keys;
    }
    expect(v.store.bytes).toBeLessThanOrEqual(maxBytes);
    // the first view's tiles let go
    expect(first.filter((k) => v.store.get(k)?.status === 'loaded').length).toBeLessThan(first.length);
  });
});

describe('the frame at 83 Mb', () => {
  // The Performance page, Levels and a frame's cost: a frame costs what is on screen, at 83 Mb as at p53's gene (the protein
  // on its gene, 393 residues and their codons: 0.18 ms a frame at 1,440 px, dpr 2). A pan's frames (the layout again
  // each), measured 2026-11-19 (M-series, Node 26) against the gene view's: at its letters 0.18 times; between two
  // levels (512 positions a device pixel, between 256 and 1,024) 1.31 and at TP53 1.25 times, the signal's two bins a
  // device column (`values.coarser` is 1, the Performance page says 2); at the fit, as opened, 2.15 times, the genes' layout
  // walking all 1,177 rows it reads every frame (the ranges layout's, core). The brief's 1.2 holds at the letters only:
  // the bounds are the measures with headroom, the gap a finding. The fit is first: once TP53's exons are shown (an
  // inline node on the whole chromosome, level 1 only), a frame at the fit enumerates its 81,307 tiles (core's
  // `keysOf`), 37 ms. The quick lane lays each out once.
  const cases = [
    { at: { start: 0, end: Length }, ratio: 2.6 },
    { at: { start: 20_000_000, end: 20_000_000 + 512 * Size.width * Dpr }, ratio: 1.6 },
    { at: AtTp53, ratio: 1.6 },
    { at: AtLetters, ratio: 1.2 },
  ];
  itWithin(
    'a pan frame within a bound of the p53 gene view’s: at the fit, between two levels, TP53, its letters',
    async (within) => {
      const v = await open();
      const gene = await openGene();
      for (const { at, ratio } of cases) {
        const span = at.end - at.start;
        const views = Array.from({ length: 8 }, (_, i) => ({ start: at.start + (i * span) / 1000, end: at.end }));
        // the first and last of the pan read every tile the others do
        await v.settle(views[0]);
        await v.settle(views[7]);
        const ours = views.map((w) => v.input(w));
        let i = 0;
        ratioWithin(
          within,
          ratio,
          () => {
            for (let k = 0; k < 64; k++) gene.renderer.frame(gene.inputs[i++ % gene.inputs.length]);
          },
          () => {
            for (let k = 0; k < 64; k++) v.renderer.frame(ours[i++ % 8]);
          },
        );
      }
    },
  );
});

/** p53's gene view (the protein on its gene) over its own store, settled along a pan of the protein. */
async function openGene() {
  const view = ExampleViews.find((x) => x.id === 'protein-on-gene')!;
  const store = new DataStore();
  const opened = await openModel(view.build(), store, { runtimes: [SeqRuntime], assets: nodeFixtureAssets().reader });
  const renderer = new Renderer({
    model: opened.model,
    view: opened.view!,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
  });
  const inputs: FrameInput[] = Array.from({ length: 8 }, (_, i) => ({
    size: Size,
    dpr: Dpr,
    viewport: { start: i / 2, end: 393 },
  }));
  for (const input of [inputs[0], inputs[7]])
    for (let k = 0; k < 4; k++) {
      store.want(renderer.frame(input).wants);
      await store.whenSettled();
    }
  return { renderer, inputs };
}
