/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import type { GridLay, MatrixLay, RangesLay } from '@molstar/seqx-core/render';
import {
  DarkTokens,
  gridCellHover,
  LightTokens,
  matrixCellHover,
  recordDrawLists,
  spanText,
} from '@molstar/seqx-core/render';
import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { GenomesViews } from '@molstar/seqx-examples';
import { AliasedMarks } from '@molstar/seqx-examples/genomes/by-region';
import { ASection, BSection } from '@molstar/seqx-examples/genomes/compare';
import { DotLevels, Inversion } from '@molstar/seqx-examples/genomes/generate/alignment';
import { NGap } from '@molstar/seqx-examples/genomes/generate/assemblies';
import { Hotspot } from '@molstar/seqx-examples/genomes/generate/contacts';
import { Duplicated } from '@molstar/seqx-examples/genomes/generate/liftover';
import { DeletionCarriers, DeletionExtent, SampleCount } from '@molstar/seqx-examples/genomes/generate/variants';
import { findGene, GenomeSection, geneRange, OpenedGene } from '@molstar/seqx-examples/genomes/tracks';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { FormatRuntime, OpenedSource, RegionRequest } from '@molstar/seqx-core/data';
import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import { withViewRegistry } from '../view-registry.js';
import { checkMarks } from './draw-checks.js';

// Examples/Genomes, headless on a fake clock as the stories compose them (a viewer of each view, the dataset's files
// read as committed, each answer after the twin's delay; the SVG export in both themes): each view's cells of the
// interaction matrix (examples.md) and the edge cases of its row, each named in a title.

let viewer: SeqViewer | undefined;

// the genome and bridge chunks landed first: their imports are real I/O, never on the fake clock
beforeAll(async () => {
  await loadGenome();
  await loadBridges();
  await import('@molstar/seqx-viewer/bridge-view');
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const view = (id: string) => GenomesViews.find((v) => v.id === id)!;

const DelayMs = 150;
const raf = DomHost.requestAnimationFrame;
/** The region reads in flight (a settled view has none), and how many began. */
const inflight = new Set<Promise<unknown>>();
const began = { n: 0 };
/** The levels the reads asked for (EC-170). */
const levelsRead: number[] = [];

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  inflight.clear();
  levelsRead.length = 0;
});

/** The view's registry with each runtime's opens and reads kept: a test waits on them, never on a count of turns. */
function watched(entries: readonly SeqRegistryEntry[], opens: Promise<unknown>[]): SeqRegistryEntry[] {
  return entries.map((e) =>
    e.kind === 'format-runtime'
      ? {
          ...e,
          runtime: {
            ...e.runtime,
            open: (source, ctx) => {
              const p = e.runtime.open(source, ctx).then((o) =>
                Object.assign(Object.create(o) as OpenedSource, {
                  read: (request: RegionRequest, signal: AbortSignal) => {
                    const r = o.read(request, signal);
                    began.n++;
                    levelsRead.push(request.level ?? 1);
                    inflight.add(r);
                    const done = () => void inflight.delete(r);
                    r.then(done, done);
                    return r;
                  },
                }),
              );
              opens.push(p);
              return p;
            },
          } satisfies FormatRuntime,
        }
      : e,
  );
}

/** The fixtures read from the examples' folder (a disk read is not on the fake clock: answered at once). */
const reader = {
  readText: async (uri: string) => readFileSync(`${PackageFolder}${uri}`, 'utf8'),
  readBytes: async (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`)),
};

/** A view composed as the story composes it, loaded and settled (its opens awaited, then its reads). */
async function open(id: string) {
  const base = DefaultSeqViewerSpec();
  const opens: Promise<unknown>[] = [];
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: withViewRegistry(base.registry ?? [], watched(view(id).registry?.({ delayMs: DelayMs }) ?? [], opens)),
    ports: { assets: () => ({ reader }) },
  });
  viewer = v;
  v.resize(900);
  void v.load(view(id).build()).then(() => v.settle());
  await vi.advanceTimersByTimeAsync(200);
  for (let seen = -1; seen !== opens.length; ) {
    seen = opens.length;
    await Promise.allSettled(opens);
    await vi.advanceTimersByTimeAsync(200);
  }
  await settle(v);
  return v;
}

/** The view's reads done: the fake clock run in steps, the reads that began awaited, until two steps begin none. */
async function settle(v: SeqViewer) {
  for (let quiet = 0; quiet < 2; ) {
    const before = began.n;
    await vi.advanceTimersByTimeAsync(100);
    v.draw();
    while (inflight.size) await Promise.all([Promise.allSettled([...inflight]), vi.advanceTimersByTimeAsync(200)]);
    quiet = began.n === before ? quiet + 1 : 0;
  }
  v.draw();
}

/** Both themes' exports settled, with no issue but the view's own, and different. */
async function bothThemes(v: SeqViewer, allowed: readonly string[] = []) {
  const out: string[] = [];
  for (const dark of [false, true]) {
    const r = v.renderSvg({ idPrefix: dark ? 'gd-' : 'gl-', tokens: dark ? DarkTokens : LightTokens });
    await vi.advanceTimersByTimeAsync(4000);
    const svg = await r;
    expect(svg.settled).toBe(true);
    expect(svg.issues.filter((i) => !allowed.includes(i.code)).map((i) => `${i.code}: ${i.message}`)).toEqual([]);
    out.push(svg.svg);
  }
  expect(out[0]).not.toBe(out[1]);
}

const layers = (v: SeqViewer) =>
  v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers.map((l) => ({ section: s, track: t, l }))));
const layer = <T>(v: SeqViewer, track: string, section?: string) =>
  layers(v).find((x) => x.track.track.label === track && (!section || x.section.section.address === section))?.l.lay as
    | T
    | undefined;
/** The element keys hit along a track's middle, a pixel at a time (what is drawn there, whatever the lays). */
function keysAlong(v: SeqViewer, track: string): string[] {
  const t = v.frame!.layout.sections.flatMap((s) => s.tracks).find((x) => x.track.label === track)!;
  const out = new Set<string>();
  for (let x = t.box.x; x < t.box.x + t.box.width; x += 0.5) {
    const k = v.hit(x, t.box.y + t.box.height / 2)?.element?.key;
    if (k !== undefined) out.add(String(k));
  }
  return [...out];
}

const bridges = (v: SeqViewer) => v.frame!.layout.bridges?.list ?? [];
const named = (v: SeqViewer, system: string) => v.view.model.value.systems.get(system)?.name;
const lengthOf = (v: SeqViewer, system: string) => v.view.model.value.systems.get(system)!.length!;

/** The centre of a grid lay's image tile `k`. */
const tileCentre = (g: GridLay, k: number): [number, number] => [
  g.ib[4 * k] + g.ib[4 * k + 2] / 2,
  g.ib[4 * k + 1] + g.ib[4 * k + 3] / 2,
];

describe('Examples/Genomes, headless', () => {
  it('synteny (EC-206 to EC-209, EC-218): a ribbon a block, the inversion crossed, the translocation onto B’s 2, no ribbon over chr1’s run of N, one map over every sequence pair, each ribbon its A chromosome’s colour; a ribbon hovered lights', async () => {
    const v = await open('synteny');
    const [b] = bridges(v);
    const ribbons = b.ribbons();
    expect(ribbons.length).toBeGreaterThan(50);
    const on = (r: (typeof ribbons)[number]) => [named(v, r.aSystem), named(v, r.bSystem)];
    const chr1 = ribbons.filter((r) => on(r)[0] === 'chr1');
    expect(chr1.some((r) => r.reversed)).toBe(true);
    expect(chr1.filter((r) => r.reversed).every((r) => r.aStart >= Inversion[0] && r.aEnd <= Inversion[1])).toBe(true);
    expect(ribbons.some((r) => on(r)[0] === 'chr3' && on(r)[1] === '2')).toBe(true);
    expect(chr1.some((r) => r.aEnd > NGap[0] && r.aStart < NGap[1])).toBe(false);
    // one map over several sequence pairs, none ambiguous (EC-207): chr1 to 1, chr2 to 2, chr3 to 3 and to 2
    const pairs = new Set(ribbons.map((r) => on(r).join(' to ')));
    for (const p of ['chr1 to 1', 'chr2 to 2', 'chr3 to 3', 'chr3 to 2']) expect(pairs.has(p), p).toBe(true);
    // (an ambiguous path would be an issue: `bothThemes` below allows none)
    // chromosome painting (EC-218): a ribbon's colour its A end's block's, one a chromosome and a different one each
    const lay = (b as unknown as { lays: { lay: { n: number; color: Uint32Array } }[] }).lays[0].lay;
    expect(lay.n).toBe(ribbons.length);
    const colours = new Map<string, Set<number>>();
    ribbons.forEach((r, i) => {
      const k = on(r)[0]!;
      if (!colours.has(k)) colours.set(k, new Set());
      colours.get(k)!.add(lay.color[i]);
    });
    for (const [k, c] of colours) expect(c.size, k).toBe(1);
    expect(new Set([...colours.values()].map((c) => [...c][0])).size).toBe(colours.size);
    // hovered: the ribbon lights (the hover pass draws)
    const r = ribbons.find((x) => x.ax1 - x.ax0 > 4)!;
    v.hover.set(v.hit((r.ax0 + r.ax1) / 2, b.box.y + 1), 'test');
    v.draw();
    expect(recordDrawLists(v.frame!.lists, ['hover']).split('\n').length).toBeGreaterThan(2);
    await bothThemes(v);
  });

  it('the dot plot (EC-109, EC-145, EC-209, EC-223): grid tiles of the alignment under its blocks; a cell is two positions, a drag selects a rectangle', async () => {
    const v = await open('dot-plot');
    const g = layers(v).find((x) => x.l.type.name === 'grid')!.l.lay as GridLay;
    expect(g.images.length).toBeGreaterThan(0);
    // its cells square where the room allows: the plot's box no wider than its sections' (EC-109)
    expect(g.box.width).toBeGreaterThan(0);
    // square cells where the room allows (EC-109): the box at most as tall as wide (aspect 1, capped by the room), and
    // both axes showing about 620 kb, a cell's shape the box's alone
    expect(g.box.height).toBeLessThanOrEqual(g.box.width);
    expect(g.axis!.scale / g.axisY!.scale).toBeCloseTo(g.box.width / g.box.height, 1);
    const [x, y] = tileCentre(g, 0);
    const h = v.hit(x, y);
    expect(h?.target).toBe('cell');
    // a cell is two positions (EC-145): one on an A sequence, one on a B sequence
    expect(['chr1', 'chr2', 'chr3']).toContain(named(v, h!.system));
    expect(['1', '2', '3']).toContain(named(v, h!.ySystem!));
    expect(Number.isFinite(h!.position) && Number.isFinite(h!.yPosition)).toBe(true);
    // a rectangle: a drag selects one range on each axis
    v.input.mac = false;
    expect(v.input.down(x - 10, y - 10, {})).toBe(true);
    v.input.move(x + 10, y + 10, {});
    v.input.up();
    expect(v.selection.value.value.location.rects.length).toBe(1);
    await bothThemes(v);
    // two levels: the fit reads the coarser, a zoom into chr1 against B's 1 the finer
    expect(g.binLevel).toBe(DotLevels[1]);
    v.view.setRange(ASection, 0, 30_000);
    v.view.setRange(BSection, 0, 30_000);
    await settle(v);
    const fine = layers(v).find((x) => x.l.type.name === 'grid')!.l.lay as GridLay;
    expect(fine.binLevel).toBe(DotLevels[0]);
    expect(fine.images.length + fine.n).toBeGreaterThan(0);
  });

  it('the signal by region (EC-169, EC-170, EC-175, EC-201, EC-202, EC-215): three tracks drawn from the indexed bedGraph at a coarse level, pending till read, the streamed tiles never in the spec', async () => {
    const v = await open('signal');
    const before = JSON.stringify(view('signal').build());
    for (const t of ['TPM (open)', 'TPM (from 0)', 'TPM (log)']) {
      const lay = layer<{ n: number; pending: number[]; domain?: [number, number] }>(v, t)!;
      expect(lay.pending, t).toEqual([]);
      expect(lay.n, t).toBeGreaterThan(100);
    }
    // the open domain is the data in view (EC-201), the one-sided one starts at 0 (EC-202)
    const openDomain = layer<{ domain?: [number, number] }>(v, 'TPM (open)')!.domain!;
    const zero = layer<{ domain?: [number, number] }>(v, 'TPM (from 0)')!.domain!;
    expect(openDomain[1]).toBeGreaterThan(100);
    expect(zero[0]).toBe(0);
    // the scale's open end a union over its layers in view (every block's), the open layer's its own block's
    expect(zero[1]).toBeGreaterThanOrEqual(openDomain[1]);
    // the log track's closed domain, TPM 1 to 10^4, its heights by the values' log10 (EC-215)
    expect(layer<{ domain?: [number, number] }>(v, 'TPM (log)')!.domain).toEqual([1, 10_000]);
    expect(layer<{ yLog: boolean }>(v, 'TPM (log)')!.yLog).toBe(true);
    expect(layer<{ yLog: boolean }>(v, 'TPM (open)')!.yLog).toBe(false);
    // read at a level above 1 at the fit (EC-170): a bin no wider than a device column's positions, never a line each
    const atFit = layer<{ axis: { scale: number } }>(v, 'TPM (open)')!.axis.scale;
    expect(Math.min(...levelsRead)).toBeGreaterThan(1);
    expect(Math.max(...levelsRead)).toBeLessThanOrEqual(Math.ceil(1 / atFit) * 4);
    // streamed data never enters the document (EC-175)
    expect(JSON.stringify(v.spec.root.value)).toBe(before);
    await bothThemes(v);
    // pending till read, never drawn as empty (EC-169): a window over chr2 not yet read shows its placeholder, then data
    v.view.setRange(ASection, 330_000, 340_000);
    v.draw();
    const waiting = layer<{ n: number; pending: number[] }>(v, 'TPM (open)')!;
    expect(waiting.pending.length).toBeGreaterThan(0);
    await settle(v);
    const read = layer<{ n: number; pending: number[] }>(v, 'TPM (open)')!;
    expect(read.pending).toEqual([]);
    expect(read.n).toBeGreaterThan(0);
  });

  it('samples as rows (EC-186, EC-220, EC-222): 40 samples by the VCF’s named index, the 600-base deletion solid in its carriers’ rows; a cell hit and hovered, every site mark hit', async () => {
    const v = await open('samples');
    const m = layer<MatrixLay>(v, 'Genotypes')!;
    expect(m.states).toBe(SampleCount);
    expect(m.n).toBeGreaterThan(12);
    const sites = layer<RangesLay>(v, 'Sites')!;
    expect(sites.pending).toEqual([]);
    // the deletion's carriers' rows: a hit across the deletion's middle in the first carrier's row is the deletion
    const siteKeys = Array.from({ length: sites.items }, (_, i) => String(sites.keys[i]));
    expect(siteKeys.some((k) => k.startsWith(`${DeletionExtent[0] + 1}:`))).toBe(true);
    const t = v.frame!.layout.sections[0].tracks.find((x) => x.track.label === 'Genotypes')!;
    const row = t.box.height / SampleCount;
    const axis = v.frame!.layout.sections[0].axis!;
    const mid = (DeletionExtent[0] + DeletionExtent[1]) / 2;
    // a cell hit is the sample's row and the allele: in a carrier's row across the deletion's middle, the deletion
    const at = (sample: number) => v.hit(axis.x(mid), t.box.y + (sample + 0.5) * row);
    const carrier = at(DeletionCarriers[0]);
    expect(carrier?.layer).toContain('/gt');
    expect(carrier?.target).toBe('cell');
    expect(carrier?.yPosition).toBe(DeletionCarriers[0]);
    // the deletion's allele across its carriers' rows (the first homozygous), the reference in a non-carrier's
    const sys = (x: string) => named(v, x) ?? x;
    const allele = (sample: number) => matrixCellHover(m, mid, sample, sys)?.rows.find((r) => r[0] === 'Value')?.[1];
    // its REF: the anchor letter and the 600 deleted (VCF's padding base)
    const ref = [DeletionExtent[0], DeletionExtent[1] + 1];
    const deletion = `chr1:${spanText(ref[0], ref[1])}`;
    for (const c of DeletionCarriers) {
      expect(matrixCellHover(m, mid, c, sys)?.rows[0][1], `sample ${c}`).toContain(`(${deletion})`);
      expect(allele(c), `sample ${c}`).toBe(c === DeletionCarriers[0] ? '1/1' : '0/1');
    }
    const other = Array.from({ length: SampleCount }, (_, i) => i).find((i) => !DeletionCarriers.includes(i))!;
    expect(allele(other)).toBe('0/0');
    // a click on a carrier's cell selects the deletion, its extent on chr1 (as a ranges click does)
    const [cx, cy] = [axis.x(mid), t.box.y + (DeletionCarriers[0] + 0.5) * row];
    v.input.move(cx, cy);
    v.input.down(cx, cy, {});
    v.input.up();
    v.input.click(cx, cy, {});
    const picked = v.selection.value.value.location;
    expect(picked.rects).toEqual([]);
    expect(picked.parts.map((p) => [named(v, p.system), ...p.ranges])).toEqual([['chr1', ...ref]]);
    v.selection.clear();
    // hovered: the cell lights (the hover pass draws)
    v.hover.set(carrier, 'test');
    v.draw();
    expect(recordDrawLists(v.frame!.lists, ['hover']).split('\n').length).toBeGreaterThan(2);
    v.hover.set(undefined, 'test');
    const { wrong, marks } = checkMarks(v);
    expect(wrong).toEqual([]);
    expect(marks).toBeGreaterThan(5);
    await bothThemes(v);
  });

  it('methylation (EC-220, EC-221): the CpGs’ matrix of 20 samples and the gene-by-sample heatmap drawn', async () => {
    const v = await open('methylation');
    const cpg = layer<MatrixLay>(v, 'CpG methylation')!;
    const genes = layer<MatrixLay>(v, 'Genes by sample')!;
    expect(cpg.states).toBe(20);
    expect(genes.states).toBe(20);
    expect(cpg.n + cpg.images.length).toBeGreaterThan(0);
    expect(genes.n + genes.images.length).toBeGreaterThan(0);
    await bothThemes(v);
  });

  it('transcripts in three levels (EC-227, EC-229): genes, transcripts and the GTF’s, read by region, the transcripts read whole alike; each mark hit and hovered; a query over a gene not read yet waiting', async () => {
    const v = await open('transcripts');
    for (const t of ['Genes (GFF3)', 'Transcripts (GFF3)', 'Transcripts (GTF)', 'Transcripts (GFF3, whole)']) {
      const lay = layer<RangesLay>(v, t)!;
      expect(lay.pending, t).toEqual([]);
      expect(lay.items, t).toBeGreaterThan(0);
    }
    const tx = layer<RangesLay>(v, 'Transcripts (GFF3)')!.items;
    expect(layer<RangesLay>(v, 'Transcripts (GTF)')!.items).toBe(tx);
    // read whole, the same transcripts in the window as by region
    const whole = layer<RangesLay>(v, 'Transcripts (GFF3, whole)')!;
    expect(whole.pending).toEqual([]);
    expect(whole.items).toBe(tx);
    expect(tx).toBeGreaterThan(layer<RangesLay>(v, 'Genes (GFF3)')!.items);
    const { wrong } = checkMarks(v);
    expect(wrong).toEqual([]);
    // a query over the genes read by region (EC-229): a gene whose tile has not arrived (far down chr1) is waiting,
    // never missing
    const far = v.addresses.locate('/data/genes/(item :id "gene:TOY0030")');
    expect(far.resolution.status).toBe('waiting');
    expect(far.location.parts).toEqual([]);
    await bothThemes(v);
  });

  it('open on a gene by name (EC-228): the GFF3’s name index finds TOY0042, its extent focused', async () => {
    const v = await open('by-name');
    const hits = await findGene(OpenedGene, reader, new AbortController().signal);
    expect(hits[0].label).toBe(OpenedGene);
    v.runCommand('address.focus', { address: geneRange(hits[0]) });
    await settle(v);
    const r = v.view.get(GenomeSection)!;
    expect(r.end - r.start).toBeLessThan(20_000);
  });

  it('names aliased (EC-186 to EC-189, EC-200): toy A’s lengths from the adapter, marks named 1, an accession and chr3 drawn, MT refused; the HLA contig from the FASTA', async () => {
    const v = await open('names');
    // the marks drawn, by their ids: every pixel along the track hit
    const drawn = keysAlong(v, 'Named elsewhere');
    // `1` by the `chr` rule, `NC_900002.1` by the alias table, `chr3` by its name, the HLA contig by its own (EC-189,
    // EC-200); `MT` refused, never guessed onto chrM (said: the layout's issue)
    expect(drawn.sort()).toEqual(['by accession', 'by name', 'by number', 'on the HLA contig']);
    expect(AliasedMarks.target).toContain('MT');
    // toy A's lengths are the adapter's (no manifest, EC-188); the HLA contig's from the FASTA, read by its name
    const systems = v.view.model.value.systems;
    const hla = [...systems.values()].filter((x) => x.name === 'HLA-A*01:01:01:01');
    // twice: toy A's (the adapter's) and the small FASTA's, each 3,503 long
    expect(hla.map((x) => x.length)).toEqual([3_503, 3_503]);
    expect([...systems.values()].some((x) => x.name === 'chr1' && x.length === 300_000)).toBe(true);
    await bothThemes(v, ['seq.resolve.assembly.name']);
  });

  it('the contact map (EC-223): every block pair’s tiles drawn, the lower triangle the upper’s mirror', async () => {
    const v = await open('contact-map');
    const g = layers(v).find((x) => x.l.type.name === 'grid')!.l.lay as GridLay;
    // read at one of the map's levels, its tiles drawn as images with nothing pending (EC-169, EC-223)
    expect([4096, 16_384]).toContain(g.binLevel);
    expect(g.pend).toEqual([]);
    expect(g.images.length).toBeGreaterThan(0);
    // counts on a log palette: an empty cell (0, no contact) left clear, every count from 1 coloured
    let clear = 0;
    let drawn = 0;
    for (const img of g.images) for (let q = 3; q < img.data.length; q += 4) img.data[q] ? drawn++ : clear++;
    expect(clear).toBeGreaterThan(0);
    expect(drawn).toBeGreaterThan(clear);
    // every block pair's lay (chr1 to chr3 on each axis), by its systems' names
    const grids = v
      .frame!.layout.sections.flatMap((s) => s.lays.flat())
      .filter((l) => l.type.name === 'grid')
      .map((l) => l.lay as GridLay);
    const pair = (x: string, y: string) => grids.find((l) => named(v, l.system) === x && named(v, l.ySystem) === y)!;
    expect(grids.length).toBe(9);
    // the lower triangle the upper's mirror: chr3 × chr1's cell is chr1 × chr3's transposed, the hotspot's among them
    const sys = (s: string) => named(v, s) ?? s;
    const value = (l: GridLay, i: number, j: number) =>
      gridCellHover(l, i, j, sys)?.rows.find((r) => r[0] === 'Value')?.[1];
    const [upper, lower] = [pair('chr1', 'chr3'), pair('chr3', 'chr1')];
    for (const [i, j] of [
      [Hotspot[0], Hotspot[1]],
      [10_000, 100_000],
      [290_000, 5_000],
    ]) {
      const u = value(upper, i, j);
      expect(u, `${i}, ${j}`).toMatch(/^[\d,.]+$/);
      expect(value(lower, j, i), `${i}, ${j}`).toBe(u);
    }
    expect(Number(value(upper, Hotspot[0], Hotspot[1])!.replace(/,/g, ''))).toBeGreaterThan(100);
    // the tiles restored at a block's end (EC-223): each pair's last cell, at both systems' ends, has its value, and
    // its tiles' clips reach the pair's box's right and bottom edges
    for (const l of grids) {
      const [nx, ny] = [lengthOf(v, l.system), lengthOf(v, l.ySystem)];
      expect(value(l, nx - 1, ny - 1), `${named(v, l.system)} × ${named(v, l.ySystem)}`).toMatch(/^[\d,.]+$/);
      let right = Number.NEGATIVE_INFINITY;
      let bottom = Number.NEGATIVE_INFINITY;
      for (let k = 0; k + 3 < l.ic.length; k += 4) {
        right = Math.max(right, l.ic[k] + l.ic[k + 2]);
        bottom = Math.max(bottom, l.ic[k + 1] + l.ic[k + 3]);
      }
      expect(right).toBeGreaterThan(l.box.x + l.box.width - 1);
      expect(bottom).toBeGreaterThan(l.box.y + l.box.height - 1);
    }
    await bothThemes(v, ['seq.render.scale.log']);
  });

  it('the liftover (EC-211, EC-212): A’s genes drawn on A′ through the map read by region, the duplicated stretch twice', async () => {
    const v = await open('liftover');
    const [b] = bridges(v);
    const ribbons = b
      .ribbons()
      .map((r) => ({ a: named(v, r.aSystem), b: named(v, r.bSystem), s: r.aStart, e: r.aEnd }));
    // chr1's duplicated 10 kb drawn twice (EC-212), chr3's lost 4 kb no ribbon
    const dup = ribbons.filter((r) => r.a === 'chr1' && r.s <= Duplicated[0] && r.e >= Duplicated[1]);
    expect(dup.length).toBe(2);
    expect(ribbons.some((r) => r.a === 'chr3' && r.s < 34_000 && r.e > 30_000)).toBe(false);
    // A's genes drawn on A′ through the map read by region (EC-211)
    expect(keysAlong(v, 'Genes (lifted from A)').length).toBeGreaterThan(0);
    await bothThemes(v);
  });
});
