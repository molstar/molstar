/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import {
  delayedSource,
  type FormatRuntime,
  openSeq,
  type OpenedSource,
  type RegionRequest,
  type SeqData,
  SourceOpenError,
  type StoreClock,
} from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome, theme } from '@molstar/seqx-core/lang';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import type { ExampleView } from '../views.js';
import { HicAssemblyName, HicChromosomes, hicOneLevel, hicSource } from './hic.js';
import { openBedpeByRegion } from './bedpe.js';
import { genotypesSource } from './genotypes.js';
import { LinkCount, LinksChromosome, RingChromosomes, linksOnChromosome, linksSource } from './links.js';
import { methylationSource } from './methylation.js';
import { SyntenyA, SyntenyAssemblyA, SyntenyAssemblyB, SyntenyB, syntenySource } from './synteny.js';

/** The blocks of the synteny the zoomed-in view reads (about 80 on chrA2: an arc each stays readable). */
const SyntenyViewBlocks = 240;

// The region reads' views (Features/Region reads): each a source generated at open (no
// file) and answered after a delay on the store's clock, the finest level later than the coarse ones, so a picture
// arrives coarse, then sharp. A source names the format `seq-generated` and the generator by `uri` (`generated:hic`).

/** The format the views' sources name, and the generators its runtime knows by `uri`. */
export const GeneratedFormat = 'seq-generated';

/** Every generator a source may name: its `seq` data (the BEDPE links are text, read by the `bed` runtime). */
const Generators: Readonly<Record<string, () => SeqData>> = {
  'generated:hic': hicSource,
  'generated:hic-one-level': hicOneLevel,
  'generated:synteny': () => syntenySource(SyntenyViewBlocks, undefined, { sided: false }),
  'generated:links': () => linksSource(),
  'generated:links-ring': () => linksSource(RingChromosomes),
  'generated:genotypes': () => genotypesSource(),
  'generated:genotypes-items': () => genotypesItemsOnly(),
  'generated:methylation': () => methylationSource(),
};

/** The BEDPE links (one chromosome) a source reads as a file's text. */
const BedpeUri = 'generated:links.bedpe';

/** The genotypes at the items' level alone: past a few kb a view cannot read it (the zoom-in state). */
function genotypesItemsOnly(): SeqData {
  const full = genotypesSource();
  return {
    coverage: { ...full.coverage, levels: [1] },
    pairs: full.pairs.filter((p) => p.request.level === 1),
  };
}

export interface GeneratedOptions {
  /** The clock the answers wait on: real timers by default. */
  readonly clock?: StoreClock;
  /** The finest level's delay, ms (coarser levels a third of it); each answer's delay by request when a function. */
  readonly delayMs?: number | ((request: RegionRequest) => number);
}

/** A delay the finest level waits and the coarser levels a third of: tiles arrive coarse, then sharp. */
function byLevel(finest: number, delayMs: number) {
  return (r: RegionRequest) => (r.level <= finest ? delayMs : Math.ceil(delayMs / 3));
}

/**
 * The runtime of `seq-generated` sources: `uri` names a generator (`generated:hic`), the data opened as a `seq` source
 * and answered after a delay. The BEDPE links open as the `bed` runtime opens a text.
 */
export function generatedRuntime(options: GeneratedOptions = {}): FormatRuntime {
  return {
    name: 'seq-generated',
    formats: [GeneratedFormat],
    async open(source, ctx): Promise<OpenedSource> {
      if (source.uri === BedpeUri) return openBedpeByRegion(source, ctx, LinksChromosome);
      const make = source.uri === undefined ? undefined : Generators[source.uri];
      if (!make) throw new SourceOpenError([]);
      const data = make();
      const delay = options.delayMs ?? 300;
      const finest = data.coverage.levels[0];
      return delayedSource(openSeq(data), {
        ...(options.clock ? { clock: options.clock } : {}),
        delayMs: typeof delay === 'number' ? byLevel(finest, delay) : delay,
      });
    },
  };
}

/** What a view needs besides the default registry: the generated runtime (its BEDPE links opened by the `bed` runtime). */
export function registry(options: { readonly delayMs?: number } = {}): SeqRegistryEntry[] {
  return [
    { kind: 'format-runtime', runtime: generatedRuntime(options) },
    { kind: 'format-runtime', runtime: BedRuntime },
  ];
}

const generated = (uri: string) => ({ format: GeneratedFormat, uri }) as never;

/** The sections' addresses (what the tests read). */
export const SyntenySection = '/views/main/synteny';
export const LinksSection = '/views/main/links';
export const RingSection = '/views/main/ring';
export const GenomeSection = '/views/main/genome';
export const HicSection = '/views/main/chr1';
export const HicPlot = '/views/main/plot';
export const GenotypesSection = '/views/main/chr';

/** The two toy genomes of the synteny as assemblies (their chromosomes' lengths inline). */
function syntenyAssemblies(b: ReturnType<typeof createSeqSpec>) {
  const a = b.data.assembly('a', assembly(SyntenyAssemblyA, SyntenyA, { order: 'source' }));
  const g = b.data.assembly('b', assembly(SyntenyAssemblyB, SyntenyB, { order: 'source' }));
  return { a, g };
}

/**
 * The synteny between two toy genomes, on one section of both assemblies kept to `chrA2` and `chrB2` (the section's
 * `include`): an arc a block of the alignment between them (the ones B reads reversed, a 6 Mb stretch, are twisted),
 * about 80 of them. The pairs are over the assemblies, read a chromosome at a time.
 */
export function pairsZoomedIn(): Node {
  const b = createSeqSpec({ title: 'Synteny, zoomed in' });
  const src = b.data.source('alignment', generated('generated:synteny'));
  const { a, g } = syntenyAssemblies(b);
  const pairs = b.data.pairs('synteny', { a, b: g, source: src });
  b.views.view('main', (v) => {
    genome(v, 'synteny', [a, g], { include: ['chrA2', 'chrB2'] }, (s) => {
      s.track('synteny', { label: 'Synteny', height: 120 }).layer('synteny', {
        type: 'pairs',
        data: pairs,
        legend: 'none',
      });
    });
  });
  return b.build();
}

/** The links' chromosome as a sequence of a view: its length from an inline chrom.sizes. */
function linksChromosome(b: ReturnType<typeof createSeqSpec>) {
  const sizes = b.data.source('sizes', {
    format: 'chrom-sizes',
    value: `${LinksChromosome.name}\t${LinksChromosome.length}\n` as never,
  });
  return { sizes, chr: b.data.sequence('chr', { alphabet: 'dna', source: sizes, name: LinksChromosome.name }) };
}

/**
 * 12,000 links on one chromosome, zoomed out: bin pairs, each an arc stroked by its count, the hover naming
 * the bin pair. The same links twice, held (the table in the spec) and as a BEDPE: one picture.
 */
export function pairsZoomedOut(): Node {
  const b = createSeqSpec({ title: `${LinkCount.toLocaleString('en-US')} links, zoomed out` });
  const { chr } = linksChromosome(b);
  const t = linksOnChromosome();
  const held = b.data.pairs('held', {
    a: chr,
    rows: {
      a_start: t.aStart,
      a_end: t.aEnd,
      b_start: t.bStart,
      b_end: t.bEnd,
      id: t.id!,
    },
  });
  const bedpe = b.data.pairs('bedpe', {
    a: chr,
    source: b.data.source('links-bedpe', generated(BedpeUri)),
  });
  b.views.view('main', (v) => {
    const s = v.section('links', { of: chr });
    s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('held', { label: 'Held', height: 120 }).layer('held', {
      type: 'pairs',
      data: held,
      legend: 'none',
    });
    s.track('bedpe', { label: 'BEDPE', height: 120 }).layer('bedpe', {
      type: 'pairs',
      data: bedpe,
      legend: 'none',
    });
  });
  return b.build();
}

/**
 * The links between a toy assembly's six chromosomes on a ring at the fit (the genome ring): bin-pair chords between
 * the chromosomes and within each, stroked by count, the strongest drawn first. The links are over the assembly, read
 * a chromosome at a time from a source with three levels.
 */
export function pairsRing(): Node {
  const b = createSeqSpec({ title: 'Links on a genome ring' });
  const g = b.data.assembly('toy', assembly(HicAssemblyName, RingChromosomes, { order: 'source' }));
  const src = b.data.source('links-source', generated('generated:links-ring'));
  const links = b.data.pairs('links', { a: g, source: src });
  b.views.view('main', (v) => {
    genome(v, 'ring', g, { layout: 'rings' }, (s) => {
      s.track('links', { label: 'Links' }).layer('links', { type: 'pairs', data: links, legend: 'none' });
    });
  });
  return b.build();
}

/** The Hi-C's chromosome (`chr1`) as a sequence read from the source. */
/** The contacts' colour: counts on the blues, the domain the range of the tiles drawn. */
const counts = () => theme.field('value', { palette: theme.continuous({ colors: 'blues' }) });

/** One chromosome's Hi-C by region: tiles of its two levels, a pending tile drawn from the coarser one. */
export function hicByRegion(): Node {
  const b = createSeqSpec({ title: 'Hi-C by region' });
  const src = b.data.source('hic', generated('generated:hic'));
  const chr = b.data.sequence('chr1', { alphabet: 'dna', source: src, name: 'chr1' });
  const contacts = b.data.grid('contacts', { x: chr, source: src });
  b.views.view('main', (v) => {
    const top = v.section('chr1', { of: chr });
    top.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' });
    v.section('plot', { layout: 'grid', of: chr, follow: top.address }).track('plot', {}, (t) => {
      t.layer('contacts', { type: 'grid', data: contacts, color: counts() });
    });
  });
  return b.build();
}

/**
 * The Hi-C over a toy assembly's six chromosomes at the fit: the plot follows the genome's section, a block pair a
 * cell of chromosomes, each read by block pair from the source's tiles (the upper triangle served, the lower its
 * mirror), with the chromosome gaps and labels on both axes.
 */
export function hicOverGenome(): Node {
  const b = createSeqSpec({ title: 'Hi-C over the genome' });
  const g = b.data.assembly('toy', assembly(HicAssemblyName, HicChromosomes, { order: 'source' }));
  const contacts = b.data.grid('contacts', { x: g, source: b.data.source('hic', generated('generated:hic')) });
  b.views.view('main', (v) => {
    const top = genome(v, 'genome', g);
    v.section('plot', { layout: 'grid', follow: top.address }).track('plot', {}, (t) => {
      t.layer('contacts', { type: 'grid', data: contacts, color: counts() });
    });
  });
  return b.build();
}

/** 40 samples over a chromosome's variants, a row a sample, a long deletion solid across its SNVs. */
export function matrixSamples(): Node {
  const b = createSeqSpec({ title: 'Samples by positions' });
  const src = b.data.source('genotypes', generated('generated:genotypes'));
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const gt = b.data.matrix('gt', { target: chr, source: src });
  b.views.view('main', (v) => {
    v.section('chr', { of: chr }, (s) => {
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' });
      s.track('gt', { label: 'Genotypes', height: 260 }).layer('gt', { type: 'matrix', data: gt });
    });
  });
  return b.build();
}

/** the same matrix on a ring, each sample a thin ring. */
export function matrixRing(): Node {
  const b = createSeqSpec({ title: 'A matrix on a ring' });
  const src = b.data.source('genotypes', generated('generated:genotypes'));
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const gt = b.data.matrix('gt', { target: chr, source: src });
  b.views.view('main', (v) => {
    v.section('chr', { of: chr, layout: 'rings' }, (s) => {
      s.track('gt', { label: 'Genotypes', height: 220 }).layer('gt', { type: 'matrix', data: gt });
    });
  });
  return b.build();
}

/** Methylation with values missing: a value a CpG and sample, missing ones left clear. */
export function matrixMethylation(): Node {
  const b = createSeqSpec({ title: 'Methylation with missing values' });
  const src = b.data.source('cpgs', generated('generated:methylation'));
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const m = b.data.matrix('methylation', { target: chr, source: src });
  b.views.view('main', (v) => {
    v.section('chr', { of: chr }, (s) => {
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' });
      s.track('methylation', { label: 'Methylation', height: 120 }).layer('methylation', { type: 'matrix', data: m });
    });
  });
  return b.build();
}

/**
 * a matrix served at its items' level alone, past what a view reads ("zoom in to read this track"), and
 * links over a stretch with no link end (nothing drawn: no arc, no placeholder).
 */
export function zoomIn(): Node {
  const b = createSeqSpec({ title: 'Zoom in to read, and bins with nothing' });
  const src = b.data.source('genotypes', generated('generated:genotypes-items'));
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const gt = b.data.matrix('gt', { target: chr, source: src });
  const links = b.data.pairs('links', { a: chr, rows: emptyLinks() });
  b.views.view('main', (v) => {
    v.section('chr', { of: chr }, (s) => {
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' });
      s.track('gt', { label: 'Genotypes', height: 100 }).layer('gt', { type: 'matrix', data: gt });
      s.track('links', { label: 'Links', height: 60 }).layer('links', {
        type: 'pairs',
        data: links,
        legend: 'none',
      });
    });
  });
  return b.build();
}

/** Forty links inside the first 200 kb of the genotypes' chromosome: none reaches 400 to 700 kb, where the view opens. */
function emptyLinks() {
  const n = 40;
  const start = Int32Array.from({ length: n }, (_, i) => 10_000 + i * 2_000);
  const far = Int32Array.from({ length: n }, (_, i) => 100_000 + i * 2_400);
  return {
    a_start: start,
    a_end: start.map((s) => s + 10),
    b_start: far,
    b_end: far.map((s) => s + 10),
  };
}

export const RegionReadsViews: readonly ExampleView[] = [
  {
    dataset: 'regions',
    id: 'pairs-zoomed-in',
    title: 'Pairs zoomed in: items as arcs',
    summary:
      'A synteny alignment between two toy chromosomes (generated, pairs in the `seq` format): chrA2 on the left, chrB2 on the right, about 80 pairs in view, far under the 4,096 an arc a pair is drawn from, so each alignment block is one arc between its two ends; the ones B reads reversed (a stretch of 6 Mb) are twisted. Each pair is found from either end.',
    build: pairsZoomedIn,
    registry,
  },
  {
    dataset: 'regions',
    id: 'pairs-zoomed-out',
    title: 'Pairs zoomed out: bin pairs by count',
    summary:
      '12,000 generated links on a 2 Mb chromosome, read twice, held in the spec and as a BEDPE file: past the 4,096 arcs a view draws, the BEDPE’s pairs are read as bin pairs (the `bed` runtime’s pairs levels, binned per tile from the rows it holds), each an arc stroked by its count; the held table, not read through the store, bins its arcs by device column. Hover a bin pair: it names both bins and the count, never a pair.',
    build: pairsZoomedOut,
    registry,
  },
  {
    dataset: 'regions',
    id: 'pairs-ring',
    title: 'A genome ring at the fit: bin-pair chords',
    summary:
      'The genome ring at the fit: 12,000 generated links over a toy assembly’s six chromosomes (38 to 10 Mb), read a chromosome at a time as bin pairs of 16 kb: chords between the chromosomes and within each, stroked and shaded by count, the strongest kept and drawn first.',
    build: pairsRing,
    registry,
  },
  {
    dataset: 'regions',
    id: 'hic-region',
    title: 'A Hi-C map by region',
    summary:
      'A toy chromosome’s symmetric contact map (counts, two levels) read in tiles by the section: a pending tile is drawn from the coarser one held (blocky, then sharp), the plot’s box fixed; the lower triangle is the upper one’s mirror.',
    build: hicByRegion,
    registry,
  },
  {
    dataset: 'regions',
    id: 'hic-genome',
    title: 'A Hi-C map over the genome',
    summary:
      'The map of a toy assembly’s six chromosomes as one system of 1,376 bins (one level served, the store’s summaries make the others): at the fit one read of the whole-assembly level.',
    build: hicOverGenome,
    registry,
  },
  {
    dataset: 'regions',
    id: 'matrix-samples',
    title: 'Samples by positions',
    summary:
      '40 generated samples over a chromosome’s variants, a row a sample, a mark a variant at its extent coloured by its genotype (uncalled in the missing colour); one long deletion is solid across the SNVs it covers; zoomed out, a bin’s majority genotype with its share as its strength.',
    build: matrixSamples,
    registry,
  },
  {
    dataset: 'regions',
    id: 'matrix-ring',
    title: 'A matrix on a ring',
    summary: 'The same 40 samples as 40 thin rings inside the chromosome’s arc: at most a run a device column a ring.',
    build: matrixRing,
    registry,
  },
  {
    dataset: 'regions',
    id: 'matrix-methylation',
    title: 'Methylation with values missing',
    summary:
      'A methylation fraction per CpG and sample (eight samples, a CpG island where it drops), about a tenth of the values missing: the missing ones are left clear, never zero.',
    build: matrixMethylation,
    registry,
  },
  {
    dataset: 'regions',
    id: 'zoom-in',
    title: 'Zoom in to read, and bins with nothing',
    summary:
      'A matrix served at its items’ level alone (past what a view reads: its band hatched and saying “zoom in to read this track”, an info issue, no wait resolving it) and links with no end in the stretch in view (nothing drawn: no arc, no placeholder). Zoom in to 20 kb and the genotypes are read.',
    build: zoomIn,
    registry,
  },
];
