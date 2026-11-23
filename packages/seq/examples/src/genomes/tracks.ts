/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { type AssetReader, DefaultOpenLimits, type NameHit } from '@molstar/seqx-core/data';
import { createSeqSpec, genome, theme } from '@molstar/seqx-core/lang';
import { GffRuntime } from '@molstar/seqx-core/runtimes/gff';
import { VcfCategories } from '@molstar/seqx-core/runtimes/vcf';
import { CpgStretch } from './generate/methylation.js';
import { DeletionExtent } from './generate/variants.js';
import { assemblyA, type Files, fileSource, GenomesFiles } from './toys.js';

// The Genomes dataset's views of toy A by region: the signal of a bedGraph read by
// region through its index (an open domain, a domain open at one end, a log scale), the samples' variants with their
// genotype matrix (a long deletion solid), methylation by region and the gene-by-sample heatmap, the genes'
// transcripts in three levels (a GFF3 by region and read whole, the GTF of the same model beside it), and a GFF3 read
// whole whose genes are found by name. Each reads the committed files (`files`, the dataset's by default).

/** The sections' addresses (what the tests read). */
export const GenomeSection = '/views/main/toy-a';
export const Chr1Section = '/views/main/chr1';

/** The signal's sources and layers (`values` over toy A, from the indexed bedGraph). */
const SignalTracks = ['open', 'floor-up', 'log'] as const;
/** The one-sided row's set low end: TPM 20, above most of the data's floor (its bars stand on it, EC-202). */
export const SignalFloor = 20;

/**
 * The signal over toy A by region (EC-169 to EC-171, EC-175): three tracks of one bedGraph read through its tabix
 * index, at the coarsest of its levels that fits a device column: its domain open (the values in view down to zero,
 * EC-201), open at its upper end only, its low end set at `SignalFloor` (`[20, null]`, EC-202: its baseline there, the
 * values below it flat), and on a log scale (TPM over 0 to 10^4, EC-215). The missing stretches (chr1's
 * run of N, 10 kb of chr2) are gaps, never zero.
 */
export function signal(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Signal by region' });
  const a = assemblyA(b, files);
  const values = b.data.values('signal', {
    target: a,
    source: fileSource(b, 'signal-file', 'bedgraph', files.signal, files.signalIndex),
  } as never);
  const floorUp = b.views.scale('floor-up', {
    label: 'TPM',
    domain: [SignalFloor, null] as unknown as [number, number],
  });
  const log = b.views.scale('log', { label: 'TPM (log)', type: 'log', domain: [1, 10_000] } as never);
  b.views.view('main', (v) =>
    genome(v, 'toy-a', a, { include: ['chr1', 'chr2', 'chr3'] }, (s) => {
      s.track(SignalTracks[0], { label: 'TPM (open)', height: 48 }).layer('values', { type: 'values', data: values });
      s.track(SignalTracks[1], { label: `TPM (from ${SignalFloor})`, height: 48 }).layer('values', {
        type: 'values',
        data: values,
        scale: floorUp,
      } as never);
      s.track(SignalTracks[2], { label: 'TPM (log)', height: 48 }).layer('values', {
        type: 'values',
        data: values,
        scale: log,
        color: theme.field('value', { scale: log }),
      } as never);
    }),
  );
  return b.build();
}

/** An allele's colour by its category: every category listed (a complete key). */
const byCategory = () => theme.field('category', { palette: theme.categorical({ values: [...VcfCategories] }) });

/**
 * Toy A's chr1 around the long deletion, read by region through the VCF's tabix index (EC-186): the sites (a row an
 * ALT allele) and the genotype matrix of 40 samples below (a row a sample, coloured by genotype class; uncalled in the
 * missing colour, EC-220), the 600-base deletion one item solid across the 12 SNVs it covers in its carriers’ rows
 * (EC-222).
 */
export function samples(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Samples’ variants as rows' });
  const a = assemblyA(b, files);
  const file = fileSource(b, 'vcf-file', 'vcf', files.vcf, files.vcfIndex);
  const sites = b.data.ranges('variants', { target: a, source: file } as never);
  const gt = b.data.matrix('gt', { target: a, source: file } as never);
  const [d0, d1] = DeletionExtent;
  b.views.view('main', (v) =>
    genome(v, 'chr1', a, { include: ['chr1'] }, (s) => {
      // chr1's block kept to the deletion and about 2 kb each side (a block overriding the derived one by its name)
      s.block('chr1', { range: [d0 - 2_000, d1 + 1_800] } as never);
      s.track('variants', { label: 'Sites', height: 48 }).layer('variants', {
        type: 'ranges',
        data: sites,
        color: byCategory(),
      });
      s.track('gt', { label: 'Genotypes', height: 280 }).layer('gt', { type: 'matrix', data: gt });
    }),
  );
  return b.build();
}

/**
 * Methylation by region (EC-220): 20 samples at 250 CpGs of chr1's 50 kb stretch, a fraction per CpG and sample
 * (missing values left clear, never 0), the island where it falls; and the same matrix layer over the 150 genes, each
 * gene's mean a sample (a gene-by-sample heatmap, EC-221), over the genome.
 */
export function methylation(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Methylation by region' });
  const a = assemblyA(b, files);
  const cpgs = b.data.matrix('cpgs', {
    target: a,
    source: fileSource(b, 'methylation-file', 'seq', files.methylation),
  } as never);
  const genes = b.data.matrix('genes', {
    target: a,
    source: fileSource(b, 'gene-methylation-file', 'seq', files.geneMethylation),
  } as never);
  const fraction = () => theme.field('value', { domain: [0, 1], palette: theme.continuous({ colors: 'viridis' }) });
  b.views.view('main', (v) => {
    genome(v, 'chr1', a, { include: ['chr1'] }, (s) => {
      s.block('chr1', { range: [...CpgStretch] } as never);
      s.track('cpgs', { label: 'CpG methylation', height: 160 }).layer('cpgs', {
        type: 'matrix',
        data: cpgs,
        color: fraction(),
      } as never);
    });
    genome(v, 'toy-a', a, { include: ['chr1', 'chr2', 'chr3'] }, (s) =>
      s.track('genes', { label: 'Genes by sample', height: 160 }).layer('genes', {
        type: 'matrix',
        data: genes,
        color: fraction(),
      } as never),
    );
  });
  return b.build();
}

/**
 * The genes in three levels (EC-227): the GFF3 read by region through its index, its genes (`features: genes`, one
 * model a gene) and its transcripts (gene → transcript → exon: a model a transcript, its exons boxes, its CDS thick, a
 * hairline with chevrons on the minus strand), the GTF of the same model beside them, and the GFF3's transcripts read
 * whole (the design's "whole and by region": the same lanes either way).
 */
export function transcripts(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Transcripts in three levels' });
  const a = assemblyA(b, files);
  const gff = (id: string, features?: string) =>
    fileSource(b, id, 'gff3', files.gff3, files.gff3Index, features ? { features } : {});
  const genes = b.data.ranges('genes', { target: a, source: gff('gff3-genes', 'genes') } as never);
  const tx = b.data.ranges('transcripts', { target: a, source: gff('gff3-transcripts') } as never);
  const gtf = b.data.ranges('gtf', {
    target: a,
    source: fileSource(b, 'gtf-transcripts', 'gtf', files.gtf, files.gtfIndex),
  } as never);
  // the GFF3's transcripts read whole too (no index: the runtime's other path), the same lanes as by region
  const whole = b.data.ranges('whole', { target: a, source: fileSource(b, 'gff3-whole', 'gff3', files.gff3) } as never);
  const model = (data: string) => ({ type: 'ranges', data, line: 'gene' }) as never;
  b.views.view('main', (v) =>
    genome(v, 'toy-a', a, { include: ['chr1'] }, (s) => {
      // chr1's first 8 kb: two genes, their transcripts in lanes (under 16 positions a pixel: items, not density)
      s.block('chr1', { range: [0, 8_000] } as never);
      s.track('genes', { label: 'Genes (GFF3)' }).layer('genes', model(genes));
      s.track('transcripts', { label: 'Transcripts (GFF3)' }).layer('transcripts', model(tx));
      s.track('gtf', { label: 'Transcripts (GTF)' }).layer('gtf', model(gtf));
      s.track('whole', { label: 'Transcripts (GFF3, whole)' }).layer('whole', model(whole));
    }),
  );
  return b.build();
}

/** The gene the by-name view opens on. */
export const OpenedGene = 'TOY0042';

/**
 * A GFF3 read whole, its genes found by name (EC-228): the runtime's name index over `Name` and `ID` (exact names
 * first, then by prefix), the view opened on `TOY0042`'s extent (the story's find box and its test ask the index and
 * focus the first hit).
 */
export function byName(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Open on a gene by name' });
  const a = assemblyA(b, files);
  const tx = b.data.ranges('transcripts', {
    target: a,
    source: fileSource(b, 'gff3-whole', 'gff3', files.gff3),
  } as never);
  b.views.view('main', (v) =>
    genome(v, 'toy-a', a, { include: ['chr1', 'chr2', 'chr3'] }, (s) =>
      s.track('transcripts', { label: 'Transcripts' }).layer('transcripts', {
        type: 'ranges',
        data: tx,
        line: 'gene',
      } as never),
    ),
  );
  return b.build();
}

/**
 * Genes of the GFF3 found by name (EC-228): the `gff` runtime's `findNames` over the file read whole (through
 * `assets`), exact names first, then by prefix. What the story's find box asks.
 */
export async function findGene(
  text: string,
  assets: AssetReader,
  signal: AbortSignal,
  files: Files = GenomesFiles,
): Promise<readonly NameHit[]> {
  const s = await GffRuntime.open(
    { address: '/data/gff3-whole', format: 'gff3', uri: files.gff3 },
    { assets, signal, limits: DefaultOpenLimits },
  );
  try {
    return (await s.findNames?.(text, { max: 8, prefix: true }, signal)) ?? [];
  } finally {
    s.dispose();
  }
}

/** The range a hit focuses: its item's extent on toy A's sequence, known from the hit without reading a tile. */
export function geneRange(hit: NameHit): string {
  return `/data/toy-a/(sequence :name ${JSON.stringify(hit.system)})/(range :from ${hit.start} :to ${hit.end})`;
}
