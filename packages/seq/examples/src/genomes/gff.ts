/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import {
  type AssetReader,
  DefaultOpenLimits,
  delayedRuntime,
  type NameHit,
  type StoreClock,
} from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { BedFormats, BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { GffFormats, GffRuntime } from '@molstar/seqx-core/runtimes/gff';
import { recordedAssets, withRecordedPages } from '../recorded-pages.js';
import { EnsemblTp53File, GencodeTp53File, MisreadGffFile, RecordedFiles, Tp53Region } from './files.js';

// GFF3's and GTF's views (Features/Genome formats/GFF3 and GTF): TP53's region from Ensembl's GFF3 and GENCODE's GTF
// (the recorded real slices) read by the opt-in `gff` runtime, its transcripts as gene models (the ranges layer's gene
// line) and, from the same file with `features: genes`, one model a gene; a misread GFF3 and its issues. GRCh38's chr17
// is `17` in Ensembl's file (its length from the file's `##sequence-region`) and `chr17` in GENCODE's (from an inline
// chrom.sizes: a GTF says no length).

/** The sections' addresses (what the tests read). */
export const Tp53Section = '/views/main/chr17';
export const MisreadGffSection = '/views/main/misread';
/** The transcripts' annotation of each TP53 view (an item's address is `<it>/(item :id …)`). */
export const TranscriptsAddress = '/data/transcripts';

/** chr17's length in GRCh38. */
const Chr17 = 83_257_441;

/** TP53's region of a file: its transcripts as gene models, then its genes (`features: genes`), labelled by name. */
function tp53(format: 'gff3' | 'gtf'): Node {
  const gff = format === 'gff3';
  const b = createSeqSpec({ title: gff ? 'TP53 from Ensembl’s GFF3' : 'TP53 from GENCODE’s GTF' });
  const common = {
    format,
    uri: gff ? EnsemblTp53File : GencodeTp53File,
    attribution: gff ? 'Ensembl release 116 (GRCh38)' : 'GENCODE release 50 (GRCh38)',
    licence: 'LicenseRef-EMBL-EBI-terms',
  };
  const file = b.data.source('tp53-file', common);
  const genesFile = b.data.source('tp53-genes', { ...common, features: 'genes' } as never);
  const name = gff ? '17' : 'chr17';
  // GFF3's chromosome length is its `##sequence-region`'s; a GTF's from an inline chrom.sizes
  const lengths = gff ? file : b.data.source('sizes', { format: 'chrom-sizes', value: `chr17\t${Chr17}\n` as never });
  const chr = b.data.sequence('chr17', { alphabet: 'dna', source: lengths, name });
  const transcripts = b.data.ranges('transcripts', { target: chr, source: file });
  const genes = b.data.ranges('genes', { target: chr, source: genesFile });
  b.views.view('main', (v) => {
    const s = v.section('chr17', { of: chr, range: [...Tp53Region] });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    // a lane a transcript: TP53's thirty-odd overlap one another (a tiled source's lanes are not known before it is read)
    s.track('genes', { label: 'Genes', height: 48 }).layer('genes', {
      type: 'ranges',
      data: genes,
      line: 'gene',
    } as never);
    s.track('transcripts', { label: 'Transcripts', height: 480 }).layer('transcripts', {
      type: 'ranges',
      data: transcripts,
      line: 'gene',
    } as never);
  });
  return b.build();
}

/** TP53's region from Ensembl GRCh38 release 116's GFF3 (the recorded slice). */
export function tp53Gff3(): Node {
  return tp53('gff3');
}

/** TP53's region from GENCODE v50's GTF (the recorded slice). */
export function tp53Gtf(): Node {
  return tp53('gtf');
}

/** A misread GFF3: what reads drawn as gene models, the rest said (a line number each) under the view. */
export function misread(): Node {
  const b = createSeqSpec({ title: 'A misread GFF3' });
  const file = b.data.source('misread-gff', { format: 'gff3', uri: MisreadGffFile });
  const chr = b.data.sequence('chrG', { alphabet: 'dna', source: file, name: 'chrG' });
  const rows = b.data.ranges('rows', { target: chr, source: file });
  b.views.view('main', (v) => {
    const s = v.section('misread', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('rows', { label: 'Read' }).layer('ranges', { type: 'ranges', data: rows, line: 'gene' } as never);
  });
  return b.build();
}

/**
 * Items of a TP53 view's file found by name (EC-228): the GFF runtime's `findNames` over the file (read as the view
 * reads it, from its recording through `assets`), exact names first, then by prefix. What a locus box asks; the story
 * shows the hits and focuses the first.
 */
export async function findTp53(
  format: 'gff3' | 'gtf',
  text: string,
  assets: AssetReader,
  signal: AbortSignal,
): Promise<readonly NameHit[]> {
  const s = await GffRuntime.open(
    { address: '/data/tp53-file', format, uri: format === 'gff3' ? EnsemblTp53File : GencodeTp53File },
    { assets: recordedAssets(assets, RecordedFiles), signal, limits: DefaultOpenLimits },
  );
  try {
    return (await s.findNames?.(text, { max: 8, prefix: true }, signal)) ?? [];
  } finally {
    s.dispose();
  }
}

/** A name's hit's item address (its transcript's, when it has an id), what a selection or a bookmark names. */
export function hitAddress(hit: NameHit): string | undefined {
  return hit.id !== undefined ? `${TranscriptsAddress}/(item :id ${JSON.stringify(hit.id)})` : undefined;
}

/** The range a name's hit focuses on chr17: the item's whole extent, known from the hit without reading a tile. */
export function hitRange(hit: NameHit): string {
  return `/data/chr17/(range :from ${hit.start} :to ${hit.end})`;
}

/**
 * What the views need besides the default registry: the `gff` runtime (and the `bed` one for the GTF view's
 * chrom.sizes), their files read from their recordings where they have one, each answer after `delayMs` (150) on `clock`
 * (the synthetic twin the stories and tests share).
 */
export function registry(options: { readonly delayMs?: number; readonly clock?: StoreClock } = {}): SeqRegistryEntry[] {
  const timing = { delayMs: options.delayMs ?? 150, ...(options.clock ? { clock: options.clock } : {}) };
  return [
    {
      kind: 'format-runtime',
      runtime: delayedRuntime(withRecordedPages(GffRuntime, RecordedFiles), {
        name: 'gff-async',
        formats: GffFormats,
        ...timing,
      }),
    },
    {
      kind: 'format-runtime',
      runtime: delayedRuntime(BedRuntime, { name: 'bed-async', formats: BedFormats, ...timing }),
    },
  ];
}
