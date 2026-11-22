/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { type AssetReader, DefaultOpenLimits, delayedRuntime, type StoreClock } from '@molstar/seqx-core/data';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { BedFormats, BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import {
  VcfCategories,
  type VcfHeader,
  VcfFormats,
  VcfRuntime,
  type VcfSource,
  vcfRuntime,
} from '@molstar/seqx-core/runtimes/vcf';
import { recordedAssets, withRecordedPages } from '../recorded-pages.js';
import {
  Chr17Length,
  ClinvarFile,
  ClinvarIndexFile,
  GeneratedGenotypesFile,
  GeneratedSystem,
  GeneratedVcfFile,
  MisreadVcfFile,
  RecordedFiles,
  Tp53Region,
} from './files.js';

// VCF's views (Features/Genome formats/VCF): sites read by the opt-in `vcf` runtime, a row an ALT allele drawn by the
// ranges layer as built, a mark at its REF extent (a long deletion a bar to its END) coloured by `category`. TP53's
// variants from NCBI ClinVar (the recorded slice, read by region through its named index), the same file's chr17 at
// the fit (past the view's tiles: the band waiting, nothing read), a generated file with long deletions over SNVs, and a
// misread file and its issues. ClinVar names GRCh38's chr17 `17`: the sequence here is named as the file names it (an
// assembly would map `17` onto `chr17` by the names rule's `chr` prefix), its length from an inline chrom.sizes (ClinVar
// writes no `##contig`).

/** The sections' addresses (what the tests read). */
export const ClinvarSection = '/views/main/chr17';
export const GeneratedVcfSection = '/views/main/variants';
export const GenotypesSection = '/views/main/genotypes';
export const MisreadVcfSection = '/views/main/misread';
/** The variants' annotation of each view (an allele's address is `<it>/(item :id "7675088:C>T")`). */
export const VariantsAddress = '/data/variants';

/** An allele's colour by its category: the next categorical colours, every category listed (a complete key). */
const byCategory = () => theme.field('category', { palette: theme.categorical({ values: [...VcfCategories] }) });

/** ClinVar's source: the file and the index it names (EC-186), its attribution and licence. */
const clinvar = {
  format: 'vcf',
  uri: ClinvarFile,
  index: ClinvarIndexFile,
  attribution: 'NCBI ClinVar, GRCh38 VCF of 2024-01-07',
  licence: 'LicenseRef-NCBI-public-data',
} as const;

/** ClinVar's chr17: TP53's region (`range`), or the whole chromosome at the fit. */
function clinvarView(title: string, range?: readonly [number, number]): Node {
  const b = createSeqSpec({ title });
  const sizes = b.data.source('sizes', { format: 'chrom-sizes', value: `17\t${Chr17Length}\n` as never });
  const file = b.data.source('clinvar', clinvar as never);
  const chr = b.data.sequence('chr17', { alphabet: 'dna', source: sizes, name: '17' });
  const variants = b.data.ranges('variants', { target: chr, source: file });
  b.views.view('main', (v) => {
    const s = v.section('chr17', range ? { of: chr, range: [...range] } : { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('variants', { label: 'ClinVar', height: 160 }).layer('variants', {
      type: 'ranges',
      data: variants,
      color: byCategory(),
    });
  });
  return b.build();
}

/** TP53's variants in ClinVar (the recorded slice, by region): an allele a row, coloured by its category. */
export function tp53Clinvar(): Node {
  return clinvarView('TP53’s variants from ClinVar', Tp53Region);
}

/** ClinVar's chr17 at the fit: more tiles than a view reads, so the track stays the waiting band and nothing is read. */
export function clinvarZoomIn(): Node {
  return clinvarView('ClinVar’s chr17 at the fit');
}

/** A generated file (read whole): SNVs, indels and two long deletions, each one item across the SNVs it covers. */
export function longDeletion(): Node {
  const b = createSeqSpec({ title: 'A long deletion over SNVs' });
  const file = b.data.source('generated-vcf', { format: 'vcf', uri: GeneratedVcfFile });
  const chr = b.data.sequence('chrG', { alphabet: 'dna', source: file, name: GeneratedSystem });
  const variants = b.data.ranges('variants', { target: chr, source: file });
  b.views.view('main', (v) => {
    const s = v.section('variants', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('variants', { label: 'Variants', height: 72 }).layer('variants', {
      type: 'ranges',
      data: variants,
      color: byCategory(),
    });
  });
  return b.build();
}

/**
 * Generated genotypes (read whole): the sites above, the matrix of 40 samples' genotype classes below, a row a sample in
 * the header's order; the long deletion one item across the SNVs it covers, solid in its carriers' rows (EC-222),
 * uncalled calls in the missing colour (EC-220), the haploid samples' calls every copy or none.
 */
export function genotypes(): Node {
  const b = createSeqSpec({ title: 'Genotypes of 40 samples' });
  const file = b.data.source('genotypes', { format: 'vcf', uri: GeneratedGenotypesFile });
  const chr = b.data.sequence('chrG', { alphabet: 'dna', source: file, name: GeneratedSystem });
  const sites = b.data.ranges('variants', { target: chr, source: file });
  const gt = b.data.matrix('gt', { target: chr, source: file } as never);
  b.views.view('main', (v) => {
    const s = v.section('genotypes', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('variants', { label: 'Sites', height: 48 }).layer('variants', {
      type: 'ranges',
      data: sites,
      color: byCategory(),
    });
    s.track('gt', { label: 'Genotypes', height: 260 }).layer('gt', { type: 'matrix', data: gt });
  });
  return b.build();
}

/** A misread VCF: the lines read drawn, the rest said (a line number each) under the view. */
export function misread(): Node {
  const b = createSeqSpec({ title: 'A misread VCF' });
  const file = b.data.source('misread-vcf', { format: 'vcf', uri: MisreadVcfFile });
  const chr = b.data.sequence('chrG', { alphabet: 'dna', source: file, name: GeneratedSystem });
  const rows = b.data.ranges('variants', { target: chr, source: file });
  b.views.view('main', (v) => {
    const s = v.section('misread', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('rows', { label: 'Read' }).layer('variants', { type: 'ranges', data: rows, color: byCategory() });
  });
  return b.build();
}

/** The `vcf` runtime over the recordings: level 1 only, as recorded (TP53's pages). */
const RecordedVcfRuntime = vcfRuntime({ maxLevel: 1 });

/** A header as a few lines: its format, contigs, INFO and FORMAT keys (`Number`, `Type`) and samples. */
export function headerSummary(h: VcfHeader): string {
  const keys = (list: VcfHeader['info']) => list.map((d) => `  ${d.id} (Number=${d.number}, ${d.type})`);
  return [
    `${h.fileformat ?? 'VCF (no ##fileformat)'}: ${h.lines} header lines, ${h.bytes} bytes`,
    `Contigs: ${h.contigs.length ? h.contigs.map((c) => (c.length ? `${c.name} (${c.length})` : c.name)).join(', ') : 'none (the names are the index’s or the lines’)'}`,
    `INFO keys: ${h.info.length}`,
    ...keys(h.info),
    `FORMAT keys: ${h.format.length}${h.format.length ? '' : ' (sites only)'}`,
    ...keys(h.format),
    `Samples: ${h.samples.length}`,
  ].join('\n');
}

/**
 * A VCF file's header summary, read as a view reads the file (by region through `index`, from its recording where it
 * has one), the source opened by the `vcf` runtime: what the story prints in place of the file's head.
 */
export async function vcfHeaderSummary(
  uri: string,
  index: string | undefined,
  assets: AssetReader,
  signal: AbortSignal,
): Promise<string> {
  const s = (await VcfRuntime.open(
    { address: '/data/file', format: 'vcf', uri, ...(index !== undefined && { index }) },
    { assets: recordedAssets(assets, RecordedFiles), signal, limits: DefaultOpenLimits },
  )) as VcfSource;
  try {
    return headerSummary(s.header);
  } finally {
    s.dispose();
  }
}

/**
 * What the views need besides the default registry: the `vcf` runtime (and the `bed` one for ClinVar's chrom.sizes),
 * their files read from their recordings where they have one, each answer after `delayMs` (150) on `clock` (the
 * synthetic twin the stories and tests share).
 */
export function registry(options: { readonly delayMs?: number; readonly clock?: StoreClock } = {}): SeqRegistryEntry[] {
  const timing = { delayMs: options.delayMs ?? 150, ...(options.clock ? { clock: options.clock } : {}) };
  return [
    {
      kind: 'format-runtime',
      // level 1 alone (`maxLevel`): ClinVar's index lists 1 to 4,096, but its recording holds TP53's pages only, so the
      // coarser tiles of chr17 at the fit (each within the view's read cost) could not be read from it
      runtime: delayedRuntime(withRecordedPages(RecordedVcfRuntime, RecordedFiles), {
        name: 'vcf-async',
        formats: VcfFormats,
        ...timing,
      }),
    },
    {
      kind: 'format-runtime',
      runtime: delayedRuntime(BedRuntime, { name: 'bed-async', formats: BedFormats, ...timing }),
    },
  ];
}
