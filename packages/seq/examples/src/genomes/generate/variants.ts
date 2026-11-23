/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { NGap, type ToyAssembly } from './assemblies.js';
import { Concern, stream } from './streams.js';

// Step 5 of the Genomes generator (examples.md, Genomes; EC-220, EC-222): a VCF of 40 samples over toy A, `GT` alone:
// 2,000 sites (SNVs at A's letters, a tenth multi-allelic), a third of the lines' diploid calls phased, a
// partly called one now and then, and no-calls; chrM's calls haploid (every sample), and one 600-base deletion on chr1
// (sequence-resolved, its REF A's letters) over 12 SNVs, carried by six samples. Sorted: what tabix indexes. Why a
// file of its own beside `generated-genotypes.vcf` (the VCF feature's 69 sites on `chrG`, kept byte for byte): the
// dataset's views read toy A's sequences by region, at a tile of sites.

export const SampleCount = 40;
export const SiteCount = 2_000;
/** The long deletion: its REF's extent on chr1 (0-based, half-open). */
export const DeletionExtent: readonly [number, number] = [20_000, 20_600];
/** The samples (0-based) carrying the deletion: the first homozygous. */
export const DeletionCarriers: readonly number[] = [3, 8, 14, 22, 29, 35];
/** SNVs inside the deletion's extent. */
export const SnvsUnderDeletion = 12;

export interface VariantsTruth {
  readonly text: string;
  /** Sites a sequence. */
  readonly sites: ReadonlyMap<string, number>;
}

export const sampleNames = (n = SampleCount) =>
  Array.from({ length: n }, (_, s) => `S${String(s + 1).padStart(2, '0')}`);

/** The VCF, its sites counted by sequence. */
export function variantsVcf(
  a: ToyAssembly,
  letters: ReadonlyMap<string, string>,
  seed: number,
  samples = SampleCount,
  siteCount = SiteCount,
): VariantsTruth {
  const names = sampleNames(samples);
  const S = names.length;
  const seqs = a.sequences.filter((s) => /^chr[123M]$/.test(s.name));
  const deletionSnvs = SnvsUnderDeletion;
  const spread = siteCount - deletionSnvs - 1;
  const total = seqs.reduce((n, s) => n + s.length, 0);
  const rows: string[] = [];
  const sites = new Map<string, number>();
  let k = 0;
  const snv = (chr: string, pos: number, haploid: boolean) => {
    const rng = stream(seed, Concern.variants, k++);
    const ref = letters.get(chr)![pos];
    const others = 'ACGT'.replace(ref, '');
    const a1 = others[Math.floor(rng.float() * 3)];
    const multi = rng.float() < 0.1;
    const alt = multi ? `${a1},${others.replace(a1, '')[Math.floor(rng.float() * 2)]}` : a1;
    const af = 0.03 + 0.45 * rng.float();
    const allele = () => (rng.float() < af ? (multi && rng.float() < 0.3 ? 2 : 1) : 0);
    // a line's diploid calls phased, or not (a third of the lines phased)
    const sep = rng.float() < 0.3 ? '|' : '/';
    const calls: string[] = [];
    for (let s = 0; s < S; s++) {
      const r = rng.float();
      if (haploid) calls.push(r < 0.03 ? '.' : `${allele()}`);
      else if (r < 0.02) calls.push('./.');
      else if (r < 0.03) calls.push(`./${allele()}`);
      else calls.push(`${allele()}${sep}${allele()}`);
    }
    rows.push([chr, pos + 1, `tv${k}`, ref, alt, '.', 'PASS', '.', 'GT', ...calls].join('\t'));
    sites.set(chr, (sites.get(chr) ?? 0) + 1);
  };
  for (const s of seqs) {
    // the last sequence takes what the others left: 2,000 sites in all
    const n = s === seqs[seqs.length - 1] ? siteCount - rows.length : Math.round((spread * s.length) / total);
    const step = s.length / n;
    for (let i = 0; i < n; i++) {
      let pos = Math.floor(i * step + step / 2);
      if (s.name === 'chr1' && pos >= NGap[0] && pos < NGap[1]) pos = NGap[1] + (pos - NGap[0]) / 4;
      if (s.name === 'chr1' && pos >= DeletionExtent[0] - 50 && pos < DeletionExtent[1] + 50) continue;
      snv(s.name, Math.floor(pos), s.name === 'chrM');
    }
    if (s.name === 'chr1') {
      // the deletion, then its SNVs a 50 bp step apart (a site's REF before its SNVs at one position: by POS)
      const [d0, d1] = DeletionExtent;
      const ref = letters.get('chr1')!.slice(d0, d1 + 1);
      const calls = names.map((_, i) =>
        DeletionCarriers.includes(i) ? (i === DeletionCarriers[0] ? '1/1' : '0|1') : '0/0',
      );
      rows.push(['chr1', d0 + 1, 'del600', ref, ref[0], '.', 'PASS', '.', 'GT', ...calls].join('\t'));
      sites.set('chr1', (sites.get('chr1') ?? 0) + 1);
      for (let j = 0; j < deletionSnvs; j++) snv('chr1', d0 + 25 + j * 50, false);
    }
  }
  // each sequence's rows by position (the deletion's 12 SNVs were written after its other sites)
  const pos = (r: string) => Number(r.split('\t', 2)[1]);
  const order = new Map(a.sequences.map((s, i) => [s.name, i]));
  const chrOf = (r: string) => order.get(r.slice(0, r.indexOf('\t')))!;
  rows.sort((x, y) => chrOf(x) - chrOf(y) || pos(x) - pos(y));
  const header = [
    '##fileformat=VCFv4.2',
    ...a.sequences.map((s) => `##contig=<ID=${s.name},length=${s.length}>`),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${names.join('\t')}`,
  ];
  return { text: `${[...header, ...rows].join('\n')}\n`, sites };
}

/** A VCF line's sequence and extent for the index (`POS − 1` to REF's end). */
export const vcfExtent = (c: readonly string[]) =>
  c.length >= 4 ? { seq: c[0], begin: Number(c[1]) - 1, end: Number(c[1]) - 1 + c[3].length } : undefined;

/** Sites of the deletion's SNVs (1-based POS), for the tests. */
export const deletionSnvPositions = () =>
  Array.from({ length: SnvsUnderDeletion }, (_, j) => DeletionExtent[0] + 25 + j * 50 + 1);

/** The sites' total. */
export const sitesTotal = (t: VariantsTruth) => [...t.sites.values()].reduce((n, v) => n + v, 0);
