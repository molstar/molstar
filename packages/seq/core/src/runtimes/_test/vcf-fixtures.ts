/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { bgzip, type IndexedLine, type IndexOptions, tabixIndex } from '@molstar/test-harness/genome/bgzf-writer';

// VCF files generated from a seed (core's PCG) for the VCF runtime's tests: SNVs (some multi-allelic), MNVs,
// insertions, deletions, symbolic deletions with `END=`, an `END=` on an SNV line, breakends, sequence-resolved SVs
// past 32 characters, a line written twice, ALTs that are no variant, and INFO keys of every `Number` and `Type`; each
// line's extent written by the generator from what it chose (not by the runtime's rule), bgzipped and tabix-indexed by
// the shared test-only writer.

export interface GeneratedVcf {
  readonly text: string;
  readonly lines: readonly IndexedLine[];
  readonly names: readonly string[];
  /** The lengths the header's `##contig` lines give (not every sequence has one). */
  readonly lengths: ReadonlyMap<string, number>;
}

export const VcfHeaderLines = [
  '##fileformat=VCFv4.2',
  '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
  '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, one an allele">',
  '##INFO=<ID=RC,Number=R,Type=Integer,Description="Read counts: the reference\'s, then each allele\'s">',
  '##INFO=<ID=DB,Number=0,Type=Flag,Description="In dbSNP">',
  '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Significance, \\"quoted\\", with commas">',
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The SV\'s type">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="The end">',
  '##INFO=<ID=GL,Number=G,Type=Float,Description="Per genotype">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
];

const enc = new TextEncoder();
const bases = 'ACGT';

function seqOf(rng: PCG, n: number): string {
  let s = '';
  for (let i = 0; i < n; i++) s += bases[Math.floor(rng.float() * 4)];
  return s;
}

/** A line's text and the extent its writer gives it (`POS − 1` to its END or REF's end). */
interface Row {
  seq: number;
  begin: number;
  end: number;
  text: string;
}

/**
 * `n` records over `names` (sorted), positions up to `span`; `contigs` names the sequences whose `##contig` says a
 * length (`span` + 1,000), and `samples` the sample columns (each `0|1`).
 */
export function generateVcf(
  n: number,
  names: readonly string[],
  seed = 1,
  span = 400_000,
  options: { readonly contigs?: readonly string[]; readonly samples?: number } = {},
): GeneratedVcf {
  const rng = new PCG(seed);
  const rows: Row[] = [];
  const samples = options.samples ?? 0;
  const gts = samples ? `\tGT${'\t0|1'.repeat(samples)}` : '';
  for (let i = 0; i < n; i++) {
    const seq = Math.floor(rng.float() * names.length);
    const pos = 1 + Math.floor(rng.float() * span);
    const begin = pos - 1;
    const kind = rng.float();
    const id = rng.float() < 0.6 ? `rs${i + 1}` : '.';
    const qual = rng.float() < 0.2 ? '.' : `${Math.floor(rng.float() * 1000) / 10}`;
    const filter = ['PASS', 'q10', '.'][Math.floor(rng.float() * 3)];
    const dp = rng.float() < 0.1 ? '' : `DP=${Math.floor(rng.float() * 200)};`;
    const db = rng.float() < 0.3 ? 'DB;' : '';
    const sig = rng.float() < 0.3 ? 'CLNSIG=Pathogenic,Likely_benign;' : '';
    let ref = seqOf(rng, 1);
    let alts: string[];
    let extra = '';
    let end = begin + 1;
    if (kind < 0.35) {
      // an SNV, some multi-allelic, some with an ALT that is no variant
      const count = rng.float() < 0.3 ? 2 + Math.floor(rng.float() * 2) : 1;
      alts = [];
      for (const b of bases) if (b !== ref && alts.length < count) alts.push(b);
      if (rng.float() < 0.1) alts.push(['*', '.', '<NON_REF>'][Math.floor(rng.float() * 3)]);
    } else if (kind < 0.45) {
      ref = seqOf(rng, 2 + Math.floor(rng.float() * 3));
      alts = [seqOf(rng, ref.length)];
      end = begin + ref.length;
    } else if (kind < 0.55) {
      alts = [ref + seqOf(rng, 1 + Math.floor(rng.float() * 6))];
    } else if (kind < 0.65) {
      ref = seqOf(rng, 2 + Math.floor(rng.float() * 8));
      alts = [ref[0]];
      end = begin + ref.length;
    } else if (kind < 0.72) {
      // a symbolic deletion to its END (EC-222: a long deletion across the SNVs it covers)
      const length = 50 + Math.floor(rng.float() * 30_000);
      alts = ['<DEL>'];
      end = pos + length;
      extra = `SVTYPE=DEL;END=${end};`;
    } else if (kind < 0.76) {
      // END= on an SNV line: the index takes it whatever the alleles
      alts = [bases[(bases.indexOf(ref) + 1) % 4]];
      end = pos + 100;
      extra = `END=${end};`;
    } else if (kind < 0.8) {
      // END= not past POS − 1: not the end (REF's)
      alts = [bases[(bases.indexOf(ref) + 2) % 4]];
      extra = `END=${pos - 1};`;
    } else if (kind < 0.84) {
      // a breakend (its mate elsewhere)
      alts = [`${ref}]${names[seq]}:${pos + 500}]`];
      extra = 'SVTYPE=BND;';
    } else if (kind < 0.88) {
      // a sequence-resolved deletion past 32 characters: its id hashed
      ref = seqOf(rng, 40 + Math.floor(rng.float() * 60));
      alts = [ref[0]];
      end = begin + ref.length;
    } else {
      alts = [bases[(bases.indexOf(ref) + 3) % 4]];
    }
    const af = alts.map(() => (Math.floor(rng.float() * 1000) / 1000).toString()).join(',');
    const rc = [ref, ...alts].map(() => Math.floor(rng.float() * 50)).join(',');
    const gl = rng.float() < 0.2 ? 'GL=-0.1,-1.5,-3;' : '';
    const info = `${dp}AF=${af};RC=${rc};${db}${sig}${gl}${extra}`.replace(/;$/, '') || '.';
    const text = `${names[seq]}\t${pos}\t${id}\t${ref}\t${alts.join(',')}\t${qual}\t${filter}\t${info}${gts}`;
    rows.push({ seq, begin, end, text });
    // a line written twice (a duplicated record: `#2` on its ids)
    if (rng.float() < 0.03) rows.push({ seq, begin, end, text });
  }
  rows.sort((a, b) => a.seq - b.seq || a.begin - b.begin);
  const lengths = new Map<string, number>();
  const header = [...VcfHeaderLines];
  for (const name of options.contigs ?? []) {
    lengths.set(name, span + 1000);
    header.splice(1, 0, `##contig=<ID=${name},length=${span + 1000}>`);
  }
  header.push(
    `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO${samples ? `\tFORMAT${Array.from({ length: samples }, (_, k) => `\tS${k + 1}`).join('')}` : ''}`,
  );
  let at = 0;
  const parts: string[] = [];
  for (const h of header) {
    parts.push(`${h}\n`);
    at += enc.encode(`${h}\n`).length;
  }
  const lines: IndexedLine[] = [];
  for (const r of rows) {
    const t = `${r.text}\n`;
    const b = enc.encode(t).length;
    lines.push({ seq: r.seq, begin: r.begin, end: r.end, start: at, stop: at + b });
    parts.push(t);
    at += b;
  }
  return { text: parts.join(''), lines, names, lengths };
}

/** VCF's index options (1-based, columns 1, 2, 0, `#` meta) over `names`. */
export const vcfIndex = (names: readonly string[]): IndexOptions => ({
  kind: 'tbi',
  format: 2,
  zeroBased: false,
  columns: [1, 2, 0],
  meta: '#',
  skip: 0,
  names,
});

/** A generated file bgzipped (blocks of `blockSize`) with its tabix index. */
export function indexedVcf(g: GeneratedVcf, blockSize = 4000) {
  const file = bgzip(enc.encode(g.text), blockSize);
  return { file: file.bytes, index: tabixIndex(g.lines, file, vcfIndex(g.names)) };
}

/** The header lines of `VcfHeaderLines` with a `#CHROM` line (sites only), as text. */
export const VcfHeaderText = `${[...VcfHeaderLines, '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'].join('\n')}\n`;

/** `n` data lines of 199 bytes without a tab (each fewer than eight columns), and no tab after them either. */
export function tablessLines(n: number): Uint8Array {
  const line = `chr1 100 . A C . . ${'x'.repeat(180)}`.slice(0, 199);
  return enc.encode(`${Array.from({ length: n }, () => line).join('\n')}\n`);
}

/** One line whose INFO is `n` flags (no `=` in any). */
export function flagsLine(n: number): Uint8Array {
  const flags = Array.from({ length: n }, (_, k) => `F${k}`).join(';');
  return enc.encode(`chr1\t100\t.\tA\tC\t.\t.\tDB;${flags}\n`);
}
