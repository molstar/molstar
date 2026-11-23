/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { NGap, type ToyAssembly } from './assemblies.js';
import { between, Concern, stream } from './streams.js';

// Step 3 of the Genomes generator (examples.md, Genomes; EC-227): 150 genes on toy A's chr1 to chr3 (none in chr1's run
// of N), three levels (gene, transcript, exon, and each transcript's CDS), both strands, named `TOY0001`…; a gene 3 to
// 12 exons, 1 to 4 transcripts (the first every exon, each other one skipping an inner exon). One model written as GFF3
// (`ID`, `Name`, `Parent`) and as GTF (`gene_id`, `transcript_id`, `gene_name`, `transcript_name`), each sorted by
// start: what tabix indexes (`genes.ts`'s lines by region read the same model as the whole file).

export interface ToyTranscript {
  readonly id: string;
  readonly name: string;
  /** Exons in genomic order, 0-based half-open. */
  readonly exons: readonly (readonly [number, number])[];
  /** The coding part (thick): from inside the first exon to inside the last. */
  readonly cds: readonly [number, number];
}

export interface ToyGene {
  readonly id: string;
  readonly name: string;
  readonly chr: string;
  readonly start: number;
  readonly end: number;
  readonly strand: '+' | '-';
  readonly transcripts: readonly ToyTranscript[];
}

/** The number of genes. */
export const GeneCount = 150;

/** The genes, in A's order and by start: slots along chr1 to chr3 (chr1's run of N skipped), a gene in each. */
export function toyGenes(a: ToyAssembly, seed: number, count = GeneCount): ToyGene[] {
  const chrs = a.sequences.filter((s) => /^chr[123]$/.test(s.name));
  const usable = chrs.reduce((n, s) => n + s.length, 0) - (NGap[1] - NGap[0]);
  const slot = Math.floor(usable / count);
  const genes: ToyGene[] = [];
  let k = 0;
  for (const c of chrs)
    for (let at = 0; at + slot <= c.length && k < count; at += slot) {
      if (c.name === 'chr1' && at + slot > NGap[0] && at < NGap[1]) continue;
      const rng = stream(seed, Concern.genes, k);
      const id = `TOY${String(k + 1).padStart(4, '0')}`;
      const length = between(rng, Math.floor(slot * 0.4), Math.floor(slot * 0.8));
      const start = at + between(rng, 100, slot - length - 100);
      const end = start + length;
      const strand = rng.float() < 0.5 ? '+' : '-';
      const n = between(rng, 3, 13);
      const step = length / n;
      const exons: [number, number][] = [];
      for (let e = 0; e < n; e++) {
        const s = start + Math.floor(e * step);
        exons.push([s, e === n - 1 ? end : s + Math.min(Math.floor(step * 0.7), between(rng, 60, 180))]);
      }
      const transcriptCount = between(rng, 1, 5);
      const transcripts: ToyTranscript[] = [];
      for (let t = 0; t < transcriptCount; t++) {
        // a transcript after the first skips one inner exon (alternative splicing)
        const skip = t > 0 && n > 2 ? 1 + ((t - 1) % (n - 2)) : -1;
        const ex = exons.filter((_, e) => e !== skip);
        const first = ex[0];
        const last = ex[ex.length - 1];
        transcripts.push({
          id: `${id}.${t + 1}`,
          name: `${id}-${201 + t}`,
          exons: ex,
          cds: [first[0] + Math.floor((first[1] - first[0]) / 3), last[1] - Math.floor((last[1] - last[0]) / 3)],
        });
      }
      genes.push({ id, name: id, chr: c.name, start, end, strand, transcripts });
      k++;
    }
  return genes;
}

/** A transcript's CDS pieces (its exons cut to its coding part), in genomic order. */
export function cdsPieces(t: ToyTranscript): [number, number][] {
  const out: [number, number][] = [];
  for (const [s, e] of t.exons) {
    const lo = Math.max(s, t.cds[0]);
    const hi = Math.min(e, t.cds[1]);
    if (hi > lo) out.push([lo, hi]);
  }
  return out;
}

/** The phase of each CDS piece (GFF3's column 8), counted along the transcript's strand. */
function phases(pieces: readonly (readonly [number, number])[], strand: '+' | '-'): number[] {
  const order = pieces.map((_, i) => i);
  if (strand === '-') order.reverse();
  const out = new Array<number>(pieces.length);
  let done = 0;
  for (const i of order) {
    out[i] = (3 - (done % 3)) % 3;
    done += pieces[i][1] - pieces[i][0];
  }
  return out;
}

/** A line of either format: its start for sorting and its rank among lines of one start (parents first). */
interface Line {
  readonly start: number;
  readonly rank: number;
  readonly text: string;
}

/** Each gene's lines sorted by start, parents before children at one start (stable otherwise). */
function sorted(lines: Line[]): string {
  return lines
    .map((l, i) => ({ l, i }))
    .sort((x, y) => x.l.start - y.l.start || x.l.rank - y.l.rank || x.i - y.i)
    .map((x) => `${x.l.text}\n`)
    .join('');
}

const cells = (...c: (string | number)[]) => c.join('\t');

/** The genes as GFF3: the header with each sequence's region, then each sequence's lines by start. */
export function genesGff3(a: ToyAssembly, genes: readonly ToyGene[]): string {
  let out = '##gff-version 3\n';
  for (const s of a.sequences) out += `##sequence-region ${s.name} 1 ${s.length}\n`;
  for (const c of a.sequences) {
    const lines: Line[] = [];
    for (const g of genes) {
      if (g.chr !== c.name) continue;
      const col = (type: string, s: number, e: number, phase: string, attrs: string) =>
        cells(g.chr, 'toy', type, s + 1, e, '.', g.strand, phase, attrs);
      lines.push({
        start: g.start,
        rank: 0,
        text: col('gene', g.start, g.end, '.', `ID=gene:${g.id};Name=${g.name};biotype=protein_coding`),
      });
      for (const t of g.transcripts) {
        const [s, e] = [t.exons[0][0], t.exons[t.exons.length - 1][1]];
        lines.push({
          start: s,
          rank: 1,
          text: col('mRNA', s, e, '.', `ID=transcript:${t.id};Parent=gene:${g.id};Name=${t.name}`),
        });
        for (const [xs, xe] of t.exons)
          lines.push({ start: xs, rank: 2, text: col('exon', xs, xe, '.', `Parent=transcript:${t.id}`) });
        const pieces = cdsPieces(t);
        const ph = phases(pieces, g.strand);
        pieces.forEach(([xs, xe], i) =>
          lines.push({
            start: xs,
            rank: 3,
            text: col('CDS', xs, xe, `${ph[i]}`, `ID=CDS:${t.id};Parent=transcript:${t.id}`),
          }),
        );
      }
    }
    out += sorted(lines);
  }
  return out;
}

/** The genes as GTF: each sequence's lines by start (GTF has no header the runtimes need). */
export function genesGtf(a: ToyAssembly, genes: readonly ToyGene[]): string {
  let out = '';
  for (const c of a.sequences) {
    const lines: Line[] = [];
    for (const g of genes) {
      if (g.chr !== c.name) continue;
      const col = (type: string, s: number, e: number, phase: string, attrs: string) =>
        cells(g.chr, 'toy', type, s + 1, e, '.', g.strand, phase, attrs);
      const gid = `gene_id "${g.id}";`;
      lines.push({ start: g.start, rank: 0, text: col('gene', g.start, g.end, '.', `${gid} gene_name "${g.name}";`) });
      for (const t of g.transcripts) {
        const tid = `${gid} transcript_id "${t.id}";`;
        const [s, e] = [t.exons[0][0], t.exons[t.exons.length - 1][1]];
        lines.push({
          start: s,
          rank: 1,
          text: col('transcript', s, e, '.', `${tid} gene_name "${g.name}"; transcript_name "${t.name}";`),
        });
        t.exons.forEach(([xs, xe], i) =>
          lines.push({ start: xs, rank: 2, text: col('exon', xs, xe, '.', `${tid} exon_number "${i + 1}";`) }),
        );
        const pieces = cdsPieces(t);
        const ph = phases(pieces, g.strand);
        pieces.forEach(([xs, xe], i) => lines.push({ start: xs, rank: 3, text: col('CDS', xs, xe, `${ph[i]}`, tid) }));
      }
    }
    out += sorted(lines);
  }
  return out;
}

/** A GFF or GTF line's sequence and extent for the index (columns 1, 4 and 5; 1-based inclusive). */
export const gffExtent = (c: readonly string[]) =>
  c.length >= 5 ? { seq: c[0], begin: Number(c[3]) - 1, end: Number(c[4]) } : undefined;
