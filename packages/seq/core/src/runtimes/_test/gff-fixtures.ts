/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { bgzip, type IndexedLine, tabixIndex } from '@molstar/test-harness/genome/bgzf-writer';

// GFF3 and GTF written by a generator (the runtime's tests read each): genes of one to three transcripts, exons with a
// CDS between UTRs on both strands, sorted by position, over a few sequences with their `##sequence-region`; and the
// recorded real slices (TP53's region from Ensembl's GFF3 and GENCODE's GTF) as bytes.

export interface GeneratedGff {
  readonly text: string;
  /** Lines of data (not directives). */
  readonly lines: number;
  readonly genes: number;
  readonly transcripts: number;
  readonly lengths: ReadonlyMap<string, number>;
}

interface Tx {
  id: string;
  start: number;
  end: number;
  exons: [number, number][];
  cds?: [number, number];
}

/** One generated gene's transcripts: exons from a shared set, a CDS inside the first and last exon chosen. */
function geneModel(rng: PCG, gid: string, start: number, count: number): { end: number; txs: Tx[] } {
  const exonCount = 2 + Math.floor(rng.float() * 6);
  const exons: [number, number][] = [];
  let at = start;
  for (let e = 0; e < exonCount; e++) {
    const len = 60 + Math.floor(rng.float() * 300);
    exons.push([at, at + len]);
    at += len + 200 + Math.floor(rng.float() * 2000);
  }
  const txs: Tx[] = [];
  for (let t = 0; t < count; t++) {
    // a transcript skips one inner exon now and then (a shared exon is a piece of each)
    const skip = exonCount > 2 && rng.float() < 0.5 ? 1 + Math.floor(rng.float() * (exonCount - 2)) : -1;
    const ex = exons.filter((_, i) => i !== skip);
    const coding = rng.float() < 0.8;
    const first = ex[0];
    const last = ex[ex.length - 1];
    const cds: [number, number] | undefined = coding
      ? [first[0] + Math.floor((first[1] - first[0]) / 3), last[1] - Math.floor((last[1] - last[0]) / 3)]
      : undefined;
    txs.push({ id: `${gid}.${t + 1}`, start: first[0], end: last[1], exons: ex, ...(cds ? { cds } : {}) });
  }
  return { end: exons[exonCount - 1][1], txs };
}

/** The pieces of a transcript: its exons, and its CDS and UTRs cut by the CDS. */
function parts(tx: Tx): { type: string; s: number; e: number; phase?: number }[] {
  const out: { type: string; s: number; e: number; phase?: number }[] = [];
  for (const [s, e] of tx.exons) out.push({ type: 'exon', s, e });
  if (!tx.cds) return out;
  const [cs, ce] = tx.cds;
  let done = 0;
  for (const [s, e] of tx.exons) {
    if (e <= cs || s >= ce) {
      out.push({ type: s < cs ? 'five_prime_UTR' : 'three_prime_UTR', s, e });
      continue;
    }
    if (s < cs) out.push({ type: 'five_prime_UTR', s, e: cs });
    const a = Math.max(s, cs);
    const b = Math.min(e, ce);
    out.push({ type: 'CDS', s: a, e: b, phase: (3 - (done % 3)) % 3 });
    done += b - a;
    if (e > ce) out.push({ type: 'three_prime_UTR', s: ce, e });
  }
  return out;
}

/**
 * A GFF3 of `genes` genes over `names` (positions to `span`): gene → mRNA or lnc_RNA → exon, CDS, UTRs; IDs
 * `gene:G1`, `transcript:G1.1`; `Name`, an `Alias` on a gene now and then.
 */
export function generateGff3(genes: number, names: readonly string[], seed = 1, span = 2_000_000): GeneratedGff {
  return generate('gff3', genes, names, seed, span, true);
}

/** A GTF of the same genes: with `transcript` lines (GENCODE's) or without (GTF 2.2's). */
export function generateGtf(
  genes: number,
  names: readonly string[],
  seed = 1,
  transcriptLines = true,
  span = 2_000_000,
): GeneratedGff {
  return generate('gtf', genes, names, seed, span, transcriptLines);
}

function generate(
  format: 'gff3' | 'gtf',
  genes: number,
  names: readonly string[],
  seed: number,
  span: number,
  transcriptLines: boolean,
): GeneratedGff {
  const rng = new PCG(seed);
  const head: string[] = format === 'gff3' ? ['##gff-version 3'] : ['##format: gtf'];
  const lengths = new Map<string, number>();
  for (const n of names) {
    lengths.set(n, span + 50_000);
    if (format === 'gff3') head.push(`##sequence-region ${n} 1 ${span + 50_000}`);
  }
  const body: { seq: number; start: number; text: string }[] = [];
  let transcripts = 0;
  const per = Math.max(1, Math.ceil(genes / names.length));
  let g = 0;
  for (let si = 0; si < names.length && g < genes; si++) {
    let at = 1000;
    for (let k = 0; k < per && g < genes; k++, g++) {
      const gid = `G${g + 1}`;
      const count = 1 + Math.floor(rng.float() * 3);
      const { end, txs } = geneModel(rng, gid, at, count);
      transcripts += txs.length;
      const strand = rng.float() < 0.5 ? '+' : '-';
      const seq = names[si];
      const line = (type: string, s: number, e: number, attrs: string, phase = '.') =>
        body.push({ seq: si, start: s, text: [seq, 'gen', type, s + 1, e, '.', strand, phase, attrs].join('\t') });
      if (format === 'gff3') {
        const alias = rng.float() < 0.3 ? `;Alias=ALT${g + 1},OLD${g + 1}` : '';
        line('gene', at, end, `ID=gene:${gid};Name=GENE${g + 1}${alias};biotype=protein_coding`);
        for (const tx of txs) {
          line(
            tx.cds ? 'mRNA' : 'lnc_RNA',
            tx.start,
            tx.end,
            `ID=transcript:${tx.id};Parent=gene:${gid};Name=GENE${g + 1}-${tx.id.split('.')[1]}`,
          );
          for (const p of parts(tx))
            line(p.type, p.s, p.e, `Parent=transcript:${tx.id}`, p.phase === undefined ? '.' : `${p.phase}`);
        }
      } else {
        const ga = `gene_id "${gid}"; gene_name "GENE${g + 1}";`;
        line('gene', at, end, ga);
        for (const tx of txs) {
          const ta = `gene_id "${gid}"; transcript_id "${tx.id}"; gene_name "GENE${g + 1}"; transcript_name "GENE${g + 1}-${tx.id.split('.')[1]}"; level 2;`;
          if (transcriptLines) line('transcript', tx.start, tx.end, ta);
          let n = 0;
          for (const p of parts(tx))
            line(
              p.type === 'five_prime_UTR' || p.type === 'three_prime_UTR' ? 'UTR' : p.type,
              p.s,
              p.e,
              `${ta} exon_number ${++n}; tag "basic"; tag "CCDS";`,
              p.phase === undefined ? '.' : `${p.phase}`,
            );
        }
      }
      at = end + 500 + Math.floor(rng.float() * 5000);
      if (at > span - 40_000) break;
    }
  }
  body.sort((a, b) => a.seq - b.seq || a.start - b.start);
  return {
    text: `${[...head, ...body.map((b) => b.text)].join('\n')}\n`,
    lines: body.length,
    genes: g,
    transcripts,
    lengths,
  };
}

/** A recorded real slice's bytes (examples' `fixtures/genomes/<name>.pages.json`, every page recorded). */
export function recordedSlice(name: string): Uint8Array {
  const path = fileURLToPath(new URL(`../../../../examples/fixtures/genomes/${name}.pages.json`, import.meta.url));
  const r = JSON.parse(readFileSync(path, 'utf8')) as { size: number; pageSize: number; pages: Record<string, string> };
  const out = new Uint8Array(r.size);
  for (const [k, b64] of Object.entries(r.pages)) out.set(Buffer.from(b64, 'base64'), Number(k) * r.pageSize);
  return out;
}

export const EnsemblTp53 = 'Homo_sapiens.GRCh38.116.chr.TP53.gff3.gz';
export const GencodeTp53 = 'gencode.v50.annotation.TP53.gtf.gz';

/**
 * A GFF3 or GTF text bgzipped (blocks of `blockSize`) with its tabix index (GFF's preset: columns 1, 4, 5, 1-based,
 * `#` meta) over `names`, each line's extent its own columns.
 */
export function indexedGff(text: string, names: readonly string[], blockSize = 4000) {
  const enc = new TextEncoder();
  const bytes = enc.encode(text);
  const lines: IndexedLine[] = [];
  for (let at = 0; at < bytes.length; ) {
    let e = bytes.indexOf(10, at);
    if (e < 0) e = bytes.length;
    const stop = Math.min(bytes.length, e + 1);
    if (bytes[at] !== 35 && e > at) {
      const cells = new TextDecoder().decode(bytes.subarray(at, e)).split('\t');
      lines.push({ seq: names.indexOf(cells[0]), begin: Number(cells[3]) - 1, end: Number(cells[4]), start: at, stop });
    }
    at = stop;
  }
  const file = bgzip(bytes, blockSize);
  const index = tabixIndex(lines, file, {
    kind: 'tbi',
    format: 0,
    zeroBased: false,
    columns: [1, 4, 5],
    meta: '#',
    skip: 0,
    names,
  });
  return { file: file.bytes, index };
}
