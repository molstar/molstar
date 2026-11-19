/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { type DynamicNode, encodeSeq, type SeqData, seqToValue } from '@molstar/seqx-core/data';
import { itemAddress, nodeAddress } from '@molstar/seqx-core/resolve';
import { Domains, GeneExons, GeneLength, geneLetters, seeded } from './data.js';

// The p53 gene's region read by region: a toy stretch of a chromosome with the toy TP53 gene in it (its letters at
// `GeneOffset`) and a neighbor on the other strand, served in regions of `RegionSize` by a `seq` source. A region's
// answer brings, as dynamic nodes, the transcripts whose coding exons it overlaps: each a coding ranges node on the
// region (minus strand for TP53) and, for TP53's, a nested source of its domains on the transcript's residues. The
// isoforms and the neighbor are made up; the domains are the protein's, cut to each isoform's length.

/** Where the gene's letters start in the region. */
export const GeneOffset = 2000;
export const RegionLength = 6000;
/** Positions a region's answer covers: each brings the transcripts it overlaps. */
export const RegionSize = 500;

/** The region's source: its address in the view's spec, its system's name and its committed file. */
export const RegionSource = '/data/region';
/** The region's sequence in the spec: the transcripts are on it. */
export const RegionSequence = '/data/tp53-region';
export const RegionSystem = 'chr17-tp53';
export const RegionFile = 'fixtures/p53/gene-region.seq';
/** The format the view names: the application gives it a runtime (a delayed `seq` one in the stories). */
export const RegionFormat = 'seq-async';

interface Transcript {
  readonly strand: '+' | '-';
  /** Coding exons on the region, genomic order. */
  readonly exons: readonly (readonly [number, number])[];
  /** Whether its protein's domains come as a nested source. */
  readonly domains: boolean;
}

const at = (e: readonly (readonly [number, number])[]) => e.map(([a, b]) => [a + GeneOffset, b + GeneOffset] as const);
const [x1, x2, x3] = GeneExons.start.map((s, i) => [s, GeneExons.end[i]] as const);

/** TP53's isoforms (T1 is the protein of the other views) and WRAP53, a neighbor on the plus strand. */
export const Transcripts: Readonly<Record<string, Transcript>> = {
  T1: { strand: '-', exons: at([x1, x2, x3]), domains: true },
  T2: { strand: '-', exons: at([x1, [x2[0], 780]]), domains: true },
  T3: { strand: '-', exons: at([x3]), domains: true },
  W1: {
    strand: '+',
    exons: [
      [3700, 3900],
      [4300, 4700],
    ],
    domains: false,
  },
};

const codingLength = (t: Transcript) => t.exons.reduce((n, [a, b]) => n + b - a, 0);

/** The region's letters: seeded, with the gene's at its offset. */
export function regionLetters(): string {
  const r = seeded(29);
  let s = '';
  for (let i = 0; i < RegionLength; i++) s += 'ACGT'[Math.floor(r() * 4)];
  return s.slice(0, GeneOffset) + geneLetters() + s.slice(GeneOffset + GeneLength);
}

/** A transcript's domains as a nested source: one region over its residues, a ranges node on them up front. */
function domainsData(id: string, format: string): SeqData {
  const length = Math.floor(codingLength(Transcripts[id]) / 3);
  const nested = nodeAddress(RegionSource, `${id}-domains`);
  const node: Node = {
    kind: 'ranges',
    params: { target: itemAddress(nodeAddress(RegionSource, id), id), source: nested },
  };
  const kept = Domains.id.map((_, i) => i).filter((i) => Domains.start[i] < length);
  return {
    coverage: {
      systems: [{ name: id, length }],
      shapes: [{ system: id, shape: 'ranges', fields: [{ name: 'label', type: 'category' }] }],
      levels: [1],
      nodes: [{ key: 'domains', node }],
    },
    pairs: [
      {
        request: { system: id, start: 0, end: length, level: 1, shapes: ['ranges'] },
        answer: {
          ok: true,
          tables: {
            ranges: {
              form: 'list',
              start: Int32Array.from(kept, (i) => Domains.start[i]),
              end: Int32Array.from(kept, (i) => Math.min(Domains.end[i], length)),
              id: categoryColumn(kept.map((i) => Domains.id[i])),
              fields: [{ name: 'label', column: categoryColumn(kept.map((i) => Domains.label[i])) }],
            },
          },
        },
      },
    ],
  };
}

/**
 * A gene's exons as a ranges node on `target` (a row an exon, each with the gene's id and strand), through a codon map
 * when `map` says so: what a by-region answer brings for each gene or transcript it overlaps.
 */
export function exonsNode(
  target: string,
  id: string,
  strand: '+' | '-',
  exons: readonly (readonly [number, number])[],
  map?: 'codon',
): Node {
  const rows = {
    start: exons.map((e) => e[0]),
    end: exons.map((e) => e[1]),
    id: exons.map(() => id),
    strand: exons.map(() => strand),
  };
  return { kind: 'ranges', params: { target, ...(map ? { map } : {}), rows } };
}

/** The nodes a transcript brings: its coding ranges node and, for TP53's, the nested source of its domains. */
function transcriptNodes(id: string, format: string): DynamicNode[] {
  const t = Transcripts[id];
  const node = exonsNode(RegionSequence, id, t.strand, t.exons, 'codon');
  const out: DynamicNode[] = [{ key: id, node }];
  if (t.domains)
    out.push({
      key: `${id}-domains`,
      node: { kind: 'source', params: { format, value: seqToValue(domainsData(id, format)) } },
    });
  return out;
}

/** The region as `seq` data: its letters and transcripts, region by region. */
export function regionData(format = RegionFormat): SeqData {
  const letters = regionLetters();
  const pairs: SeqData['pairs'][number][] = [];
  for (let start = 0; start < RegionLength; start += RegionSize) {
    const end = Math.min(RegionLength, start + RegionSize);
    const ids = Object.keys(Transcripts).filter((id) => Transcripts[id].exons.some(([a, b]) => a < end && b > start));
    const nodes = ids.flatMap((id) => transcriptNodes(id, format));
    pairs.push({
      request: { system: RegionSystem, start, end, level: 1, shapes: ['letters'] },
      answer: {
        ok: true,
        tables: { letters: { form: 'letters', start, letters: letters.slice(start, end) } },
        ...(nodes.length ? { nodes } : {}),
      },
    });
  }
  return {
    coverage: {
      systems: [{ name: RegionSystem, length: RegionLength, alphabet: 'dna' }],
      shapes: [{ system: RegionSystem, shape: 'letters' }],
      levels: [1],
    },
    pairs,
  };
}

/** The committed file's text. */
export function regionSeqText(): string {
  return `${encodeSeq(regionData())}\n`;
}
