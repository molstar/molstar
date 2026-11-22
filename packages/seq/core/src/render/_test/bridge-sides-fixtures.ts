/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge between two genomes over a sided pairs source (generated): assemblies `P` and `G` of two chromosomes each,
// named alike (`chr1`, `chr2` on both: a sided source lists each name once a side), alignments of blocks between them
// served by region from either side at levels 1, 1 kb and 16 kb, and two genome sections with a bridge drawing them.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { PairsList, SeqData, SeqPair } from '@molstar/seqx-core/data';
import { binPairs, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';

export const ChrLength = 1_000_000;
export const Chromosomes = ['chr1', 'chr2'] as const;
const Levels = [1, 1024, 16_384];

/** A block between the genomes: its chromosome and start on each, its length, its strand on `G`. */
export interface SideBlock {
  readonly a: string;
  readonly as: number;
  readonly b: string;
  readonly bs: number;
  readonly length: number;
  readonly strand: 1 | -1;
}

/** `count` alignments of 2 to 4 blocks (200 to 600 long) between random chromosomes, seeded. */
export function sideBlocks(count: number, seed: number): SideBlock[] {
  const r = new PCG(seed);
  const out: SideBlock[] = [];
  for (let i = 0; i < count; i++) {
    const a = Chromosomes[r.int() % 2];
    const b = Chromosomes[r.int() % 2];
    let as = r.int() % (ChrLength - 5000);
    let bs = r.int() % (ChrLength - 5000);
    const strand = r.int() % 3 ? 1 : -1;
    const m = 2 + (r.int() % 3);
    for (let j = 0; j < m; j++) {
      const length = 200 + (r.int() % 400);
      out.push({ a, as, b, bs, length, strand });
      as += length + 1 + (r.int() % 99);
      bs += length + 1 + (r.int() % 99);
    }
  }
  return out;
}

function listOf(blocks: readonly SideBlock[]): PairsList {
  const n = blocks.length;
  return {
    form: 'list',
    aSystem: categoryColumn(blocks.map((x) => x.a)),
    bSystem: categoryColumn(blocks.map((x) => x.b)),
    aStart: Int32Array.from(blocks, (x) => x.as),
    aEnd: Int32Array.from(blocks, (x) => x.as + x.length),
    bStart: Int32Array.from(blocks, (x) => x.bs),
    bEnd: Int32Array.from(blocks, (x) => x.bs + x.length),
    aStrand: new Int8Array(n).fill(1),
    bStrand: Int8Array.from(blocks, (x) => x.strand),
  };
}

/** The blocks served by region on each side (a name listed once a side), the whole list a request, bins above 1. */
export function sidedSource(blocks: readonly SideBlock[]): SeqData {
  const list = listOf(blocks);
  const pairs: SeqPair[] = [];
  for (const side of ['a', 'b'] as const)
    for (const system of Chromosomes)
      for (const level of Levels)
        pairs.push({
          request: { system, start: 0, end: ChrLength, level, shapes: ['pairs'], side },
          answer: { ok: true, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
        });
  return {
    coverage: {
      systems: (['a', 'b'] as const).flatMap((side) => Chromosomes.map((name) => ({ name, length: ChrLength, side }))),
      shapes: Chromosomes.map((system) => ({ system, shape: 'pairs' as const })),
      levels: Levels,
    },
    pairs,
  };
}

export const TopSection = '/views/v/top';
export const BottomSection = '/views/v/bottom';

/** `P` above, `G` below (genome sections), a bridge between them drawing the blocks read by region. */
export function sidesBridgeSpec(blocks: readonly SideBlock[]) {
  const b = createSeqSpec({ title: 'two genomes, sided' });
  const seqs = Chromosomes.map((name) => ({ name, length: ChrLength }));
  const P = b.data.assembly('P', assembly('P', seqs, { order: 'source' } as never));
  const G = b.data.assembly('G', assembly('G', seqs, { order: 'source' } as never));
  const src = b.data.source('src', { format: 'seq', value: seqToValue(sidedSource(blocks)) as never });
  const aln = b.data.pairs('aln', { a: P, b: G, source: src } as never);
  b.views.view('v', (v) => {
    genome(v, 'top', P);
    v.bridge('br', { a: TopSection, b: BottomSection } as never, (br) =>
      br.layer('aln', { type: 'pairs', data: aln } as never),
    );
    genome(v, 'bottom', G);
  });
  return b.build();
}
