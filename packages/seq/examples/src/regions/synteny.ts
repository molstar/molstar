/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { Coverage, PairsList, SeqData, SeqPair, SystemInfo } from '@molstar/seqx-core/data';

// A generated synteny alignment between two toy genomes (Features/Region reads): aligned
// blocks as sided pairs in the `seq` format, a pair the block's extent on genome A and on genome B (the `a` side names
// A's chromosomes, the `b` side B's), found from either end. B is A with an inversion on its second chromosome and 6 Mb
// inserted before it, so most blocks are colinear and some are reversed. Deterministic on core's `PCG`.

export interface SyntenyChromosome {
  readonly name: string;
  readonly length: number;
}

/** Genome A's chromosomes and genome B's (B's `chrB2` is A's `chrA2` after 6 Mb of its own). */
export const SyntenyA: readonly SyntenyChromosome[] = [
  { name: 'chrA1', length: 40_000_000 },
  { name: 'chrA2', length: 30_000_000 },
  { name: 'chrA3', length: 20_000_000 },
];
export const SyntenyB: readonly SyntenyChromosome[] = [
  { name: 'chrB1', length: 40_000_000 },
  { name: 'chrB2', length: 36_000_000 },
  { name: 'chrB3', length: 20_000_000 },
];

/** The assemblies' names as a view writes them. */
export const SyntenyAssemblyA = 'Toy A';
export const SyntenyAssemblyB = 'Toy B';

/** The blocks: about 6,000, past the pairs layer's gate at the genome and a few dozen in a few Mb. */
export const SyntenyBlocks = 6_000;

interface Block {
  readonly a: string;
  readonly aStart: number;
  readonly aEnd: number;
  readonly b: string;
  readonly bStart: number;
  readonly bEnd: number;
  readonly strand: 1 | -1;
}

/** An inverted stretch of A's `chrA2` (positions), where B reads it reversed. */
export const SyntenyInversion: readonly [number, number] = [8_000_000, 14_000_000];

/** The blocks walking A's chromosomes in order, B's end by the map above. */
function blocks(n: number, seed: number): Block[] {
  const rng = new PCG(seed);
  const total = SyntenyA.reduce((a, c) => a + c.length, 0);
  const pitch = Math.floor(total / n);
  const out: Block[] = [];
  const [i0, i1] = SyntenyInversion;
  for (const c of SyntenyA) {
    for (let at = 20_000; at + pitch < c.length; at += pitch) {
      const aStart = at + Math.floor(rng.float() * pitch * 0.2);
      const aEnd = Math.min(c.length, aStart + Math.floor(pitch * (0.35 + 0.4 * rng.float())));
      let b: string;
      let bStart: number;
      let bEnd: number;
      let strand: 1 | -1 = 1;
      if (c.name === 'chrA1') {
        b = 'chrB1';
        bStart = aStart;
        bEnd = aEnd;
      } else if (c.name === 'chrA2') {
        b = 'chrB2';
        bStart = aStart + 6_000_000;
        bEnd = aEnd + 6_000_000;
        if (aStart >= i0 && aEnd <= i1) {
          // read reversed: the stretch mirrored on itself
          strand = -1;
          bStart = i0 + i1 - aEnd + 6_000_000;
          bEnd = i0 + i1 - aStart + 6_000_000;
        }
      } else {
        b = 'chrB3';
        bStart = aStart;
        bEnd = aEnd;
      }
      out.push({ a: c.name, aStart, aEnd, b, bStart, bEnd: Math.min(bEnd, lengthOf(b)), strand });
    }
  }
  return out;
}

const lengthOf = (name: string) => SyntenyB.find((c) => c.name === name)!.length;

/** The blocks as the table both sides serve, rows ordered by A's end. */
export function syntenyPairs(n = SyntenyBlocks, seed = 13): PairsList {
  const bl = blocks(n, seed);
  return {
    form: 'list',
    aSystem: categoryColumn(bl.map((b) => b.a)),
    bSystem: categoryColumn(bl.map((b) => b.b)),
    aStart: Int32Array.from(bl, (b) => b.aStart),
    aEnd: Int32Array.from(bl, (b) => b.aEnd),
    bStart: Int32Array.from(bl, (b) => b.bStart),
    bEnd: Int32Array.from(bl, (b) => b.bEnd),
    aStrand: Int8Array.from(bl, () => 1),
    bStrand: Int8Array.from(bl, (b) => b.strand),
    row: Int32Array.from(bl, (_, i) => i),
  };
}

/** The rows of `t` with an end on `system`, on `side` only when it names one (`a`: its first end), as a table of their own. */
function endingOn(t: PairsList, system: string, side?: 'a' | 'b'): PairsList {
  const rows: number[] = [];
  for (let i = 0; i < t.aStart.length; i++) {
    const onA = t.aSystem!.get(i) === system;
    const onB = t.bSystem!.get(i) === system;
    if (side === 'a' ? onA : side === 'b' ? onB : onA || onB) rows.push(i);
  }
  const take = <A extends Int32Array | Int8Array>(c: A) =>
    (c.constructor as unknown as { from(a: number[]): A }).from(rows.map((i) => c[i]));
  const cat = (c: NonNullable<PairsList['aSystem']>) => categoryColumn(rows.map((i) => c.get(i)));
  return {
    form: 'list',
    aSystem: cat(t.aSystem!),
    bSystem: cat(t.bSystem!),
    aStart: take(t.aStart),
    aEnd: take(t.aEnd),
    bStart: take(t.bStart),
    bEnd: take(t.bEnd),
    aStrand: take(t.aStrand!),
    bStrand: take(t.bStrand!),
    row: take(t.row!),
  };
}

export interface SyntenyOptions {
  /**
   * Sided (the default): A's chromosomes on side `a`, B's on side `b`, a request naming its side and answering the pairs
   * with that end (a map between genomes, whose names may collide). Not sided: one list of systems, a request on one
   * answering the pairs with an end on it either way, what a pairs layer over a section's sequences reads.
   */
  readonly sided?: boolean;
}

/**
 * The alignment as pairs in the `seq` format, level 1 alone (the store's summaries give the pairs their levels): a
 * pair is found from either end.
 */
export function syntenySource(n = SyntenyBlocks, seed = 13, options: SyntenyOptions = {}): SeqData {
  const sided = options.sided !== false;
  const all = syntenyPairs(n, seed);
  const systems: SystemInfo[] = [
    ...SyntenyA.map((c) => ({ name: c.name, length: c.length, ...(sided && { side: 'a' as const }) })),
    ...SyntenyB.map((c) => ({ name: c.name, length: c.length, ...(sided && { side: 'b' as const }) })),
  ];
  const coverage: Coverage = {
    systems,
    shapes: [...SyntenyA, ...SyntenyB].map((c) => ({ system: c.name, shape: 'pairs' as const })),
    levels: [1],
  };
  const pairs: SeqPair[] = systems.map((s) => ({
    request: { system: s.name, start: 0, end: s.length!, level: 1, shapes: ['pairs'], ...(s.side && { side: s.side }) },
    answer: { ok: true, tables: { pairs: endingOn(all, s.name, s.side) } },
  }));
  return { coverage, pairs };
}
