/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import { binPairs, type GridSparse, type PairsList, type SeqData, type SeqPair } from '@molstar/seqx-core/data';
import { NGap, type ToyAssembly } from './assemblies.js';
import { between, Concern, stream } from './streams.js';

// Step 2 of the Genomes generator: toy B made from toy A by one alignment truth, and the
// truth written down. B names its sequences unprefixed (`1`, `2`, `3`, `MT`: the alias rule exercised, EC-189); B's
// `1` holds A's chr1 with an 80 kb inversion (EC-206) and its own letters where A has its run of N (a gap: nothing
// aligned), B's `2` is A's chr2 with A's chr3 tail after it (a 60 kb translocation), B's `3` the rest of A's chr3.
// Each piece is aligned in gapless blocks of 3 to 9 kb with a small indel between two (an insertion in B or a
// deletion), its letters A's with 1 % substituted; the truth is the blocks, a run of equal length on each genome with
// its strand, the blocks of one piece one alignment (`row`). Written as sided pairs in the `seq` format (side `a` A's
// sequences, `b` B's, EC-208), and as the dot plot's grid: the aligned positions counted in cells of two levels.

/** A piece of B: A's `[a0, a1)` on `a` (absent: B's own letters, `length` of them), placed on B read `strand`. */
interface Piece {
  readonly a?: string;
  readonly a0: number;
  readonly a1: number;
  readonly strand: 1 | -1;
}

/** A's chr1 stretch B holds reversed (an inversion, 80 kb). */
export const Inversion: readonly [number, number] = [40_000, 120_000];
/** The length of A's chr3 tail carried to the end of B's `2` (a translocation). */
export const TranslocationLength = 60_000;

/** B's sequences as pieces of A, in order (lengths from A's). */
function pieces(a: ToyAssembly): readonly (readonly [name: string, pieces: readonly Piece[]])[] {
  const len = (n: string) => a.sequences.find((s) => s.name === n)!.length;
  const [l1, l2, l3] = [len('chr1'), len('chr2'), len('chr3')];
  const t = l3 - TranslocationLength;
  return [
    [
      '1',
      [
        { a: 'chr1', a0: 0, a1: Inversion[0], strand: 1 },
        { a: 'chr1', a0: Inversion[0], a1: Inversion[1], strand: -1 },
        { a: 'chr1', a0: Inversion[1], a1: NGap[0], strand: 1 },
        { a0: NGap[0], a1: NGap[1], strand: 1 },
        { a: 'chr1', a0: NGap[1], a1: l1, strand: 1 },
      ],
    ],
    [
      '2',
      [
        { a: 'chr2', a0: 0, a1: l2, strand: 1 },
        { a: 'chr3', a0: t, a1: l3, strand: 1 },
      ],
    ],
    ['3', [{ a: 'chr3', a0: 0, a1: t, strand: 1 }]],
    ['MT', [{ a: 'chrM', a0: 0, a1: len('chrM'), strand: 1 }]],
  ];
}

/** An aligned block: A's `[aStart, aEnd)` and B's `[bStart, bEnd)`, equal lengths, B read `strand`; `row` its piece. */
export interface Block {
  readonly a: string;
  readonly aStart: number;
  readonly aEnd: number;
  readonly b: string;
  readonly bStart: number;
  readonly bEnd: number;
  readonly strand: 1 | -1;
  readonly row: number;
}

export interface AlignmentTruth {
  /** Toy B: its sequences' names and lengths. */
  readonly b: ToyAssembly;
  /** The blocks, in B's order. */
  readonly blocks: readonly Block[];
  /** Toy B's letters, by name. */
  readonly letters: ReadonlyMap<string, string>;
}

const Complement: Readonly<Record<string, string>> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };

/** Toy B, its letters and the blocks aligning it to A. */
export function alignToyB(a: ToyAssembly, aLetters: ReadonlyMap<string, string>, seed: number): AlignmentTruth {
  const blocks: Block[] = [];
  const sequences: { name: string; length: number; cyclic?: boolean }[] = [];
  const letters = new Map<string, string>();
  let row = 0;
  for (const [name, list] of pieces(a)) {
    const out: string[] = [];
    for (const p of list) {
      const rng = stream(seed, Concern.alignment, row);
      const random = () => 'ACGT'[Math.floor(rng.float() * 4)];
      const at = out.length;
      if (!p.a) {
        for (let i = p.a0; i < p.a1; i++) out.push(random());
        row++;
        continue;
      }
      const src = aLetters.get(p.a)!;
      // the piece laid forward: blocks of A with B's offsets, an indel between two
      const laid: { a0: number; a1: number; b0: number }[] = [];
      const piece: string[] = [];
      for (let x = p.a0; x < p.a1; ) {
        const n = Math.min(p.a1 - x, between(rng, 3_000, 9_000));
        laid.push({ a0: x, a1: x + n, b0: piece.length });
        for (let i = x; i < x + n; i++) piece.push(rng.float() < 0.01 ? random() : src[i]);
        x += n;
        if (x >= p.a1) break;
        const d = between(rng, 1, 26);
        if (rng.float() < 0.5) x = Math.min(p.a1, x + d);
        else for (let i = 0; i < d; i++) piece.push(random());
      }
      const total = piece.length;
      if (p.strand > 0) for (const c of piece) out.push(c);
      else for (let i = total - 1; i >= 0; i--) out.push(Complement[piece[i]] ?? piece[i]);
      for (const l of laid) {
        const n = l.a1 - l.a0;
        const b0 = p.strand > 0 ? at + l.b0 : at + total - l.b0 - n;
        blocks.push({ a: p.a, aStart: l.a0, aEnd: l.a1, b: name, bStart: b0, bEnd: b0 + n, strand: p.strand, row });
      }
      row++;
    }
    sequences.push({ name, length: out.length, ...(name === 'MT' ? { cyclic: true } : {}) });
    letters.set(name, out.join(''));
  }
  return { b: { name: 'Toy B', sequences }, blocks, letters };
}

/** B's position of A's `p` on `a` through the blocks (undefined where `p` is not aligned). */
export function toB(blocks: readonly Block[], a: string, p: number): { b: string; p: number } | undefined {
  for (const k of blocks) {
    if (k.a !== a || p < k.aStart || p >= k.aEnd) continue;
    return { b: k.b, p: k.strand > 0 ? k.bStart + (p - k.aStart) : k.bEnd - 1 - (p - k.aStart) };
  }
  return undefined;
}

/** The blocks as a pairs list (`a` on A, `b` on B), those of `keep`. */
function pairsList(blocks: readonly Block[]): PairsList {
  return {
    form: 'list',
    aSystem: categoryColumn(blocks.map((k) => k.a)),
    bSystem: categoryColumn(blocks.map((k) => k.b)),
    aStart: Int32Array.from(blocks, (k) => k.aStart),
    aEnd: Int32Array.from(blocks, (k) => k.aEnd),
    bStart: Int32Array.from(blocks, (k) => k.bStart),
    bEnd: Int32Array.from(blocks, (k) => k.bEnd),
    aStrand: new Int8Array(blocks.length).fill(1),
    bStrand: Int8Array.from(blocks, (k) => k.strand),
    row: Int32Array.from(blocks, (k) => k.row),
  };
}

/**
 * Two assemblies' runs as a sided `seq` source: `a`'s sequences on side `a`, `b`'s on side `b`, a request on a side
 * answered with the runs whose end on that side is on its sequence (found from either end, EC-208). Level 1 alone by
 * default (the store's summaries give the pairs their levels), else bin pairs at each coarser level of `levels`.
 */
export function sidedPairs(
  a: ToyAssembly,
  b: ToyAssembly,
  runs: readonly Pick<Block, 'a' | 'aStart' | 'aEnd' | 'b' | 'bStart' | 'bEnd' | 'strand' | 'row'>[],
  levels: readonly number[] = [1],
): SeqData {
  const sides = [
    ['a', a.sequences],
    ['b', b.sequences],
  ] as const;
  const pairs: SeqPair[] = [];
  for (const [side, seqs] of sides)
    for (const s of seqs) {
      const list = pairsList(runs.filter((r) => (side === 'a' ? r.a : r.b) === s.name) as Block[]);
      for (const level of levels)
        pairs.push({
          request: { system: s.name, start: 0, end: s.length, level, shapes: ['pairs'], side },
          answer: { ok: true, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
        });
    }
  return {
    coverage: {
      systems: sides.flatMap(([side, seqs]) => seqs.map((s) => ({ name: s.name, length: s.length, side }))),
      shapes: a.sequences.map((s) => ({ system: s.name, shape: 'pairs' as const })),
      levels: [...levels],
    },
    pairs,
  };
}

/** The dot plot's levels: cells of 1,024 and 4,096 positions a side. */
export const DotLevels: readonly number[] = [1024, 4096];

/**
 * The dot plot's grid (x A's sequences, y B's): at each level each cell the positions aligned inside it, sparse (the
 * rest 0); a pair of sequences with none aligned answers an empty grid.
 */
export function dotGrid(a: ToyAssembly, b: ToyAssembly, blocks: readonly Block[]): SeqData {
  const pairs: SeqPair[] = [];
  for (const level of DotLevels)
    for (const sx of a.sequences)
      for (const sy of b.sequences) {
        const cells = new Map<number, number>();
        const ny = Math.ceil(sy.length / level);
        for (const k of blocks) {
          if (k.a !== sx.name || k.b !== sy.name) continue;
          for (let p = k.aStart; p < k.aEnd; p++) {
            const q = k.strand > 0 ? k.bStart + (p - k.aStart) : k.bEnd - 1 - (p - k.aStart);
            const key = Math.floor(p / level) * ny + Math.floor(q / level);
            cells.set(key, (cells.get(key) ?? 0) + 1);
          }
        }
        const keys = [...cells.keys()].sort((x, y) => x - y);
        const grid: GridSparse = {
          form: 'sparse',
          x: 0,
          y: 0,
          xBin: Int32Array.from(keys, (k) => Math.floor(k / ny)),
          yBin: Int32Array.from(keys, (k) => k % ny),
          value: Float32Array.from(keys, (k) => cells.get(k)!),
          fill: 0,
        };
        pairs.push({
          request: {
            system: sx.name,
            start: 0,
            end: sx.length,
            level,
            shapes: ['grid'],
            y: { system: sy.name, start: 0, end: sy.length },
          },
          answer: { ok: true, tables: { grid } },
        });
      }
  return {
    coverage: {
      systems: [...a.sequences, ...b.sequences].map((s) => ({ name: s.name, length: s.length })),
      shapes: a.sequences.map((s) => ({
        system: s.name,
        shape: 'grid' as const,
        y: b.sequences.map((o) => o.name),
        aggregate: 'sum' as const,
      })),
      levels: [...DotLevels],
    },
    pairs,
  };
}
