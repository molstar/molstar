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
import type { FormatRuntime, PairsList, SeqData, SeqPair, StoreClock } from '@molstar/seqx-core/data';
import { binPairs, delayedSource, openSeq, seqToValue } from '@molstar/seqx-core/data';
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

/** The name side `b` lists a chromosome by (`aliased`: `chr1` as `1`). */
const bName = (name: string, aliased: boolean) => (aliased ? name.replace('chr', '') : name);

function listOf(blocks: readonly SideBlock[], aliased = false): PairsList {
  const n = blocks.length;
  return {
    form: 'list',
    aSystem: categoryColumn(blocks.map((x) => x.a)),
    bSystem: categoryColumn(blocks.map((x) => bName(x.b, aliased))),
    aStart: Int32Array.from(blocks, (x) => x.as),
    aEnd: Int32Array.from(blocks, (x) => x.as + x.length),
    bStart: Int32Array.from(blocks, (x) => x.bs),
    bEnd: Int32Array.from(blocks, (x) => x.bs + x.length),
    aStrand: new Int8Array(n).fill(1),
    bStrand: Int8Array.from(blocks, (x) => x.strand),
  };
}

/**
 * The blocks served by region on each side (a name listed once a side), the whole list a request, bins above 1; `served`
 * the sides answered (another listed but not served: its reads fail); `aliased`, side `b` names `chr1` as `1`.
 */
export function sidedSource(
  blocks: readonly SideBlock[],
  served: readonly ('a' | 'b')[] = ['a', 'b'],
  aliased = false,
): SeqData {
  const list = listOf(blocks, aliased);
  const pairs: SeqPair[] = [];
  for (const side of served)
    for (const chr of Chromosomes)
      for (const level of Levels)
        pairs.push({
          request: {
            system: side === 'b' ? bName(chr, aliased) : chr,
            start: 0,
            end: ChrLength,
            level,
            shapes: ['pairs'],
            side,
          },
          answer: { ok: true, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
        });
  return {
    coverage: {
      systems: (['a', 'b'] as const).flatMap((side) =>
        Chromosomes.map((name) => ({ name: side === 'b' ? bName(name, aliased) : name, length: ChrLength, side })),
      ),
      shapes: [...new Set([...Chromosomes, ...Chromosomes.map((c) => bName(c, aliased))])].map((system) => ({
        system,
        shape: 'pairs' as const,
      })),
      levels: Levels,
    },
    pairs,
  };
}

export const TopSection = '/views/v/top';
export const BottomSection = '/views/v/bottom';

/** The format of a source whose side `b` answers on a clock of the test's (`slowSide`). */
export const SlowFormat = 'seq-slow-b';

/** `sidedSource` with side `b`'s reads waiting on `clock` (a timer of 1 ms each), side `a`'s at once. */
export function slowSide(blocks: readonly SideBlock[], clock: StoreClock): FormatRuntime {
  return {
    name: SlowFormat,
    formats: [SlowFormat],
    open: async () => delayedSource(openSeq(sidedSource(blocks)), { clock, delayMs: (r) => (r.side === 'b' ? 1 : 0) }),
  };
}

/** A clock whose timers wait until `flush` is called, and fire at once from then on. */
export function manualClock(): StoreClock & { flush(): void } {
  let due: (() => void)[] | undefined = [];
  return {
    now: () => 0,
    setTimeout: (fn) => (due ? due.push(fn) : queueMicrotask(fn)),
    clearTimeout: () => {},
    flush() {
      const now = due ?? [];
      due = undefined;
      for (const f of now) f();
    },
  };
}

/**
 * `P` above, `G` below (genome sections), a bridge between them drawing the blocks read by region (`served`: the sides
 * the source answers; `slow`: the source of `slowSide`'s format instead).
 */
export function sidesBridgeSpec(
  blocks: readonly SideBlock[],
  options: { served?: readonly ('a' | 'b')[]; slow?: boolean; self?: boolean; aliased?: boolean } = {},
) {
  const b = createSeqSpec({ title: 'two genomes, sided' });
  const seqs = Chromosomes.map((name) => ({ name, length: ChrLength }));
  const P = b.data.assembly('P', assembly('P', seqs, { order: 'source' } as never));
  const G = b.data.assembly('G', assembly('G', seqs, { order: 'source' } as never));
  const src = options.slow
    ? b.data.source('src', { format: SlowFormat, uri: 'slow' } as never)
    : b.data.source('src', {
        format: 'seq',
        value: seqToValue(sidedSource(blocks, options.served, options.aliased)) as never,
      });
  // `self`: one genome against itself (`P` on both ends and both sides)
  const aln = b.data.pairs('aln', { a: P, b: options.self ? P : G, source: src } as never);
  b.views.view('v', (v) => {
    genome(v, 'top', P);
    v.bridge('br', { a: TopSection, b: BottomSection } as never, (br) =>
      br.layer('aln', { type: 'pairs', data: aln } as never),
    );
    genome(v, 'bottom', options.self ? P : G);
  });
  return b.build();
}
