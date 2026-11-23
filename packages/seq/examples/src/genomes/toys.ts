/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import {
  delayedSource,
  FastaRuntime,
  type FormatRuntime,
  type OpenedSource,
  openSeq,
  type SeqData,
  SourceOpenError,
  type StoreClock,
} from '@molstar/seqx-core/data';
import { assembly, type createSeqSpec } from '@molstar/seqx-core/lang';
import { BedFormats, BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { GffFormats, GffRuntime } from '@molstar/seqx-core/runtimes/gff';
import { VcfFormats, VcfRuntime } from '@molstar/seqx-core/runtimes/vcf';
import { alignToyB } from './generate/alignment.js';
import { toyA, toyALetters } from './generate/assemblies.js';
import { GenomesFiles } from './generate/index.js';
import { GenomesSeed } from './generate/streams.js';

// The Genomes dataset's views share this: the committed files' sources (each view's builder
// takes the files as an argument, these by default), the assemblies read from their `chrom-sizes` files (no
// manifest: the lengths are the file's, EC-188), and the registry: the opt-in genome runtimes and the FASTA one, each
// answering after a delay on a clock (the synthetic twin the stories and tests share), and the in-memory
// letters adapter of both toys (`toy-letters`: letters by region from the seed, never stored).

export { GenomesFiles };

export type Builder = ReturnType<typeof createSeqSpec>;

/** The files a view reads: the committed ones by default (a builder takes another set of the same shape). */
export type Files = Readonly<Record<keyof typeof GenomesFiles, string>>;

/** The format the letters adapter reads, and the `uri` naming each toy. */
export const LettersFormat = 'toy-letters';
export const LettersA = 'generated:toy-a';
export const LettersB = 'generated:toy-b';

/** A source of `format` at `uri` (with its tabix `index`, named: never `<uri>.tbi`, EC-186). */
export const fileSource = (b: Builder, id: string, format: string, uri: string, index?: string, extra = {}) =>
  b.data.source(id, { format, uri, ...(index !== undefined ? { index } : {}), ...extra } as never);

/** Toy A, its sequences and lengths from its chrom.sizes. */
export const assemblyA = (b: Builder, f: Files, id = 'toy-a') =>
  b.data.assembly(id, { name: 'Toy A', source: fileSource(b, `${id}-sizes`, 'chrom-sizes', f.sizesA) } as never);
/** Toy B (unprefixed names), from its chrom.sizes. */
export const assemblyB = (b: Builder, f: Files, id = 'toy-b') =>
  b.data.assembly(id, { name: 'Toy B', source: fileSource(b, `${id}-sizes`, 'chrom-sizes', f.sizesB) } as never);
/** Toy A′ (A revised), from its chrom.sizes. */
export const assemblyA2 = (b: Builder, f: Files, id = 'toy-a2') =>
  b.data.assembly(id, { name: 'Toy A′', source: fileSource(b, `${id}-sizes`, 'chrom-sizes', f.sizesA2) } as never);

/**
 * Toy B's sequences as the generator makes them by default (a test checks them against it): inline where a view
 * needs the lengths before any file opens (a grid section following a genome's).
 */
export const ToyBSequences: readonly { name: string; length: number; cyclic?: boolean }[] = [
  { name: '1', length: 299_757 },
  { name: '2', length: 259_952 },
  { name: '3', length: 60_032 },
  { name: 'MT', length: 16_594, cyclic: true },
];
/** Toy A inline (its default lengths). */
export const inlineA = (b: Builder, id = 'toy-a', match?: 'exact') =>
  b.data.assembly(id, assembly('Toy A', toyA().sequences, match ? { match } : {}));
/** Toy B inline. */
export const inlineB = (b: Builder, id = 'toy-b', match?: 'exact') =>
  b.data.assembly(id, assembly('Toy B', ToyBSequences, match ? { match } : {}));

/** Both toys' letters, generated once per seed (toy B's through the alignment truth). */
const lettersCache = new Map<number, { a: SeqData; b: SeqData }>();
function lettersOf(seed: number) {
  let got = lettersCache.get(seed);
  if (!got) {
    const a = toyA();
    const la = toyALetters(a, seed);
    const aligned = alignToyB(a, la, seed);
    const data = (
      seqs: readonly { name: string; length: number; cyclic?: boolean }[],
      l: ReadonlyMap<string, string>,
    ) =>
      ({
        coverage: {
          systems: seqs.map((s) => ({
            name: s.name,
            length: s.length,
            alphabet: 'dna',
            ...(s.cyclic && { cyclic: true }),
          })),
          shapes: seqs.map((s) => ({ system: s.name, shape: 'letters' as const })),
          levels: [1],
        },
        pairs: seqs.map((s) => ({
          request: { system: s.name, start: 0, end: s.length, level: 1, shapes: ['letters' as const] },
          answer: {
            ok: true as const,
            tables: { letters: { form: 'letters' as const, start: 0, letters: l.get(s.name)! } },
          },
        })),
      }) as SeqData;
    got = { a: data(a.sequences, la), b: data(aligned.b.sequences, aligned.letters) };
    lettersCache.set(seed, got);
  }
  return got;
}

/** The letters adapter: toy A's or B's sequences, lengths and letters by region, from the seed (`genomes.seed`). */
export function lettersRuntime(seed = GenomesSeed): FormatRuntime {
  return {
    name: 'toy-letters',
    formats: [LettersFormat],
    async open(source): Promise<OpenedSource> {
      const which = source.uri === LettersA ? 'a' : source.uri === LettersB ? 'b' : undefined;
      if (!which) throw new SourceOpenError([]);
      return openSeq(lettersOf(seed)[which]);
    },
  };
}

export interface GenomesRegistryOptions {
  /** Each answer's delay, ms (150; 0 in tests: no timer). */
  readonly delayMs?: number;
  /** The clock the answers wait on: real timers by default. */
  readonly clock?: StoreClock;
}

/** What the dataset's views need besides the default registry: each runtime delayed on the clock. */
export function registry(options: GenomesRegistryOptions = {}): SeqRegistryEntry[] {
  const timing = { delayMs: options.delayMs ?? 150, ...(options.clock ? { clock: options.clock } : {}) };
  // each source delayed, keeping the inner one's name index (a GFF3 read whole finds its genes by name, EC-228)
  const delayed = (runtime: FormatRuntime, name: string, formats: readonly string[]): SeqRegistryEntry => ({
    kind: 'format-runtime',
    runtime: {
      name,
      formats,
      async open(source, ctx) {
        const inner = await runtime.open(source, ctx);
        const out = delayedSource(inner, timing);
        return inner.findNames ? { ...out, findNames: inner.findNames.bind(inner) } : out;
      },
    },
  });
  return [
    delayed(BedRuntime, 'bed-async', BedFormats),
    delayed(VcfRuntime, 'vcf-async', VcfFormats),
    delayed(GffRuntime, 'gff-async', GffFormats),
    delayed(FastaRuntime, 'fasta-async', ['fasta']),
    delayed(lettersRuntime(), 'toy-letters-async', [LettersFormat]),
  ];
}
