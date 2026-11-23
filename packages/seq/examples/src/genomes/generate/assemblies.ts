/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqData } from '@molstar/seqx-core/data';
import { between, Concern, stream } from './streams.js';

// Step 1 of the Genomes generator (examples.md, Genomes, The generator: outline): toy A, prefixed (`chr1` to `chr3`,
// the cyclic `chrM` and an HLA contig whose name carries `*` and `:`, EC-200), its letters a 4 kb tile at a time from
// a stream per `(seed, sequence, tile)` (a GC share waving along the sequence, a few CA repeats, a run of N inside
// chr1), its lengths as `chrom-sizes`, a small FASTA of the two short sequences for the built FASTA reader, and an alias
// table in the `seq` format (each sequence's bare number and an `NC_`-style accession; `MT` left out, so a file naming
// chrM `MT` is refused, EC-189). Toy B (unprefixed) is the alignment's (`alignment.ts`), A′ the liftover's.

export interface ToySequence {
  readonly name: string;
  readonly length: number;
  readonly cyclic?: boolean;
}

export interface ToyAssembly {
  readonly name: string;
  readonly sequences: readonly ToySequence[];
}

/** The HLA contig's name (GRCh38's alt contigs' spelling). */
export const HlaName = 'HLA-A*01:01:01:01';

/** The toy lengths (`genomes.toy.lengths`): chr1, chr2, chr3. */
export const DefaultToyLengths: readonly [number, number, number] = [300_000, 200_000, 120_000];

/** Toy A for the given chromosome lengths: chr1 to chr3, then `chrM` (16,569, cyclic) and the HLA contig (3,503). */
export function toyA(lengths: readonly [number, number, number] = DefaultToyLengths): ToyAssembly {
  return {
    name: 'Toy A',
    sequences: [
      { name: 'chr1', length: lengths[0] },
      { name: 'chr2', length: lengths[1] },
      { name: 'chr3', length: lengths[2] },
      { name: 'chrM', length: 16_569, cyclic: true },
      { name: HlaName, length: 3_503 },
    ],
  };
}

/** The run of N inside chr1 (a gap of the assembly: 20,000 letters, nothing aligned to it). */
export const NGap: readonly [number, number] = [140_000, 160_000];

/** A letters tile. */
const Tile = 4096;

/** The letters of toy A's sequences, by name: tile by tile from their own streams. */
export function toyALetters(a: ToyAssembly, seed: number): Map<string, string> {
  const out = new Map<string, string>();
  a.sequences.forEach((s, i) => {
    const chars = new Array<string>(s.length);
    for (let t = 0; t * Tile < s.length; t++) {
      const rng = stream(seed, Concern.letters, i * 4096 + t);
      for (let p = t * Tile; p < Math.min(s.length, (t + 1) * Tile); p++) {
        const gc = 0.41 + 0.14 * Math.sin((2 * Math.PI * p) / 47_000 + i);
        const r = rng.float();
        const k = rng.float() < 0.5 ? 0 : 1;
        chars[p] = r < gc ? 'GC'[k] : 'AT'[k];
      }
    }
    // a few CA repeats: 60 to 200 letters each, at places from the sequence's own stream
    const rep = stream(seed, Concern.letters, 1_000_000 + i);
    for (let k = 0; k < Math.max(1, Math.floor(s.length / 60_000)); k++) {
      const at = between(rep, 0, s.length - 200);
      const n = between(rep, 30, 100);
      for (let j = 0; j < 2 * n; j++) chars[at + j] = j % 2 ? 'A' : 'C';
    }
    if (s.name === 'chr1') for (let p = NGap[0]; p < Math.min(NGap[1], s.length); p++) chars[p] = 'N';
    out.set(s.name, chars.join(''));
  });
  return out;
}

/** A `chrom-sizes` file (UCSC's): a name and a length a line, in the assembly's order. */
export function chromSizes(a: ToyAssembly): string {
  return a.sequences.map((s) => `${s.name}\t${s.length}\n`).join('');
}

/** A FASTA of `names`' letters, 60 a line. */
export function fasta(letters: ReadonlyMap<string, string>, names: readonly string[]): string {
  let out = '';
  for (const n of names) {
    const l = letters.get(n)!;
    out += `>${n}\n`;
    for (let i = 0; i < l.length; i += 60) out += `${l.slice(i, i + 60)}\n`;
  }
  return out;
}

/** Toy A's aliases: chr1 to chr3 their bare numbers and accessions, chrM its accession alone (no `MT`). */
export function toyAAliases(a: ToyAssembly): Map<string, readonly string[]> {
  const out = new Map<string, readonly string[]>();
  a.sequences.forEach((s, i) => {
    const accession = `NC_9000${String(i + 1).padStart(2, '0')}.1`;
    const m = /^chr(\d+)$/.exec(s.name);
    if (m) out.set(s.name, [m[1], accession]);
    else if (s.name === 'chrM') out.set(s.name, [accession]);
  });
  return out;
}

/**
 * The alias table as a `seq` source (a chromAlias-like table: its coverage's systems, their lengths and aliases, no
 * shapes): what an assembly's `alias_source` reads.
 */
export function aliasTable(a: ToyAssembly): SeqData {
  const aliases = toyAAliases(a);
  return {
    coverage: {
      systems: a.sequences.map((s) => ({
        name: s.name,
        length: s.length,
        ...(aliases.has(s.name) ? { aliases: [...aliases.get(s.name)!] } : {}),
      })),
      shapes: [],
      levels: [1],
    },
    pairs: [],
  };
}
