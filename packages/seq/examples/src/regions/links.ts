/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { Coverage, PairsList, SeqData, SeqPair } from '@molstar/seqx-core/data';
import { binPairs } from '@molstar/seqx-core/data';
import { HicChromosomes } from './hic.js';

// 12,000 generated links (Features/Region reads): the same links as a table held
// whole (one system), as a BEDPE file's text, and as a `seq` source over a toy assembly's chromosomes, so one picture
// is drawn from either. A link joins two of a few dozen hotspots (a chromosome's arms, a translocation's pair), its
// ends jittered, so bins at a coarse level hold counts of very different strengths. Deterministic on core's `PCG`.

/** The links' number. */
export const LinkCount = 12_000;

/** One system the links are on: its name and length. */
export interface LinkSystem {
  readonly name: string;
  readonly length: number;
}

/** The single chromosome the held links are on. */
export const LinksChromosome: LinkSystem = { name: 'chr', length: 2_000_000 };

/** The toy assembly's chromosomes the ring's links join (the Hi-C's). */
export const LinksChromosomes: readonly LinkSystem[] = HicChromosomes;

interface Links {
  readonly aSystem: string[];
  readonly bSystem: string[];
  readonly aStart: Int32Array;
  readonly bStart: Int32Array;
  readonly aStrand: Int8Array;
  readonly bStrand: Int8Array;
}

/** An end's size: ten positions. */
const EndLength = 10;

/**
 * `n` links between hotspots of `systems` (20 hotspots on each system, a link joining two of them, the stronger ones
 * chosen more often, a quarter within a hotspot's chromosome only), each end jittered by 2 kb at most.
 */
function generate(systems: readonly LinkSystem[], n: number, seed: number): Links {
  const rng = new PCG(seed);
  const spots = systems.map((s) =>
    Array.from({ length: 20 }, (_, i) => Math.floor((s.length * (i + 0.5)) / 20 - EndLength)),
  );
  // a few strong pairs of hotspots draw most links
  const strong = Array.from({ length: 12 }, () => [
    Math.floor(rng.float() * systems.length),
    Math.floor(rng.float() * 20),
    Math.floor(rng.float() * systems.length),
    Math.floor(rng.float() * 20),
  ]);
  const out: Links = {
    aSystem: [],
    bSystem: [],
    aStart: new Int32Array(n),
    bStart: new Int32Array(n),
    aStrand: new Int8Array(n),
    bStrand: new Int8Array(n),
  };
  for (let i = 0; i < n; i++) {
    let sa: number;
    let ha: number;
    let sb: number;
    let hb: number;
    if (rng.float() < 0.6) {
      [sa, ha, sb, hb] = strong[Math.floor(rng.float() * strong.length)];
    } else {
      sa = Math.floor(rng.float() * systems.length);
      ha = Math.floor(rng.float() * 20);
      sb = rng.float() < 0.25 ? sa : Math.floor(rng.float() * systems.length);
      hb = Math.floor(rng.float() * 20);
    }
    const jitter = () => Math.floor(rng.float() * 2000);
    out.aSystem.push(systems[sa].name);
    out.bSystem.push(systems[sb].name);
    out.aStart[i] = Math.min(systems[sa].length - EndLength, spots[sa][ha] + jitter());
    out.bStart[i] = Math.min(systems[sb].length - EndLength, spots[sb][hb] + jitter());
    out.aStrand[i] = 1;
    out.bStrand[i] = rng.float() < 0.85 ? 1 : -1;
  }
  return out;
}

function listOf(l: Links, sided: boolean): PairsList {
  const n = l.aStart.length;
  return {
    form: 'list',
    ...(sided && { aSystem: categoryColumn(l.aSystem), bSystem: categoryColumn(l.bSystem) }),
    aStart: l.aStart,
    aEnd: l.aStart.map((s) => s + EndLength),
    bStart: l.bStart,
    bEnd: l.bStart.map((s) => s + EndLength),
    aStrand: l.aStrand,
    bStrand: l.bStrand,
    id: categoryColumn(Array.from({ length: n }, (_, i) => `l${i}`)),
  };
}

/** The links on one chromosome (`LinksChromosome`), as the table a held `pairs` node's rows hold. */
export function linksOnChromosome(n = LinkCount, seed = 9): PairsList {
  return listOf(generate([LinksChromosome], n, seed), false);
}

/** The same links as a BEDPE file's text (tab separated: both ends, a name, a score, both strands). */
export function linksBedpe(n = LinkCount, seed = 9, system: LinkSystem = LinksChromosome): string {
  const l = generate([system], n, seed);
  const lines: string[] = [];
  const strand = (s: number) => (s < 0 ? '-' : '+');
  for (let i = 0; i < n; i++)
    lines.push(
      [
        system.name,
        l.aStart[i],
        l.aStart[i] + EndLength,
        system.name,
        l.bStart[i],
        l.bStart[i] + EndLength,
        `l${i}`,
        1,
        strand(l.aStrand[i]),
        strand(l.bStrand[i]),
      ].join('\t'),
    );
  return `${lines.join('\n')}\n`;
}

/** The levels a generated links source serves: the items, and bins of 1 kb and 16 kb. */
export const LinkLevels: readonly number[] = [1, 1024, 16_384];

/**
 * Links between the toy assembly's chromosomes as a `seq` source: each chromosome answers the links with an end on it,
 * at the items' level and at two coarser ones (bin pairs by count: the pairs level, served as a source with levels).
 */
export function linksSource(systems: readonly LinkSystem[] = LinksChromosomes, n = LinkCount, seed = 9): SeqData {
  const all = generate(systems, n, seed);
  const full = listOf(all, true);
  const coverage: Coverage = {
    systems: systems.map((s) => ({ name: s.name, length: s.length })),
    shapes: systems.map((s) => ({ system: s.name, shape: 'pairs' as const })),
    levels: [...LinkLevels],
  };
  const pairs: SeqPair[] = [];
  for (const s of systems) {
    const rows: number[] = [];
    for (let i = 0; i < n; i++) if (all.aSystem[i] === s.name || all.bSystem[i] === s.name) rows.push(i);
    const take = <A extends Int32Array | Int8Array>(c: A) =>
      (c.constructor as unknown as { from(a: number[]): A }).from(rows.map((i) => c[i]));
    const cat = (c: NonNullable<PairsList['aSystem']>) => categoryColumn(rows.map((i) => c.get(i)));
    const own: PairsList = {
      form: 'list',
      aSystem: cat(full.aSystem!),
      bSystem: cat(full.bSystem!),
      aStart: take(full.aStart),
      aEnd: take(full.aEnd),
      bStart: take(full.bStart),
      bEnd: take(full.bEnd),
      aStrand: take(full.aStrand!),
      bStrand: take(full.bStrand!),
      id: cat(full.id!),
    };
    for (const level of LinkLevels)
      pairs.push({
        request: { system: s.name, start: 0, end: s.length, level, shapes: ['pairs'] },
        answer: { ok: true, tables: { pairs: level === 1 ? own : binPairs(own, level) } },
      });
  }
  return { coverage, pairs };
}
