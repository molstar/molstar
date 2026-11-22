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

/**
 * The genome ring's chromosomes: the Hi-C's six at 100,000 times their length (38 Mb down to 10 Mb), so the ring at the
 * fit reads the links' bins (16 kb), never the links themselves.
 */
export const RingChromosomes: readonly LinkSystem[] = HicChromosomes.map((c) => ({ ...c, length: c.length * 100_000 }));

interface Links {
  readonly aSystem: string[];
  readonly bSystem: string[];
  readonly aStart: Int32Array;
  readonly bStart: Int32Array;
  readonly aStrand: Int8Array;
  readonly bStrand: Int8Array;
}

/** The least distance between a local link's ends: four bins of 16 kb. */
const LocalGap = 65_536;

/** The coarsest level's bin. */
const LastBin = 16_384;

/** An end's size: ten positions. */
const EndLength = 10;

/**
 * `n` links between hotspots of `systems` (20 hotspots on each system, a link joining two of them, the stronger ones
 * chosen more often, a quarter within a hotspot's chromosome only), each end jittered by 2 kb at most (a twentieth of a
 * short chromosome), never before its start.
 */
function generate(systems: readonly LinkSystem[], n: number, seed: number): Links {
  const rng = new PCG(seed);
  const spots = systems.map((s) =>
    Array.from({ length: 20 }, (_, i) => Math.max(0, Math.floor((s.length * (i + 0.5)) / 20 - EndLength))),
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
    // within a twentieth of a short system (a toy assembly's chromosomes of a few hundred positions)
    const jitter = (s: LinkSystem) => Math.floor(rng.float() * Math.min(2000, s.length / 20));
    out.aSystem.push(systems[sa].name);
    out.bSystem.push(systems[sb].name);
    out.aStart[i] = Math.min(systems[sa].length - EndLength, spots[sa][ha] + jitter(systems[sa]));
    out.bStart[i] = Math.min(systems[sb].length - EndLength, spots[sb][hb] + jitter(systems[sb]));
    out.aStrand[i] = 1;
    out.bStrand[i] = rng.float() < 0.85 ? 1 : -1;
  }
  return out;
}

/**
 * `n` links over a genome ring as structure, not a disc (synteny blocks and structural variants): a handful of strong
 * blocks between two chromosomes (a stretch of one joined to a stretch of another, the ends running together or
 * inverted), local links within a chromosome (short, a tandem or a loop), and a thin scatter. Each strong block is a
 * region pair a few hundredths of the chromosomes' lengths wide, so its bins at a coarse level are a band of chords.
 */
function generateStructured(systems: readonly LinkSystem[], n: number, seed: number): Links {
  const rng = new PCG(seed);
  const out: Links = {
    aSystem: [],
    bSystem: [],
    aStart: new Int32Array(n),
    bStart: new Int32Array(n),
    aStrand: new Int8Array(n),
    bStrand: new Int8Array(n),
  };
  const pick = (k: number) => Math.min(k - 1, Math.floor(rng.float() * k));
  /** A stretch of a system: its start and width as shares of its length. */
  const stretch = (s: LinkSystem, width: number) => {
    const w = Math.floor(s.length * width);
    return { at: Math.floor(rng.float() * (s.length - w - EndLength)), width: w };
  };
  // the strong blocks: pairs of chromosomes (never alike), a stretch each, their ends in order or inverted
  const blocks: {
    a: number;
    b: number;
    sa: { at: number; width: number };
    sb: { at: number; width: number };
    inverted: boolean;
  }[] = [];
  for (let k = 0; k < 7; k++) {
    const a = pick(systems.length);
    let b = pick(systems.length - 1);
    if (b >= a) b++;
    const width = 0.015 + 0.02 * rng.float();
    blocks.push({ a, b, sa: stretch(systems[a], width), sb: stretch(systems[b], width), inverted: rng.float() < 0.3 });
  }
  // the local links: tandem duplications and loops, each a stretch of one chromosome read again a fixed span on (so
  // their chords run together as a band, never a fan)
  const locals = Array.from({ length: 10 }, () => {
    const s = pick(systems.length);
    const span = stretch(systems[s], 0.02 + 0.02 * rng.float());
    return { s, span, offset: Math.max(LocalGap, Math.floor(systems[s].length * (0.05 + 0.1 * rng.float()))) };
  });
  for (let i = 0; i < n; i++) {
    const r = rng.float();
    let sa: number;
    let sb: number;
    let pa: number;
    let pb: number;
    if (r < 0.7) {
      // a strong block: a position along its stretch on one chromosome joined to the same share on the other
      const bl = blocks[pick(blocks.length)];
      const t = rng.float();
      sa = bl.a;
      sb = bl.b;
      pa = bl.sa.at + Math.floor(t * bl.sa.width);
      pb = bl.sb.at + Math.floor((bl.inverted ? 1 - t : t) * bl.sb.width);
    } else if (r < 0.985) {
      // a local link: a position in the stretch joined to the one a span on (the ends are two bins at least apart)
      const l = locals[pick(locals.length)];
      sa = sb = l.s;
      pa = l.span.at + Math.floor(rng.float() * l.span.width);
      pb = pa + l.offset;
    } else {
      // a thin scatter
      sa = pick(systems.length);
      sb = pick(systems.length);
      pa = Math.floor(rng.float() * (systems[sa].length - EndLength));
      pb = Math.floor(rng.float() * (systems[sb].length - EndLength));
    }
    out.aSystem.push(systems[sa].name);
    out.bSystem.push(systems[sb].name);
    // never in a chromosome's last bin (it ends inside it: a chord's end must be placed)
    out.aStart[i] = Math.max(0, Math.min(systems[sa].length - LastBin - EndLength, pa));
    out.bStart[i] = Math.max(0, Math.min(systems[sb].length - LastBin - EndLength, pb));
    out.aStrand[i] = 1;
    out.bStrand[i] = 1;
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
 * at the items' level and at two coarser ones (bin pairs by count: the pairs level, served as a source with levels);
 * `structured`: synteny blocks and local links (`generateStructured`), what a genome ring reads as structure.
 */
export function linksSource(
  systems: readonly LinkSystem[] = LinksChromosomes,
  n = LinkCount,
  seed = 9,
  structured = false,
): SeqData {
  const all = structured ? generateStructured(systems, n, seed) : generate(systems, n, seed);
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
