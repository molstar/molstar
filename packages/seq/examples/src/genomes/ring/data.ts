/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';

// The genome ring's data (the showcase of `@molstar/seqx-ext-rondo`): GRCh38's chromosomes by their lengths (facts from
// the GRC's assembly report, typed here: no file, no fetch), a seeded generator of links between them (60 % within a
// chromosome with a distance decay, 40 % between, a dozen hub pairs and hot spots, a `type` and a `q`), written held
// (inline rows) and as BEDPE (one picture), and synthetic cytobands for tests. Randomness from core's `PCG` only:
// byte-identical across runs.

/** GRCh38's primary chromosomes (1 to 22, X) and Y, the 24th, by their lengths (GRCh38.p14's assembly report). */
export const Grch38Chromosomes: readonly { readonly name: string; readonly length: number }[] = [
  ...[
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
    50818468,
  ].map((length, i) => ({ name: `chr${i + 1}`, length })),
  { name: 'chrX', length: 156040895 },
  { name: 'chrY', length: 57227415 },
];

/** The ring's 23 (1 to 22, X): Y joins as the 24th by `include`. */
export const RingInclude: readonly string[] = Grch38Chromosomes.slice(0, 23).map((c) => c.name);

/** A link: its ends (a chromosome's index and a start), its kind and its q-value. */
export interface RingLink {
  readonly a: number;
  readonly as: number;
  readonly b: number;
  readonly bs: number;
  readonly type: 'contact' | 'other';
  readonly q: number;
}

/** A link's end's length in bases. */
export const LinkEnd = 5000;

/** The hub pairs (indices into the chromosomes): strong clusters between them. */
const Hubs: readonly (readonly [number, number])[] = [
  [0, 16],
  [3, 6],
  [1, 10],
  [8, 21],
  [4, 4],
  [5, 12],
  [2, 19],
  [7, 14],
  [9, 22],
  [11, 17],
  [13, 18],
  [15, 20],
];

/**
 * `n` links over the first `count` chromosomes (`seed`): 60 % within a chromosome, the partner at a distance that
 * decays (most within a few Mb), 40 % between; half of the between ones on the hub pairs, each hub at two hot spots
 * (a cluster a few Mb wide), so stems and strong ribbons exist. `type` `contact` (70 %) or `other`; `q` in (0, 0.05].
 */
export function ringLinks(n = 12_000, seed = 38, count = 23): RingLink[] {
  const chroms = Grch38Chromosomes.slice(0, count);
  const rng = new PCG(seed);
  const f = () => rng.float();
  const spot = (c: number, k: number) => Math.floor(((k * 0.37 + 0.2) % 0.9) * (chroms[c].length - 5e6));
  const out: RingLink[] = [];
  for (let i = 0; i < n; i++) {
    let a: number;
    let b: number;
    let as: number;
    let bs: number;
    const u = f();
    if (u < 0.6) {
      a = b = Math.floor(f() * count);
      const len = chroms[a].length - LinkEnd;
      as = Math.floor(f() * len);
      // the distance decays: a power of a uniform, a few Mb at most most of the time
      const d = Math.floor(f() ** 4 * 2e7) * (f() < 0.5 ? -1 : 1);
      bs = Math.max(0, Math.min(len, as + d));
    } else if (u < 0.8) {
      const h = Hubs[Math.floor(f() * Hubs.length)];
      [a, b] = h[0] < count && h[1] < count ? h : [0, 1];
      const k = Math.floor(f() * 2);
      as = Math.min(chroms[a].length - LinkEnd, spot(a, k) + Math.floor(f() * 3e6));
      bs = Math.min(chroms[b].length - LinkEnd, spot(b, k + 1) + Math.floor(f() * 3e6));
    } else {
      a = Math.floor(f() * count);
      b = Math.floor(f() * count);
      as = Math.floor(f() * (chroms[a].length - LinkEnd));
      bs = Math.floor(f() * (chroms[b].length - LinkEnd));
    }
    const type = f() < 0.7 ? 'contact' : 'other';
    const q = Math.max(1e-6, Math.round(f() * 0.05 * 1e6) / 1e6);
    out.push({ a, as, b, bs, type, q });
  }
  return out;
}

/** The links as BEDPE (`chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 type q`): read whole. */
export function ringBedpe(links: readonly RingLink[]): string {
  const name = (i: number) => Grch38Chromosomes[i].name;
  const lines = links.map((l, i) =>
    [name(l.a), l.as, l.as + LinkEnd, name(l.b), l.bs, l.bs + LinkEnd, `link${i + 1}`, 0, '+', '+', l.type, l.q].join(
      '\t',
    ),
  );
  return `${lines.join('\n')}\n`;
}

/** The links as a pairs node's inline rows over the assembly (`a_name`, `b_name`): held. */
export function ringRows(links: readonly RingLink[]) {
  const name = (i: number) => Grch38Chromosomes[i].name;
  return {
    a_name: links.map((l) => name(l.a)),
    a_start: Int32Array.from(links, (l) => l.as),
    a_end: Int32Array.from(links, (l) => l.as + LinkEnd),
    b_name: links.map((l) => name(l.b)),
    b_start: Int32Array.from(links, (l) => l.bs),
    b_end: Int32Array.from(links, (l) => l.bs + LinkEnd),
    type: links.map((l) => l.type),
    q: Float64Array.from(links, (l) => l.q),
  };
}

/** The committed BEDPE of 1,000 links (seed 38): the reader's tests read it as a file. */
export const RingLinks1000File = 'fixtures/genomes/ring-links-1000.bedpe';

/**
 * Synthetic cytobands (`chrom start end name stain`, BED 4+1): alternating `gneg` and `gpos25`…`gpos100` of lengths
 * drawn from the seed, an `acen` pair at 40 % of each length, a `gvar` on the last third of a chromosome named `Y`
 * (`chrY`) when present: every test and story without the recording is deterministic.
 */
export function syntheticBands(sequences: readonly { readonly name: string; readonly length: number }[], seed = 5) {
  const rng = new PCG(seed);
  const pos = ['gpos25', 'gpos50', 'gpos75', 'gpos100'];
  const lines: string[] = [];
  for (const s of sequences) {
    const cen = Math.floor(s.length * 0.4);
    const half = Math.max(1, Math.floor(s.length * 0.01));
    const gvar = /(^|chr)Y$/.test(s.name) ? Math.floor((s.length * 2) / 3) : s.length;
    let at = 0;
    let k = 0;
    const band = (end: number, stain: string) => {
      if (end <= at) return;
      lines.push([s.name, at, end, `${s.name.replace(/^chr/, '')}b${++k}`, stain].join('\t'));
      at = end;
    };
    let neg = true;
    while (at < s.length) {
      if (at < cen - half && at + 1 > 0) {
        const len = Math.floor((0.01 + rng.float() * 0.03) * s.length);
        const end = Math.min(cen - half, at + len);
        band(end, neg ? 'gneg' : pos[rng.int() % 4]);
        neg = !neg;
        if (end === cen - half) {
          band(cen, 'acen');
          band(cen + half, 'acen');
        }
        continue;
      }
      if (at >= gvar) {
        band(s.length, 'gvar');
        break;
      }
      const len = Math.floor((0.01 + rng.float() * 0.03) * s.length);
      band(Math.min(gvar, at + len), neg ? 'gneg' : pos[rng.int() % 4]);
      neg = !neg;
    }
  }
  return `${lines.join('\n')}\n`;
}
