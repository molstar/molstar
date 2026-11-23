/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';

// The genome ring's data (the showcase of `@molstar/seqx-ext-rondo`): GRCh38's chromosomes by their lengths (facts from
// the GRC's assembly report, typed here: no file, no fetch), a seeded generator of links between them (60 % within a
// chromosome with a distance decay, 40 % between, a few named translocations, made up, and the contacts about them; a
// `type`, strands, a `q`, a fusion's genes), written held
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

/** A kind of link: a 3D contact, or a translocation (a structural variant's junction). */
export type RingLinkType = 'contact' | 'translocation';

/** A link: its ends (a chromosome's index, a start and a strand), its kind, its q-value, its genes (a fusion's). */
export interface RingLink {
  readonly a: number;
  readonly as: number;
  readonly b: number;
  readonly bs: number;
  readonly type: RingLinkType;
  readonly q: number;
  /** Each end's strand: `+`, `-`, or `.` (a contact's: none). */
  readonly sa: '+' | '-' | '.';
  readonly sb: '+' | '-' | '.';
  /** A fusion's two genes (made up); absent elsewhere. */
  readonly genes?: readonly [string, string];
}

/** A link's end's length in bases. */
export const LinkEnd = 5000;

/**
 * A translocation of the showcase, **made up**: its name, its two junctions (a chromosome's index and a position), the
 * strands its reads join, and for a fusion its two genes. The coordinates are invented, not any real breakpoint's, and
 * the genes invented names (`X…`): the first is in the style of BCR-ABL1 (a chr9–chr22 junction), named so it is never
 * read as the real one.
 */
export interface RingTranslocation {
  readonly name: string;
  readonly a: number;
  readonly at: number;
  readonly b: number;
  readonly bt: number;
  readonly sa: '+' | '-';
  readonly sb: '+' | '-';
  readonly genes?: readonly [string, string];
}

/** The showcase's named translocations (invented): where the generator's hub pairs were. */
export const RingTranslocations: readonly RingTranslocation[] = [
  {
    name: 't(9;22) XABL::XBCR (made up)',
    a: 8,
    at: 101_300_000,
    b: 21,
    bt: 31_800_000,
    sa: '+',
    sb: '-',
    genes: ['XABL', 'XBCR'],
  },
  {
    name: 't(1;17) XETS::XRUN (made up)',
    a: 0,
    at: 61_200_000,
    b: 16,
    bt: 44_700_000,
    sa: '-',
    sb: '+',
    genes: ['XETS', 'XRUN'],
  },
  { name: 't(4;7) (made up, no fusion)', a: 3, at: 92_500_000, b: 6, bt: 120_400_000, sa: '+', sb: '+' },
  {
    name: 't(12;18) XPAX::XFOX (made up)',
    a: 11,
    at: 38_900_000,
    b: 17,
    bt: 52_100_000,
    sa: '-',
    sb: '-',
    genes: ['XPAX', 'XFOX'],
  },
  { name: 't(2;11) (made up, no fusion)', a: 1, at: 170_600_000, b: 10, bt: 77_300_000, sa: '+', sb: '-' },
];

/**
 * `n` links over the first `count` chromosomes (`seed`): 60 % within a chromosome, the partner at a distance that
 * decays (most within a few Mb), 40 % between; half of the between ones about the named translocations' junctions
 * (`RingTranslocations`): about 6 in a hundred the translocation's own reads within 20 kb of it (`type`
 * `translocation`, its strands, a fusion's genes), the rest contacts enriched about it over a few Mb, so stems and
 * strong ribbons exist. Every other link a `contact`, unstranded; `q` in (0, 0.05].
 */
export function ringLinks(n = 12_000, seed = 38, count = 23): RingLink[] {
  const chroms = Grch38Chromosomes.slice(0, count);
  const rng = new PCG(seed);
  const f = () => rng.float();
  const clamp = (c: number, x: number) => Math.max(0, Math.min(chroms[c].length - LinkEnd, Math.floor(x)));
  const out: RingLink[] = [];
  for (let i = 0; i < n; i++) {
    let a: number;
    let b: number;
    let as: number;
    let bs: number;
    let t: RingTranslocation | undefined;
    const u = f();
    if (u < 0.6) {
      a = b = Math.floor(f() * count);
      const len = chroms[a].length - LinkEnd;
      as = Math.floor(f() * len);
      // the distance decays: a power of a uniform, a few Mb at most most of the time
      const d = Math.floor(f() ** 4 * 2e7) * (f() < 0.5 ? -1 : 1);
      bs = Math.max(0, Math.min(len, as + d));
    } else if (u < 0.8) {
      const x = RingTranslocations[Math.floor(f() * RingTranslocations.length)];
      const on = x.a < count && x.b < count;
      [a, b] = on ? [x.a, x.b] : [0, 1];
      // its own reads within 20 kb of the junction; contacts about it, decaying over a few Mb
      const own = on && f() < 0.06;
      if (own) t = x;
      const spread = own ? 2e4 : 3e6 * f() ** 2;
      as = clamp(a, (on ? x.at : chroms[a].length / 2) + (f() - 0.5) * 2 * spread);
      bs = clamp(b, (on ? x.bt : chroms[b].length / 2) + (f() - 0.5) * 2 * spread);
    } else {
      a = Math.floor(f() * count);
      b = Math.floor(f() * count);
      as = Math.floor(f() * (chroms[a].length - LinkEnd));
      bs = Math.floor(f() * (chroms[b].length - LinkEnd));
    }
    const q = Math.max(1e-6, Math.round(f() * 0.05 * 1e6) / 1e6);
    if (t)
      out.push({ a, as, b, bs, type: 'translocation', q, sa: t.sa, sb: t.sb, ...(t.genes ? { genes: t.genes } : {}) });
    else out.push({ a, as, b, bs, type: 'contact', q, sa: '.', sb: '.' });
  }
  return out;
}

/**
 * The links as BEDPE (`chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 type q gene1 gene2`, `.` where
 * none): read whole.
 */
export function ringBedpe(links: readonly RingLink[]): string {
  const name = (i: number) => Grch38Chromosomes[i].name;
  const lines = links.map((l, i) =>
    [
      name(l.a),
      l.as,
      l.as + LinkEnd,
      name(l.b),
      l.bs,
      l.bs + LinkEnd,
      `link${i + 1}`,
      0,
      l.sa,
      l.sb,
      l.type,
      l.q,
      l.genes?.[0] ?? '.',
      l.genes?.[1] ?? '.',
    ].join('\t'),
  );
  return `${lines.join('\n')}\n`;
}

/** The links as a pairs node's inline rows over the assembly (`a_name`, `b_name`; strands; genes where any): held. */
export function ringRows(links: readonly RingLink[]) {
  const name = (i: number) => Grch38Chromosomes[i].name;
  return {
    a_name: links.map((l) => name(l.a)),
    a_start: Int32Array.from(links, (l) => l.as),
    a_end: Int32Array.from(links, (l) => l.as + LinkEnd),
    b_name: links.map((l) => name(l.b)),
    b_start: Int32Array.from(links, (l) => l.bs),
    b_end: Int32Array.from(links, (l) => l.bs + LinkEnd),
    a_strand: links.map((l) => l.sa),
    b_strand: links.map((l) => l.sb),
    type: links.map((l) => l.type),
    q: Float64Array.from(links, (l) => l.q),
    gene_a: links.map((l) => l.genes?.[0] ?? null),
    gene_b: links.map((l) => l.genes?.[1] ?? null),
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
