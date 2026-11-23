/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Drawing against the data (_spec/README.md, Stems, Ribbons, Stubs): every pair read with an end in view is in exactly
// one ribbon, stub or fan (or within its stem); a neck is `k ×` its stem's count and the sum of its ribbons, stubs and
// fan (a hub of 300 partners too); a stem's count the pairs with an end in it that leave it; a ribbon's ends inside
// its stems’ necks; stubs only for ends with no place, one a stem a partner; a reversed block, pairs across a tile edge.

import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { bundledTextMetrics } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { rondo } from '../index.js';
import { Kind, type RondoLay } from '../lay.js';
import {
  genomeChromosomes,
  lay,
  type Link,
  linksBedpe,
  ringSpec,
  Size,
  settle,
  toyChromosomes,
  toyLinks,
} from './helpers.js';

/** Each row's two stems (by the lay's stem keys), its strand, its weight; rows the filter leaves out skipped. */
function rows(l: RondoLay) {
  const read = l.read!;
  const d = read.data!;
  const index = new Map<string, number>();
  for (let s = 0; s < l.ns; s++) index.set(`${l.stemSystem[s]}:${l.stemBin[s]}`, s);
  const sysA = (k: number) => (d.aUnit ? d.aUnit[k] : 0);
  const sysB = (k: number) => (d.bUnit ? d.bUnit[k] : l.systems.length > 1 ? 1 : 0);
  // a row whose end is in no stem is a stem missed: the test fails on it, never drops the row from a total
  const stemOf = (k: number, sys: number, start: number) => {
    const key = `${sys}:${Math.floor(start / l.stemLevel)}`;
    const s = index.get(key);
    if (s === undefined) throw new Error(`row ${k}: no stem at ${key}`);
    return s;
  };
  const out: { a: number; b: number; strand: number; w: number }[] = [];
  for (let k = 0; k < d.n; k++) {
    if (!read.shown[k]) continue;
    const a = stemOf(k, sysA(k), d.aStart[k]);
    const b = stemOf(k, sysB(k), d.bStart[k]);
    out.push({ a, b, strand: d.aStrand[k] * d.bStrand[k] < 0 ? -1 : 1, w: read.weight ? read.weight[k] : 1 });
  }
  return out;
}

/**
 * Each pair's one place, by rule (within its stem; the ribbon of its stem pair and strand, drawn or in the fan; the stub
 * of its placed stem toward the partner), summed by place and compared with what the layout holds there: a double
 * count and a miss cannot cancel. Every ribbon, stub and within count the layout holds is some pairs' exactly.
 */
function accounted(l: RondoLay) {
  const ribbons = new Map<string, number>();
  const stubs = new Map<string, number>();
  const within = new Float64Array(l.ns);
  let total = 0;
  for (const r of rows(l)) {
    const pa = l.stemPlaced[r.a];
    const pb = l.stemPlaced[r.b];
    if (!pa && !pb) continue;
    total += r.w;
    if (r.a === r.b) within[r.a] += r.w;
    else if (pa && pb) {
      const key = `${Math.min(r.a, r.b)}|${Math.max(r.a, r.b)}|${r.strand}`;
      ribbons.set(key, (ribbons.get(key) ?? 0) + r.w);
    } else {
      const key = `${pa ? r.a : r.b}|${l.stemSystem[pa ? r.b : r.a]}`;
      stubs.set(key, (stubs.get(key) ?? 0) + r.w);
    }
  }
  const wrong: string[] = [];
  const held = new Set<string>();
  for (let r = 0; r < l.nr; r++) {
    const key = `${Math.min(l.ribA[r], l.ribB[r])}|${Math.max(l.ribA[r], l.ribB[r])}|${l.ribStrand[r]}`;
    if (held.has(key)) wrong.push(`ribbon ${key} twice`);
    held.add(key);
    if (ribbons.get(key) !== l.ribCount[r]) wrong.push(`ribbon ${key}: ${l.ribCount[r]}, pairs ${ribbons.get(key)}`);
  }
  for (const key of ribbons.keys()) if (!held.has(key)) wrong.push(`pairs of ${key} in no ribbon`);
  for (let b = 0; b < l.nb; b++) {
    const key = `${l.stubStem[b]}|${l.stubPartner[b]}`;
    if (held.has(`stub ${key}`)) wrong.push(`stub ${key} twice`);
    held.add(`stub ${key}`);
    if (stubs.get(key) !== l.stubCount[b]) wrong.push(`stub ${key}: ${l.stubCount[b]}, pairs ${stubs.get(key)}`);
  }
  for (const key of stubs.keys()) if (!held.has(`stub ${key}`)) wrong.push(`pairs of ${key} in no stub`);
  for (let s = 0; s < l.ns; s++)
    if (l.stemPlaced[s] && within[s] !== l.stemWithin[s]) wrong.push(`within ${s}: ${l.stemWithin[s]}, ${within[s]}`);
  return { total, wrong };
}

/** Each placed stem's members' counts (drawn ribbons and stubs at it) and its fan, against its count. */
function necksAdd(l: RondoLay) {
  const sum = new Float64Array(l.ns);
  for (let r = 0; r < l.nr; r++)
    if (l.ribDrawn[r]) {
      sum[l.ribA[r]] += l.ribCount[r];
      sum[l.ribB[r]] += l.ribCount[r];
    }
  for (let b = 0; b < l.nb; b++) if (l.stubDrawn[b]) sum[l.stubStem[b]] += l.stubCount[b];
  const bad: number[] = [];
  for (let s = 0; s < l.ns; s++)
    if (l.stemPlaced[s] && Math.abs(sum[s] + l.stemFan[s] - l.stemCount[s]) > 1e-6) bad.push(s);
  return bad;
}

/** chr1 at 1 Mb linked once each to 300 places spread over the other chromosomes, and 40 times to chr3 at 2 Mb. */
function hubSpec() {
  const chroms = toyChromosomes();
  const links: Link[] = [];
  for (let i = 0; i < 300; i++) {
    const b = 1 + (i % 5);
    links.push({
      a: 0,
      b,
      as: 1_000_000 + (i % 50),
      bs: Math.floor((i * 997_331) % (chroms[b].length - 1000)),
      type: 'contact',
    });
  }
  for (let i = 0; i < 40; i++) links.push({ a: 0, b: 2, as: 1_000_100, bs: 2_000_000 + i, type: 'contact' });
  return ringSpec({ chroms, links, held: true });
}

describe('drawing against the data', () => {
  it('every pair with an end in view is in exactly one ribbon, stub or fan, or within its stem (EC-171)', async () => {
    for (const o of [{}, { held: true }, { include: ['chr1'] }, { include: ['chr2', 'chr5'] }]) {
      const s = await settle(ringSpec(o));
      const l = lay(await s.frame({ size: Size }))!;
      const a = accounted(l);
      expect(a.total).toBeGreaterThan(0);
      expect(a.wrong).toEqual([]);
    }
  });

  it('a neck is the sum of its ribbons, stubs and fan, `k ×` its stem’s count wide: the genome, a pair, one, a hub', async () => {
    const specs = [ringSpec(), ringSpec({ include: ['chr2', 'chr5'] }), ringSpec({ include: ['chr1'] }), hubSpec()];
    for (const [i, spec] of specs.entries()) {
      const l = lay(await (await settle(spec)).frame({ size: Size }))!;
      expect(necksAdd(l), `spec ${i}`).toEqual([]);
      const perPx = l.frame.perPx(l.prep!.config.stemDepth * l.frame.room);
      let widths = 0;
      for (let st = 0; st < l.ns; st++) {
        if (!l.stemPlaced[st] || l.k * l.stemCount[st] < 1) continue;
        expect(2 * l.stemNeck[st], `spec ${i}`).toBeCloseTo(l.k * l.stemCount[st] * perPx, 9);
        widths++;
      }
      expect(widths, `spec ${i}`).toBeGreaterThan(0);
      // the widest stem in view is `maxWidthPx` wide
      let max = 0;
      for (let st = 0; st < l.ns; st++) if (l.stemPlaced[st]) max = Math.max(max, l.stemCount[st]);
      expect(l.k * max, `spec ${i}`).toBeCloseTo(l.prep!.config.maxWidthPx, 9);
    }
  });

  it('a hub of 300 partners keeps its width and count: its weak partners its fan', async () => {
    const l = lay(await (await settle(hubSpec())).frame({ size: Size }))!;
    let hub = -1;
    for (let st = 0; st < l.ns; st++) if (hub < 0 || l.stemCount[st] > l.stemCount[hub]) hub = st;
    expect(l.stemCount[hub]).toBe(340);
    expect(necksAdd(l)).toEqual([]);
    expect(l.stemFan[hub]).toBeGreaterThan(0);
    // the strongest ribbon is the 40 to chr3, drawn
    expect(l.ribCount[l.strongest]).toBeGreaterThanOrEqual(40);
  });

  it('a ribbon’s ends are inside its stems’ necks; a stem’s count the pairs leaving it', async () => {
    const s = await settle(ringSpec());
    const l = lay(await s.frame({ size: Size }))!;
    const inNeck = (st: number, t0: number, t1: number) => {
      const mid = (l.stemT0[st] + l.stemT1[st]) / 2;
      return t0 >= mid - l.stemNeck[st] - 1e-9 && t1 <= mid + l.stemNeck[st] + 1e-9;
    };
    for (let r = 0; r < l.nr; r++) {
      if (!l.ribDrawn[r]) continue;
      expect(inNeck(l.ribA[r], l.ribA0[r], l.ribA1[r])).toBe(true);
      expect(inNeck(l.ribB[r], l.ribB0[r], l.ribB1[r])).toBe(true);
    }
    const count = new Float64Array(l.ns);
    for (const r of rows(l)) {
      if (r.a === r.b || (!l.stemPlaced[r.a] && !l.stemPlaced[r.b])) continue;
      if (l.stemPlaced[r.a]) count[r.a] += r.w;
      if (l.stemPlaced[r.b]) count[r.b] += r.w;
    }
    for (let st = 0; st < l.ns; st++) if (l.stemPlaced[st]) expect(l.stemCount[st]).toBe(count[st]);
  });

  it('stubs only for ends with no place, one a stem a partner, in the partner’s hue', async () => {
    const s = await settle(ringSpec({ include: ['chr2'] }));
    const l = lay(await s.frame({ size: Size }))!;
    expect(l.nb).toBeGreaterThan(0);
    const seen = new Set<string>();
    for (let b = 0; b < l.nb; b++) {
      const st = l.stubStem[b];
      expect(l.stemPlaced[st]).toBe(1);
      expect(l.names[l.stemSystem[st]]).toBe('chr2');
      expect(l.names[l.stubPartner[b]]).not.toBe('chr2');
      const key = `${st}|${l.stubPartner[b]}`;
      expect(seen.has(key)).toBe(false);
      seen.add(key);
      if (l.stubPoly[b] >= 0)
        expect(l.polyColor[l.stubPoly[b]] >>> 8).toBe(l.hues[l.stubPartner[b] % l.hues.length] >>> 8);
    }
    // a ring of every chromosome: no stub (every end has a place)
    const all = lay(await (await settle(ringSpec())).frame({ size: Size }))!;
    expect(all.nb).toBe(0);
    for (let i = 0; i < all.polys.n; i++) expect(all.polyKind[i]).not.toBe(Kind.stub);
  });
});

/** A ring of explicit blocks, chr2 reversed, its links read from BEDPE. */
function reversedSpec() {
  const chroms = toyChromosomes(6);
  const links = toyLinks(chroms, 1500, 5);
  const b = createSeqSpec({ title: 'a reversed block', extensions: [rondo] });
  const g = b.data.assembly('toy', assembly('Toy', chroms));
  const pairs = b.data.pairs('links', {
    a: g,
    source: b.data.source('src', { format: 'bedpe', value: linksBedpe(chroms, links) as never }),
  } as never);
  const seq = (name: string) => `/data/toy/(sequence :name ${name})`;
  b.views.view('v', (v) =>
    v.section('ring', { layout: 'rings', center: 'none' } as never, (s) => {
      for (const c of chroms)
        s.block(c.name, {
          of: seq(c.name),
          label: c.name,
          ...(c.name === 'chr2' ? { direction: 'reversed' } : {}),
        } as never);
      s.track('links', { height: 4 }).layer('ribbons', {
        type: 'rondo:ribbons',
        fallback: 'pairs',
        data: pairs,
      } as never);
    }),
  );
  return b.build();
}

describe('the stubs’ labels', () => {
  it('never overlap: the strongest labelled first, a label only where its box is free', async () => {
    const text = bundledTextMetrics();
    // one chromosome of a genome with many partners: a stub a stem a partner, many tips near each other
    const chroms = toyChromosomes(12);
    const links = toyLinks(chroms, 3000, 9, 0.2);
    for (const include of [['chr1'], ['chr3'], ['chr2', 'chr7']]) {
      const l = lay(await (await settle(ringSpec({ chroms, links, include }))).frame({ size: Size }))!;
      const size = Number.parseFloat(l.font);
      const boxes: [number, number, number, number, string][] = [];
      for (let s = 0; s < l.nb; s++) {
        const t = l.stubLabel[s];
        if (!t) continue;
        const w = text.measure(l.font, t);
        expect(l.stubLabelW[s]).toBeCloseTo(w, 6);
        const y = l.stubLabelY[s] + 4;
        boxes.push([l.stubLabelX[s] - w / 2, y - size, l.stubLabelX[s] + w / 2, y + size * 0.3, t]);
      }
      expect(boxes.length, include.join()).toBeGreaterThan(0);
      for (let i = 0; i < boxes.length; i++)
        for (let j = i + 1; j < boxes.length; j++) {
          const [a, b] = [boxes[i], boxes[j]];
          const overlap = a[0] < b[2] && a[2] > b[0] && a[1] < b[3] && a[3] > b[1];
          expect(overlap, `${include.join()}: ${a[4]} and ${b[4]}`).toBe(false);
        }
      // the strongest stub is labelled
      let strongest = -1;
      for (let s = 0; s < l.nb; s++)
        if (l.stubPoly[s] >= 0 && (strongest < 0 || l.stubCount[s] > l.stubCount[strongest])) strongest = s;
      expect(l.stubLabel[strongest]).not.toBe('');
    }
  });
});

describe('drawing against the data: a reversed block, a tile edge', () => {
  it('a reversed block: its stems run counterclockwise (`stemDir` −1), every pair in one place, the necks add', async () => {
    const l = lay(await (await settle(reversedSpec())).frame({ size: Size }))!;
    const chr2 = l.names.indexOf('chr2');
    const mine: number[] = [];
    for (let s = 0; s < l.ns; s++) if (l.stemPlaced[s] && l.stemSystem[s] === chr2) mine.push(s);
    expect(mine.length).toBeGreaterThan(2);
    for (const s of mine) expect(l.stemDir[s]).toBe(-1);
    for (let s = 0; s < l.ns; s++) if (l.stemPlaced[s] && l.stemSystem[s] !== chr2) expect(l.stemDir[s]).toBe(1);
    // later positions at earlier places on chr2's arc
    mine.sort((a, b) => l.stemBin[a] - l.stemBin[b]);
    for (let i = 1; i < mine.length; i++) expect(l.stemT0[mine[i]]).toBeLessThan(l.stemT0[mine[i - 1]]);
    const a = accounted(l);
    expect(a.total).toBeGreaterThan(0);
    expect(a.wrong).toEqual([]);
    expect(necksAdd(l)).toEqual([]);
  });

  it('pairs across a tile edge: read from both tiles, each counted once', async () => {
    // chr1 of human size alone on the ring: a tile at its level is shorter than it
    const chroms = genomeChromosomes(6);
    const first = await settle(ringSpec({ chroms, links: toyLinks(chroms, 6000), include: ['chr1'] }));
    const l0 = lay(await first.frame({ size: Size }))!;
    const tile = first.store.settings.tileBins * l0.level;
    const edge = Math.floor(chroms[0].length / 2 / tile) * tile;
    expect(edge).toBeGreaterThan(0);
    // 120 links from either side of chr1's tile edge to chr3, 60 a side, beside the default's
    const links: Link[] = [...toyLinks(chroms, 6000)];
    for (let i = 0; i < 120; i++)
      links.push({ a: 0, b: 2, as: edge + (i < 60 ? -600 - i : i), bs: 30_000_000 + i * 10, type: 'contact' });
    const l = lay(await (await settle(ringSpec({ chroms, links, include: ['chr1'] }))).frame({ size: Size }))!;
    expect(l.level).toBe(l0.level);
    expect(l.items).toBe(false);
    const a = accounted(l);
    expect(a.wrong).toEqual([]);
    // each stem either side of the edge counts the links with an end in it (an end in each of two stems: both)
    const stem = (pos: number) => Math.floor(pos / l.stemLevel);
    const sides = [stem(edge - 1), stem(edge)];
    for (const side of sides) {
      let expected = 0;
      for (const k of links) {
        const inA = k.a === 0 && stem(k.as) === side;
        const inB = k.b === 0 && stem(k.bs) === side;
        if (inA && inB) continue;
        if (inA || inB) expected++;
      }
      let counted = -1;
      for (let s = 0; s < l.ns; s++)
        if (l.names[l.stemSystem[s]] === 'chr1' && l.stemBin[s] === side) counted = l.stemCount[s];
      expect(counted).toBe(expected);
    }
  });
});

describe('the toy links', () => {
  it('are seeded: the same links twice', () => {
    const c = toyChromosomes();
    expect(toyLinks(c, 50, 3)).toEqual(toyLinks(c, 50, 3));
  });
});
