/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Drawing against the data (_spec/README.md, Stems, Ribbons, Stubs): every pair read with an end in view is in exactly
// one ribbon, stub or fan (or within its stem); a neck is `k ×` its stem's count and the sum of its ribbons, stubs and
// fan (a hub of 300 partners too); a stem's count the pairs with an end in it that leave it; a ribbon's ends inside
// its stems' necks; stubs only for ends with no place, one a stem a partner.

import { describe, expect, it } from 'vitest';
import { Kind, type RondoLay } from '../lay.js';
import { lay, type Link, ringSpec, Size, settle, toyChromosomes, toyLinks } from './helpers.js';

/** Each row's two stems (by the lay's stem keys), its weight; rows the filter leaves out skipped. */
function rows(l: RondoLay) {
  const read = l.read!;
  const d = read.data!;
  const index = new Map<string, number>();
  for (let s = 0; s < l.ns; s++) index.set(`${l.stemSystem[s]}:${l.stemBin[s]}`, s);
  const sysA = (k: number) => (d.aUnit ? d.aUnit[k] : 0);
  const sysB = (k: number) => (d.bUnit ? d.bUnit[k] : l.systems.length > 1 ? 1 : 0);
  const out: { a: number; b: number; w: number }[] = [];
  for (let k = 0; k < d.n; k++) {
    if (!read.shown[k]) continue;
    const a = index.get(`${sysA(k)}:${Math.floor(d.aStart[k] / l.stemLevel)}`)!;
    const b = index.get(`${sysB(k)}:${Math.floor(d.bStart[k] / l.stemLevel)}`)!;
    out.push({ a, b, w: read.weight ? read.weight[k] : 1 });
  }
  return out;
}

function accounted(l: RondoLay) {
  let total = 0;
  for (const r of rows(l)) if (l.stemPlaced[r.a] || l.stemPlaced[r.b]) total += r.w;
  let drawn = 0;
  let fanned = 0;
  for (let r = 0; r < l.nr; r++) l.ribDrawn[r] ? (drawn += l.ribCount[r]) : (fanned += l.ribCount[r]);
  let stubs = 0;
  for (let b = 0; b < l.nb; b++) stubs += l.stubCount[b];
  let within = 0;
  for (let s = 0; s < l.ns; s++) within += l.stemWithin[s];
  return { total, drawn, fanned, stubs, within };
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

describe('drawing against the data', () => {
  it('every pair with an end in view is in exactly one ribbon, stub or fan, or within its stem (EC-171)', async () => {
    for (const o of [{}, { held: true }, { include: ['chr1'] }, { include: ['chr2', 'chr5'] }]) {
      const s = await settle(ringSpec(o));
      const l = lay(await s.frame({ size: Size }))!;
      const a = accounted(l);
      expect(a.total).toBeGreaterThan(0);
      expect(a.drawn + a.fanned + a.stubs + a.within).toBeCloseTo(a.total, 6);
    }
  });

  it('a neck is the sum of its ribbons, stubs and fan, `k ×` its stem’s count wide', async () => {
    const s = await settle(ringSpec());
    const l = lay(await s.frame({ size: Size }))!;
    expect(necksAdd(l)).toEqual([]);
    const perPx = l.frame.perPx(l.prep!.config.stemDepth * l.frame.room);
    for (let st = 0; st < l.ns; st++) {
      if (!l.stemPlaced[st] || l.k * l.stemCount[st] < 1) continue;
      expect(2 * l.stemNeck[st]).toBeCloseTo(l.k * l.stemCount[st] * perPx, 9);
    }
    // the widest stem in view is `maxWidthPx` wide
    let max = 0;
    for (let st = 0; st < l.ns; st++) if (l.stemPlaced[st]) max = Math.max(max, l.stemCount[st]);
    expect(l.k * max).toBeCloseTo(l.prep!.config.maxWidthPx, 9);
  });

  it('a hub of 300 partners keeps its width and count: its weak partners its fan', async () => {
    const chroms = toyChromosomes();
    const links: Link[] = [];
    // chr1 at 1 Mb linked once each to 300 places spread over the other chromosomes, and 40 times to chr3 at 2 Mb
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
    const s = await settle(ringSpec({ chroms, links, held: true }));
    const l = lay(await s.frame({ size: Size }))!;
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

describe('the toy links', () => {
  it('are seeded: the same links twice', () => {
    const c = toyChromosomes();
    expect(toyLinks(c, 50, 3)).toEqual(toyLinks(c, 50, 3));
  });
});
