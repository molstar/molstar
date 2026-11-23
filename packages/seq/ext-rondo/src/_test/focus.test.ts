/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { DefaultRondoConfig } from '../config.js';
import { rondoWith } from '../index.js';
import { foldedOf, Kind, type RondoLay } from '../lay.js';
import { describeRondo } from '../words.js';
import { type Link, lay, ringSpec, Size, settle, toyChromosomes } from './helpers.js';

// On one chromosome or a pair, the information first (_spec/README.md, Stems, Stubs): the trumpets lighter
// (`rondo.focusStems`), a stem's stubs under `rondo.stubMin` links folded into one "other chromosomes" stub, its count
// their sum (a neck still the sum of its parts), its card listing them.

const chroms = toyChromosomes(6);

/** One place of chr1 (one stem): 30 links to chr2, one each to chr3, chr4 and chr5, and 5 to chr6; a few within. */
function links(): Link[] {
  const out: Link[] = [];
  const at = (_i: number) => 100_000;
  for (let i = 0; i < 30; i++) out.push({ a: 0, b: 1, as: at(i), bs: 2_000_000 + i * 10_000, type: 'contact' });
  for (const b of [2, 3, 4]) out.push({ a: 0, b, as: at(b), bs: 1_000_000, type: 'contact' });
  for (let i = 0; i < 5; i++) out.push({ a: 0, b: 5, as: at(i + 10), bs: 3_000_000 + i * 10_000, type: 'contact' });
  for (let i = 0; i < 4; i++) out.push({ a: 0, b: 0, as: at(i), bs: at(i) + 400_000, type: 'contact' });
  return out;
}

const alphaOf = (c: number) => c & 0xff;

async function ring(include?: string[], extension = rondoWith()) {
  const s = await settle(ringSpec({ chroms, links: links(), held: true, ...(include ? { include } : {}) }), [
    extension,
  ]);
  return lay(await s.frame({ size: Size, dpr: 2 }))!;
}

/** The trumpets' alphas. */
const trumpets = (l: RondoLay) => {
  const out = new Set<number>();
  for (let i = 0; i < l.polys.n; i++) if (l.polyKind[i] === Kind.trumpet) out.add(alphaOf(l.polyColor[i]));
  return [...out];
};

describe('one chromosome or a pair: the information first', () => {
  it('the trumpets lighter on one chromosome or a pair, as they were on the genome', async () => {
    const op = Math.round(DefaultRondoConfig.focusStems * 0.35 * 255);
    const whole = await ring();
    expect(whole.focus).toBe(false);
    expect(trumpets(whole)).toEqual([Math.round(0.35 * 255)]);
    for (const include of [['chr1'], ['chr1', 'chr2']]) {
      const l = await ring(include);
      expect(l.focus, include.join()).toBe(true);
      expect(trumpets(l), include.join()).toEqual([op]);
    }
  });

  it('stubs under `rondo.stubMin` fold into one “other chromosomes” stub a stem, its count their sum', async () => {
    const l = await ring(['chr1']);
    const names = (b: number) => (l.stubPartner[b] < 0 ? 'other' : l.names[l.stubPartner[b]]);
    const stubs = Array.from({ length: l.nb }, (_, b) => `${names(b)} ${l.stubCount[b]}`).sort();
    expect(new Set(Array.from({ length: l.nb }, (_, b) => l.stubStem[b])).size).toBe(1);
    // chr2's 30 and chr6's 5 stay; chr3, chr4 and chr5's one each fold
    expect(stubs).toEqual(['chr2 30', 'chr6 5', 'other 3']);
    const folded = Array.from({ length: l.nb }, (_, b) => b).find((b) => l.stubPartner[b] < 0)!;
    expect(foldedOf(l, folded).map((f) => [l.names[f.partner], f.count])).toEqual([
      ['chr3', 1],
      ['chr4', 1],
      ['chr5', 1],
    ]);
    // the neck the sum of its parts: the stem's count its stubs' (the folded one's the sum) and its ribbons'
    const stem = l.stubStem[folded];
    let sum = 0;
    for (let b = 0; b < l.nb; b++) if (l.stubStem[b] === stem) sum += l.stubCount[b];
    for (let r = 0; r < l.nr; r++) if (l.ribA[r] === stem || l.ribB[r] === stem) sum += l.ribCount[r];
    expect(sum).toBe(l.stemCount[stem]);
    // its label and its card
    if (l.stubLabel[folded]) expect(l.stubLabel[folded]).toBe('other chromosomes · 3');
    const card = describeRondo({ element: { key: l.stubKey[folded] } } as never, l)!;
    expect(card[0]).toMatch(/→ other chromosomes$/);
    expect(card[2]).toBe('Partners: chr3 (1), chr4 (1), chr5 (1)');
  });

  it('`rondo.stubMin` 1 folds none; a lone small stub is never folded', async () => {
    const none = await ring(['chr1'], rondoWith({ stubMin: 1 }));
    expect(Array.from({ length: none.nb }, (_, b) => none.stubPartner[b]).every((p) => p >= 0)).toBe(true);
    expect(none.nb).toBe(5);
    // at 6, chr6's 5 joins them; chr2's 30 stays alone
    const more = await ring(['chr1'], rondoWith({ stubMin: 6 }));
    const folded = Array.from({ length: more.nb }, (_, b) => b).find((b) => more.stubPartner[b] < 0)!;
    expect(more.stubCount[folded]).toBe(8);
    expect(more.nb).toBe(2);
  });
});
