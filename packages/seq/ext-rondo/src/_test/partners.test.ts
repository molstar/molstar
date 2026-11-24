/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { RondoLay } from '../lay.js';
import { capped } from '../menu.js';
import { partnerChromosomes } from '../partners.js';

// A stem's partners by chromosome on a synthetic lay (no layout): every partner listed, the menu's cap its own.

/** A lay of one stem (0, on system 0) with a ribbon to each of `n` other systems, `k` stems a system, and links within. */
function hub(n: number, k: number): RondoLay {
  const lay = new RondoLay();
  lay.names = Array.from({ length: n + 1 }, (_, i) => `chr${i + 1}`);
  const stems = 1 + n * k;
  lay.ns = stems;
  lay.stemSystem = new Int32Array(stems);
  lay.withinLinks = new Float64Array(stems);
  lay.stemWithin = new Float64Array(stems);
  lay.withinLinks[0] = 5;
  lay.stemWithin[0] = 5;
  lay.nr = n * k;
  lay.ribA = new Int32Array(lay.nr);
  lay.ribB = new Int32Array(lay.nr);
  lay.ribCount = new Float64Array(lay.nr);
  lay.ribLinks = new Float64Array(lay.nr);
  lay.ribStrand = new Int8Array(lay.nr);
  for (let s = 1; s < stems; s++) {
    const sys = 1 + Math.floor((s - 1) / k);
    lay.stemSystem[s] = sys;
    const r = s - 1;
    lay.ribA[r] = 0;
    lay.ribB[r] = s;
    // more links to the earlier systems, a stem's pieces by count
    lay.ribCount[r] = lay.ribLinks[r] = 1000 - 10 * sys - (s % k);
    lay.ribStrand[r] = s % 2 ? 1 : -1;
  }
  return lay;
}

describe('partnerChromosomes', () => {
  it('lists every partner chromosome and every piece, the stem’s own chromosome last (the keys walk them all)', () => {
    const rows = partnerChromosomes(hub(30, 12), 0);
    expect(rows).toHaveLength(31);
    expect(rows.slice(0, 30).map((r) => r.system)).toEqual(Array.from({ length: 30 }, (_, i) => i + 1));
    for (const r of rows.slice(0, 30)) expect(r.pieces).toHaveLength(12);
    expect(rows[30]).toMatchObject({ system: 0, within: true, links: 5 });
  });
});

describe('the stem menu’s cap', () => {
  const rows = (n: number) => Array.from({ length: n }, (_, i) => ({ label: `r${i}`, command: 'c', args: { i } }));

  it('at most 8 rows: past it the first 8 and “n more”, its submenu the rest in order', () => {
    expect(capped(rows(8)).map((r) => r.label)).toEqual(rows(8).map((r) => r.label));
    const c = capped(rows(12));
    expect(c.map((r) => r.label)).toEqual(['r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', '4 more']);
    expect(c[8].items!.map((r) => r.label)).toEqual(['r8', 'r9', 'r10', 'r11']);
    expect(c[8].args).toEqual({ i: 8 });
  });

  it('every row reachable: a long list nests “n more” in turn, none left out', () => {
    const flat = (xs: readonly { label: string; items?: readonly { label: string }[] }[]): string[] =>
      xs.flatMap((x) => (x.items && / more$/.test(x.label) ? flat(x.items as never) : [x.label]));
    const c = capped(rows(30));
    expect(c).toHaveLength(9);
    expect(flat(c)).toEqual(rows(30).map((r) => r.label));
  });
});
