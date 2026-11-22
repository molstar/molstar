/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { openSeq, type PairsList } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { SyntenyA, SyntenyB, SyntenyBlocks, syntenyPairs, syntenySource } from '../regions/synteny.js';

const signal = new AbortController().signal;

describe('the generated synteny (Features/Region reads; EC-208)', () => {
  it('is about 6,000 blocks inside their chromosomes, deterministic, mostly colinear with reversed ones', () => {
    const t = syntenyPairs();
    expect(Math.abs(t.aStart.length - SyntenyBlocks)).toBeLessThan(200);
    expect(Array.from(syntenyPairs().bStart.subarray(0, 40))).toEqual(Array.from(t.bStart.subarray(0, 40)));
    const lengthA = new Map(SyntenyA.map((c) => [c.name, c.length]));
    const lengthB = new Map(SyntenyB.map((c) => [c.name, c.length]));
    let reversed = 0;
    for (let i = 0; i < t.aStart.length; i++) {
      expect(t.aStart[i]).toBeLessThan(t.aEnd[i]);
      expect(t.bStart[i]).toBeLessThan(t.bEnd[i]);
      expect(t.aEnd[i]).toBeLessThanOrEqual(lengthA.get(t.aSystem!.get(i)!)!);
      expect(t.bEnd[i]).toBeLessThanOrEqual(lengthB.get(t.bSystem!.get(i)!)!);
      if (t.bStrand![i] < 0) reversed++;
    }
    expect(reversed).toBeGreaterThan(50);
    expect(reversed).toBeLessThan(t.aStart.length / 4);
  });

  it('is sided: each side answers the pairs with an end on its system, a pair found from either end', async () => {
    const data = syntenySource();
    expect(data.coverage.systems.filter((s) => s.side === 'a').map((s) => s.name)).toEqual(SyntenyA.map((c) => c.name));
    expect(data.coverage.systems.filter((s) => s.side === 'b').map((s) => s.name)).toEqual(SyntenyB.map((c) => c.name));
    const src = openSeq(data);
    const all = syntenyPairs();
    const read = async (system: string, length: number, side: 'a' | 'b') => {
      const a = await src.read({ system, start: 0, end: length, level: 1, shapes: ['pairs'], side }, signal);
      if (!a.ok) throw new Error(a.issue.message);
      return a.tables.pairs as PairsList;
    };
    // the rows with an end on chrB2 are found by side b, the same rows by side a on the A chromosomes holding their other end
    const onB2 = await read('chrB2', SyntenyB[1].length, 'b');
    let expected = 0;
    for (let i = 0; i < all.aStart.length; i++) if (all.bSystem!.get(i) === 'chrB2') expected++;
    expect(onB2.aStart.length).toBe(expected);
    let fromA = 0;
    for (const c of SyntenyA) {
      const t = await read(c.name, c.length, 'a');
      for (let i = 0; i < t.aStart.length; i++) if (t.bSystem!.get(i) === 'chrB2') fromA++;
    }
    expect(fromA).toBe(expected);
  });

  it('unsided is one list of systems, each answering the pairs with an end on it either way', async () => {
    const data = syntenySource(undefined, undefined, { sided: false });
    expect(data.coverage.systems.every((c) => c.side === undefined)).toBe(true);
    const src = openSeq(data);
    const all = syntenyPairs();
    for (const c of [SyntenyA[1], SyntenyB[1]]) {
      const a = await src.read({ system: c.name, start: 0, end: c.length, level: 1, shapes: ['pairs'] }, signal);
      if (!a.ok) throw new Error(a.issue.message);
      let expected = 0;
      for (let i = 0; i < all.aStart.length; i++)
        if (all.aSystem!.get(i) === c.name || all.bSystem!.get(i) === c.name) expected++;
      expect((a.tables.pairs as PairsList).aStart.length).toBe(expected);
    }
  });
});
