/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Properties (_spec/README.md, Stems, Ribbons, Hits): stems a partition of the non-empty bins, a level `2^j` × the
// read's, the same at DPR 1, 2 and in the SVG; items read only at or under `rondo.maxRibbons`; the twist by geometry
// the rule over every direction combination; the closed hit the drawn cover's topmost polygon (sampled); held and
// region pairs the same stems.

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { describe, expect, it } from 'vitest';
import { crossed, ribbonCorners, segmentsCross } from '../geometry.js';
import { rondo } from '../index.js';
import type { RondoLay } from '../lay.js';
import { lay, ringSpec, rng, Size, settle } from './helpers.js';

/** The placed stems by key: system, bin, count. */
const stems = (l: RondoLay) => {
  const out: string[] = [];
  for (let s = 0; s < l.ns; s++)
    if (l.stemPlaced[s] && l.stemCount[s] > 0)
      out.push(`${l.names[l.stemSystem[s]]}:${l.stemBin[s]}:${l.stemCount[s]}`);
  return out.sort();
};

describe('stems', () => {
  it('are a partition of the non-empty bins: disjoint ranges a system, the level 2^j × the read’s', async () => {
    const l = lay(await (await settle(ringSpec())).frame({ size: Size }))!;
    const j = Math.log2(l.stemLevel / l.level);
    expect(Number.isInteger(j) && j >= 0).toBe(true);
    const by = new Map<number, [number, number][]>();
    for (let s = 0; s < l.ns; s++) {
      const r = by.get(l.stemSystem[s]) ?? [];
      r.push([l.stemStart[s], l.stemEnd[s]]);
      by.set(l.stemSystem[s], r);
    }
    for (const ranges of by.values()) {
      ranges.sort((a, b) => a[0] - b[0]);
      for (let i = 1; i < ranges.length; i++) expect(ranges[i][0]).toBeGreaterThanOrEqual(ranges[i - 1][1]);
      for (const [a, b] of ranges) {
        expect(a % l.stemLevel).toBe(0);
        expect(b).toBeGreaterThan(a);
      }
    }
    // every non-empty bin (an end of a pair read) in exactly one stem of its system
    const read = l.read!;
    const d = read.data!;
    const sysOf = (k: number, end: 'a' | 'b') =>
      end === 'a' ? (d.aUnit ? d.aUnit[k] : 0) : d.bUnit ? d.bUnit[k] : l.systems.length > 1 ? 1 : 0;
    let ends = 0;
    for (let k = 0; k < d.n; k++) {
      if (!read.shown[k]) continue;
      for (const end of ['a', 'b'] as const) {
        const sys = sysOf(k, end);
        const pos = end === 'a' ? d.aStart[k] : d.bStart[k];
        let holders = 0;
        for (let s = 0; s < l.ns; s++)
          if (l.stemSystem[s] === sys && pos >= l.stemStart[s] && pos < l.stemEnd[s]) holders++;
        expect(holders, `row ${k} ${end}`).toBe(1);
        ends++;
      }
    }
    expect(ends).toBeGreaterThan(0);
  });

  it('are the same at DPR 1 and 2, and in the SVG (EC-108)', async () => {
    const s = await settle(ringSpec({ include: ['chr1', 'chr2'] }));
    const f1 = await s.frame({ size: Size, dpr: 1 });
    const one = lay(f1)!;
    const a = { level: one.stemLevel, stems: stems(one) };
    const two = lay(await s.frame({ size: Size, dpr: 2 }))!;
    expect({ level: two.stemLevel, stems: stems(two) }).toEqual(a);
    // the SVG: the composition's export draws the canvas's paths (a polygon a trumpet, ribbon, stub and fan)
    const ctx = createSeqContext({ extensions: [rondo], registry: [{ kind: 'format-runtime', runtime: BedRuntime }] });
    const svg = await ctx.renderSvg(ringSpec({ include: ['chr1', 'chr2'] }), { size: Size });
    expect(svg.issues.filter((i) => i.severity === 'error')).toEqual([]);
    expect(svg.lists.passes.data.paths.n).toBe(f1.lists.passes.data.paths.n);
    expect(svg.svg.split('<path').length - 1).toBeGreaterThanOrEqual(one.polys.n - 1);
  });

  it('held and region pairs give the same stems', async () => {
    const region = lay(await (await settle(ringSpec({ include: ['chr3'] }))).frame({ size: Size }))!;
    const held = lay(await (await settle(ringSpec({ include: ['chr3'], held: true }))).frame({ size: Size }))!;
    expect([region.items, held.items]).toEqual([true, true]);
    expect(held.stemLevel).toBe(region.stemLevel);
    expect(stems(held)).toEqual(stems(region));
  });
});

describe('the items gate', () => {
  it('items are read only at or under `rondo.maxRibbons`: bins past it, never culled below it (EC-169)', async () => {
    const fit = lay(await (await settle(ringSpec())).frame({ size: Size }))!;
    expect(fit.items).toBe(false);
    let n = 0;
    for (let k = 0; k < fit.read!.data!.n; k++) n += fit.read!.weight![k];
    expect(n).toBeGreaterThan(fit.prep!.config.maxRibbons);
    const one = lay(await (await settle(ringSpec({ include: ['chr4'] }))).frame({ size: Size }))!;
    expect(one.items).toBe(true);
    // under the gate every ribbon is drawn unless under a device px: the bound culls none
    let drawn = 0;
    for (let r = 0; r < one.nr; r++) drawn += one.ribDrawn[r];
    expect(drawn).toBeLessThanOrEqual(one.prep!.config.maxRibbons);
    for (let r = 0; r < one.nr; r++) if (!one.ribDrawn[r]) expect(one.k * one.ribCount[r]).toBeLessThan(1);
  });
});

describe('the twist (EC-206)', () => {
  it('drawn equals the rule over every direction combination: crossed exactly when strand × dirA × dirB is −1', () => {
    const out = new Float64Array(4);
    // two arcs facing each other across a ring: A at the top, B at the bottom
    const pt = (t: number) => [Math.sin(t), -Math.cos(t)] as const;
    for (const strand of [-1, 0, 1])
      for (const dirA of [-1, 1])
        for (const dirB of [-1, 1]) {
          ribbonCorners(-0.1, 0.1, dirA, Math.PI - 0.1, Math.PI + 0.1, dirB, strand, out);
          const [x0, y0] = pt(out[0]);
          const [x1, y1] = pt(out[1]);
          const [x2, y2] = pt(out[2]);
          const [x3, y3] = pt(out[3]);
          // the two sides (A → B and B → A) cross through the hole exactly when the ribbon is crossed
          expect(segmentsCross(x0, y0, x1, y1, x3, y3, x2, y2)).toBe(crossed(strand, dirA, dirB));
          expect(crossed(strand, dirA, dirB)).toBe(strand * dirA * dirB === -1);
        }
  });
});

describe('the hit', () => {
  it('the closed hit is the drawn cover’s topmost polygon, sampled', async () => {
    const f = await (await settle(ringSpec())).frame({ size: Size });
    const l = lay(f)!;
    const type = f.layout.sections[0].sectionLays.find((x) => x.type.name === 'rondo:ribbons')!.type;
    const r = rng(5);
    const at = new Float64Array(2);
    let hits = 0;
    for (let i = 0; i < 3000; i++) {
      at[0] = l.frame.cx + (r() * 2 - 1) * l.frame.r;
      at[1] = l.frame.cy + (r() * 2 - 1) * l.frame.r;
      // the brute force: the last drawn polygon (not the hole) containing the point
      let top = -1;
      for (let j = 0; j < l.polys.n; j++) {
        const p = l.order[j];
        if (p !== l.hole && l.polys.contains(p, at[0], at[1])) top = p;
      }
      const h = type.hit(l, at);
      expect(h).toBe(top < 0 ? undefined : l.hits[top]);
      if (h) hits++;
    }
    expect(hits).toBeGreaterThan(100);
  });
});
