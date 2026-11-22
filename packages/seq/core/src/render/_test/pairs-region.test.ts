/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Pairs read by region (_spec/README.md "Render", The pairs layer): zoomed out, bin pairs each an arc stroked by count
// and hit as its bin; zoomed in past the gate, the items themselves; every bin pair in view its arc (EC-208, EC-171).

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import { binPairs, type PairsList, type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, pairsBinHover } from '../index.js';
import type { PairsLay } from '../pairs.js';
import { settle } from './helpers.js';
import { LinksLength as Length, linksSource, linksSpec as spec } from './region-render-fixtures.js';

const Size = { width: 720, height: 300 };
const layOf = (f: { layout: { sections: readonly { sectionLays: readonly { lay: unknown }[] }[] } }) =>
  f.layout.sections[0].sectionLays[0].lay as PairsLay;

describe('pairs by region', () => {
  it('zoomed out: bin pairs, each an arc stroked by its count and hit as its bin (Default 2)', async () => {
    const s = await settle(spec(6000), [SeqRuntime]);
    const f = await s.frame({ size: Size, viewport: { start: 0, end: Length } });
    const lay = layOf(f);
    const prep = lay.prep!;
    expect(prep.weight).toBeDefined();
    expect(prep.level).toBe(1024);
    // every bin pair in view has its arc (fewer than `maxArcs`), its count summed into the weights
    expect(lay.n).toBe(prep.data!.n);
    expect(Array.from(prep.weight!).reduce((a, b) => a + b, 0)).toBe(6000);
    const widths = Array.from(lay.width.subarray(0, lay.n));
    expect(Math.max(...widths)).toBeCloseTo(6);
    expect(Math.min(...widths)).toBeGreaterThanOrEqual(1);
    // the widest arc's middle hits a bin pair: its level said, its extent both bins
    const i = widths.indexOf(Math.max(...widths));
    const c = lay.curve.subarray(8 * i, 8 * i + 8);
    const x = (c[0] + 3 * c[2] + 3 * c[4] + c[6]) / 8;
    const y = (c[1] + 3 * c[3] + 3 * c[5] + c[7]) / 8;
    const h = hitView(f.layout, x, y);
    expect(h?.bin).toBe(1024);
    expect(String(h?.element?.key)).toMatch(/^bin:/);
    // its hover: two bin ranges, its count, its level, never a pair
    const k = prep.index.get(h!.element!.key)!;
    const card = pairsBinHover(lay, h!.element!.key, () => 'chr')!;
    expect(card.title).toMatch(/^chr:[\d,]+–[\d,]+ ↔ chr:[\d,]+–[\d,]+$/);
    const n = prep.weight![k];
    expect(card.rows.find((r) => r[0] === 'Pairs')?.[1]).toMatch(new RegExp(`^${n} pairs?( \\(.+\\))?$`));
    expect(card.rows).toContainEqual(['Bins', '1 kb bins']);
    expect(card.text).toMatch(/zoom in to see them$/);
  });

  it('zoomed in under the gate: the items themselves, with no bins (Default 1)', async () => {
    const s = await settle(spec(6000), [SeqRuntime]);
    // three hotspots in view: the links with an end there fit `maxArcs`
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 400_000 } });
    const lay = layOf(f);
    expect(lay.prep!.weight).toBeUndefined();
    const { list } = linksSource(6000);
    let inView = 0;
    for (let i = 0; i < list.aStart.length; i++) if (list.aStart[i] < 400_000 || list.bStart[i] < 400_000) inView++;
    expect(inView).toBeLessThanOrEqual(4096);
    // the items tile is read uncut: every row of it, each once; the arcs drawn those over the view
    expect(lay.prep!.data!.n).toBe(6000);
    expect(lay.n).toBe(inView);
  });

  it('a pan inside the items tiles draws the pairs that enter the view, nothing built again (EC-208)', async () => {
    const s = await settle(spec(4000), [SeqRuntime]);
    const at = (start: number) => ({ size: Size, viewport: { start, end: start + 400_000 } });
    const first = layOf(await s.frame(at(0)));
    const prep = first.prep;
    const lay = layOf(await s.frame(at(50_000)));
    // the same tiles read: the same data, not built again
    expect(lay.prep).toBe(prep);
    const { list } = linksSource(4000);
    // an arc is drawn unless both its ends are past one side of the view
    let shown = 0;
    let entering = 0;
    for (let i = 0; i < list.aStart.length; i++) {
      const a = list.aStart[i];
      const b = list.bStart[i];
      if ((a < 50_000 && b < 50_000) || (a >= 450_000 && b >= 450_000)) continue;
      shown++;
      if (a >= 400_000 && a < 450_000 && b >= 400_000) entering++;
    }
    expect(entering).toBeGreaterThan(0);
    // within `maxArcs`: each its own arc
    expect(shown).toBeLessThanOrEqual(4096);
    expect(lay.n).toBe(shown);
  });
});
