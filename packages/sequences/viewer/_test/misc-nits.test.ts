/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The misc gap review's nits (2026-10-14): a reversed chromosome's centromere pinches (GM-06); a data column reads a
// blank, a hex or `Infinity` as text, not a number (GM-07); a raw column's plurality past 8 colors (GM-09).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import type { SectionAxisView } from '../src/axis.js';
import { Dirty } from '../src/draw/index.js';
import { type CytoBands, ideogramNavigator, type PresetEntry } from '../src/index.js';
import { NumberColumn } from '../src/renderers/columns.js';
import { blocksDoc } from './axis-fixtures.js';
import { composed, recording } from './composer-fixtures.js';
import { directed } from './conformance-fixtures.js';

const BANDS: CytoBands = {
  c0: [
    [0, 180, 'p1', 'gneg'],
    [180, 220, 'p11', 'acen'],
    [220, 260, 'q11', 'acen'],
    [260, 400, 'q1', 'gneg'],
  ],
};
const KARYO: PresetEntry = { name: 'karyo', geometry: 'linear', render: {}, navigator: 'ideogram' };

describe("GM-06: the ideogram's centromere on a reversed chromosome", () => {
  for (const direction of [1, -1] as const)
    it(`direction ${direction}: both triangles' apexes meet at the centromere`, async () => {
      const t = await composed(directed(blocksDoc(1, 400), direction), {
        width: 840,
        height: 200,
        dpr: 1,
        registry: [{ navigators: [ideogramNavigator(BANDS)], presets: [KARYO] }],
      });
      t.ctx.ui.setStyle(t.section, 'karyo');
      await t.ctx.graph.update();
      const be = t.frame(recording(), Dirty.nav);
      const av = t.ctx.layout.axis(t.section) as SectionAxisView;
      // The centromere's middle (index 220) is shown at 220 forward, 180 reversed.
      const k = direction > 0 ? 220 : 180;
      const mid = Math.round(av.navX0s[0] + (k * (av.navX1s[0] - av.navX0s[0])) / 400);
      const c = Array.from(be.pathItems.coords.subarray(0, be.pathItems.nc));
      // Each triangle: its base's top, its apex, its base's bottom.
      expect([c[2], c[8]]).toEqual([mid, mid]);
      t.ctx.dispose();
    });
});

describe('GM-07: a number column over a string field', () => {
  it('reads decimals, and a blank, a hex or `Infinity` as text', () => {
    const strings = ['', ' ', '\t', '0x1A', 'Infinity', '12.5', ' 3 ', '-1e3', 'abc'];
    const set = { strings: { get: (k: number) => strings[k] } } as never;
    const field = { kind: 'string', data: Int32Array.from(strings, (_, k) => k) } as never;
    const ann = { count: strings.length } as never;
    const s = NumberColumn.prepare?.(field, new Int32Array(0), ann, set) as unknown as { value: Float64Array };
    expect(Array.from(s.value, (v) => (Number.isNaN(v) ? null : v))).toEqual([
      null,
      null,
      null,
      null,
      null,
      12.5,
      3,
      -1000,
      null,
    ]);
  });
});

describe("GM-09: a raw column counts a theme's colors past 8 (LC-11)", () => {
  it("a filtered protein whose columns' plurality letter is the tenth color seen (Taylor: a hue each) draws it", async () => {
    // Per 20 residues from 10 (where this view's columns of 20 start): nine of a color each, then eleven L.
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: `${'L'.repeat(10)}ACDEFGHIKL`.repeat(2000) });
    // Every position kept but one, far right: filtered, so its columns are read raw.
    const value = { start: [0, 30_001], stop: [30_000, 40_000] };
    const keep = b.data.annotation({ ref: 'keep', geometry: 'ranges', target: a, value });
    const sec = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.s', positions: keep });
    sec.block({ sequence: a });
    sec.track({ label: 'Letters' }).letters({ color: { theme: 'scheme', params: { name: 'taylor' } } } as never);
    const t = await composed(b.getDocument(), { section: 'v.s', width: 1000, height: 200 });
    const vp = t.vp;
    if (!vp) throw new Error('no viewport');
    const cells = (cell: number) => {
      vp.setView(cell, 0, 'api');
      const be = recording();
      t.frame(be, Dirty.all);
      const r = be.rectItems;
      return Array.from({ length: r.n }, (_, i) => i)
        .filter((i) => r.h[i] > 4)
        .map((i) => [r.x[i], r.color[i]]);
    };
    // L's color: residue 0's cell at 16 px a position, none of the next nine's.
    const near = cells(16);
    const L = near[0][1];
    expect(new Set(near.slice(10, 19).map((c) => c[1])).size).toBe(9);
    expect(new Set(near.slice(10, 19).map((c) => c[1])).has(L)).toBe(false);
    const far = cells(1 / 20);
    expect(far.length).toBeGreaterThan(50);
    expect(far.filter((c) => c[1] === L).length).toBeGreaterThan(0.9 * far.length);
  });
});
