/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The layout on edge geometry (the review gate's findings, 2026-10-14): stacked blocks, reversed and filtered pieces.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { Lod } from '../src/index.js';
import { DefaultRenderProps, type RenderClamp, sectionProps } from '../src/layout.js';
import { DefaultPresets } from '../src/presets.js';
import { driven } from './input-fixtures.js';
import { opened } from './layout-fixtures.js';

describe('VL-01: a stacked section anchors on the block under the point', () => {
  async function stacked() {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(50) });
    b.data.sequence({ ref: 'b', value: 'MKVLA'.repeat(2) });
    b.view({ ref: 'v', label: 'V', arrangement: 'stacked' }).track({ ref: 'v.letters' }).letters();
    const t = await opened(b.getDocument());
    const v = t.vp('v.main');
    v.zoomTo(0, 100, 300);
    return v;
  }

  it("a zoom keeps the position under the cursor, past the short block's end", async () => {
    const v = await stacked();
    const s = v.live.current;
    const before = (s.scrollX + 600) / s.cell;
    v.zoomAt(1.5, 600, 'api');
    // Within the scroll's snap to a device px (before VL-01: 243.9 → 95.9).
    expect((s.scrollX + 600) / s.cell).toBeCloseTo(before, 0);
  });

  it('the center element is on the first block reaching the center, here the long one', async () => {
    const v = await stacked();
    const io = new Float64Array(2);
    expect(v.centerElement(io)).toBe(0);
    expect(io[0]).toBeCloseTo(200, 0);
  });
});

describe('VL-02: a reversed filtered block with many breaks builds its axis', () => {
  function filtered(n: number, direction: 1 | -1) {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', alphabet: 'dna', value: 'ACGT'.repeat(n / 2) });
    const start = Array.from({ length: n }, (_, i) => 2 * i);
    b.data.annotation({ ref: 'shown', geometry: 'ranges', target: a, value: { start, stop: start.map((s) => s + 1) } });
    const v = b.view({ ref: 'v', label: 'V' });
    v.section({ ref: 'v.s', positions: 'shown', direction }).block({ sequence: a }).track().letters();
    return b.getDocument();
  }

  for (const direction of [1, -1] as const) {
    it(`150,000 single-position runs, direction ${direction}: every one shown and laid out`, async () => {
      const t = await opened(filtered(150_000, direction));
      const ax = t.vp('v.s').axis;
      expect(ax?.shown[0]).toBe(150_000);
      expect(t.ctx.layout.get('v.s')?.current.entryCount).toBeGreaterThan(0);
    }, 60_000);
  }
});

describe("VL-03: a follower's data area or DPR changing at the same width leads nothing", () => {
  function doc() {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(60) });
    const v = b.view({ ref: 'v', label: 'V' });
    v.section({ ref: 'v.lead' }).block({ sequence: a }).track().letters();
    const f = v.section({ ref: 'v.follow', follow: 'v.lead', height: 'box', style: { vertical_nav: 'overview' } });
    f.block({ sequence: a });
    for (let i = 0; i < 12; i++) f.track().letters();
    return b.getDocument();
  }
  async function moved(dpr = 1) {
    const t = await driven(doc(), { size: { 'v.lead': [416, 200], 'v.follow': [1000, 120] } });
    t.vp('v.follow').setSize(1000, 120, dpr);
    const lead = t.vp('v.lead');
    lead.zoomTo(0, 520, 570, { origin: 'pointer' });
    t.scheduler.advance(1000);
    expect(t.layout('v.follow').vnavWidth).toBeGreaterThan(0);
    return { t, lead, fol: t.vp('v.follow'), before: lead.live.current.scrollX };
  }

  it("the follower's overview leaves (taller, same width): the leader stays", async () => {
    const { t, lead, fol, before } = await moved();
    fol.setSize(1000, 2000, 1);
    t.layout('v.follow');
    t.scheduler.advance(1000);
    expect(t.layout('v.follow').vnavWidth).toBe(0);
    expect(lead.live.current.scrollX).toBe(before);
  });

  it('a DPR change of the follower (its scroll snapped again, before VL-03 a lead): the leader stays', async () => {
    // From 1.25 to 1 the follower's scroll re-snaps (LC-31: from 1.5 at an integer scroll nothing moved).
    const { t, lead, fol, before } = await moved(1.25);
    fol.setSize(1000, 120, 1);
    t.layout('v.follow');
    t.scheduler.advance(1000);
    expect(lead.live.current.scrollX).toBe(before);
  });
});

describe("VL-11: the host's clamp reaches the sizes that follow others", () => {
  it("a clamp of the derived sizes holds; a clamp of the labels' size reaches the letters' cap it gives", () => {
    const capped = sectionProps({ ...DefaultRenderProps, labelFontPx: 20 }, undefined, null, {
      max: { maxLetterPx: 14, wrappedFontPx: 10 },
    });
    expect([capped.maxLetterPx, capped.wrappedFontPx]).toEqual([14, 10]);
    const labels = sectionProps({ ...DefaultRenderProps, labelFontPx: 20 }, undefined, null, {
      max: { labelFontPx: 12 },
    });
    expect([labels.labelFontPx, labels.maxLetterPx, labels.wrappedFontPx]).toEqual([12, 12, 12]);
  });
  it("a cap's clamp reaches the text layout's letters; a `min` of a size left to follow is no cap (LC-30)", () => {
    const at = (clamp: RenderClamp) => {
      const p = sectionProps({ ...DefaultRenderProps, labelFontPx: 20 }, undefined, null, clamp);
      return [p.maxLetterPx, p.wrappedFontPx];
    };
    expect(at({ max: { maxLetterPx: 14 } })).toEqual([14, 14]);
    expect(at({ min: { maxLetterPx: 14 } })).toEqual([20, 20]);
    expect(at({ min: { wrappedFontPx: 10 } })).toEqual([20, 20]);
    // A `min` past what it follows still holds.
    expect(at({ min: { maxLetterPx: 24 } })).toEqual([24, 24]);
  });
});

describe("VL-17: the rings' LOD follows the section's `minLetterCell`", () => {
  async function ring(minLetterCell: number) {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: 's', value: 'ACGTACGTTA'.repeat(6), alphabet: 'dna', cyclic: true });
    b.view({ ref: 'v', label: 'V', style: 'circular' }).block({ sequence: s }).track({ label: 'L' }).letters();
    const t = await opened(b.getDocument(), {
      registry: [{ presets: DefaultPresets.filter((p) => p.name === 'circular') }],
      render: { minLetterCell },
    });
    const v = t.ctx.layout.axis('v.main');
    return { cell: v?.cell ?? 0, lod: v?.lod ?? -1 };
  }

  it("letters where a position's arc reaches it, cells below it", async () => {
    const at8 = await ring(8);
    expect(at8.cell).toBeGreaterThan(8);
    expect(at8.lod).toBe(Lod.letters);
    const higher = await ring(Math.ceil(at8.cell) + 1);
    expect(higher.lod).toBe(Lod.cells);
  });
});
