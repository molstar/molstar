/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Blocks of a genome in the tracks and rings layouts (_spec/README.md "Render", Blocks and Rings layout): a block's own
// direction and scale (its transform, its label's mark and `×N`), the gap share and runs of small blocks in both
// layouts (EC-102), groups of two assemblies (EC-210), the built-in fields and their block-ordered domain (EC-218),
// drawing against the data (at the fit, across a block gap, inside a run, on a reversed block, EC-206), the export's
// block attributes and its viewport by `{ system, position }`, and equal block lists that follow.

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { blockElementAt, placesOf } from '@molstar/seqx-core/view';
import { beforeAll, describe, expect, it } from 'vitest';
import {
  type BlockPlacement,
  bundledTextMetrics,
  hitView,
  LightTokens,
  PlacePoint,
  type RingPlacement,
} from '../index.js';
import { renderSvg } from '../render-svg.js';
import type { RangesLay } from '../ranges-lay.js';
import { dna, genomeSpec, MixedFeatures, mixedSpec, toyGenome, twoGenomesSpec } from './genome-fixtures.js';
import { settle } from './helpers.js';

beforeAll(async () => {
  await loadGenome();
});

const Size = { width: 900, height: 400 };

type Layout = Awaited<ReturnType<Awaited<ReturnType<typeof settle>>['frame']>>['layout'];
const placementOf = (layout: Layout) => layout.sections[0].placement as BlockPlacement;

/**
 * Hit-testing is placement's inverse: a point on a strip's data reads a unit that places back at its x on that strip,
 * and the position the strip's transform reads there is the view's (`blockPositionAt`).
 */
function inverse(layout: Layout, rng: PCG, samples = 300): number {
  const s = layout.sections[0];
  const p = s.placement as BlockPlacement;
  const blocks = s.section.blocks!;
  const at = new PlacePoint();
  const back = new PlacePoint();
  let checked = 0;
  for (let n = 0; n < samples; n++) {
    const i = rng.int() % p.strips.n;
    if (p.strips.chip[i]) continue;
    const a = p.strips.axes[i];
    const x = a.x0 + ((rng.int() % 10_000) / 10_000) * a.width * 0.999;
    const band = p.strips.firstBand(i);
    const k = i * p.strips.bands + band;
    const y = p.strips.bandY[k] + p.strips.bandH[k] / 2;
    expect(p.at([x, y], at)).toBe(true);
    if (at.kind !== 'data') continue;
    expect(at.strip).toBe(i);
    expect(p.place(at.offset, back)).toBe(true);
    expect(back.strip).toBe(i);
    expect(back.x).toBeCloseTo(x, 6);
    const block = blocks[p.strips.block[i]];
    expect(blockElementAt(block, at.offset)).toBe(a.positionAt(x));
    checked++;
  }
  return checked;
}

describe('a block’s own direction and scale (EC-206)', () => {
  it('hit-testing is the inverse of placement over mixed directions and scales (property)', async () => {
    const s = await settle(mixedSpec());
    const rng = new PCG(5);
    for (const viewport of [undefined, { start: 120, end: 520 }, { start: 380, end: 460 }]) {
      const f = await s.frame({ size: Size, ...(viewport ? { viewport } : {}) });
      expect(inverse(f.layout, rng)).toBeGreaterThan(100);
    }
  });

  it('a reversed block reads right to left, a scaled one takes its scale’s share; labels say so', async () => {
    const s = await settle(mixedSpec());
    const f = await s.frame({ size: Size });
    const p = placementOf(f.layout);
    const a = p.strips.axes;
    // block 1 [150, 400) reversed: its left edge shows 399, its right 150
    expect(a[1].dir).toBe(-1);
    expect(a[1].positionAt(a[1].x0 + 0.01)).toBe(399);
    expect(a[1].positionAt(a[1].x0 + a[1].width - 0.01)).toBe(150);
    // block 2: 100 positions at scale 3 against block 0's 150 at 1: twice its width
    expect(a[2].width / a[0].width).toBeCloseTo(2, 1);
    expect(p.labels).toEqual(['b0', 'b1 ←', 'b2 ×3', 'b3 ← ×0.5']);
  });

  it('draws against the data on a reversed and a scaled block: each feature hit at its place (EC-206)', async () => {
    const s = await settle(mixedSpec());
    const f = await s.frame({ size: Size });
    const sec = f.layout.sections[0];
    const p = sec.placement as BlockPlacement;
    const band = 1;
    const pt = new PlacePoint();
    let hits = 0;
    for (const feat of MixedFeatures) {
      const mid = (feat.start + feat.end) / 2;
      for (const place of placesOf(sec.section, `/data/s${feat.seq}`, mid)) {
        expect(p.place(place.unit, pt)).toBe(true);
        const k = pt.strip * p.strips.bands + band;
        const h = hitView(f.layout, pt.x, p.strips.bandY[k] + 6);
        expect((h?.element as { key?: string } | undefined)?.key, feat.id).toBe(feat.id);
        expect(h?.block).toBe(place.block);
        hits++;
      }
    }
    expect(hits).toBe(4);
  });

  it('draws against the data inside a run: a block read while a device pixel wide shows its feature at its place', async () => {
    const tiny = Array.from({ length: 40 }, (_, k) => ({ seq: 1, range: [k * 5, k * 5 + 5] as const }));
    const s = await settle(mixedSpec([{ seq: 0, range: [0, 400], scale: 10 }, ...tiny]));
    const f = await s.frame({ size: Size, dpr: 2 });
    const sec = f.layout.sections[0];
    const p = sec.placement as BlockPlacement;
    const [place] = placesOf(sec.section, '/data/s1', 52);
    expect(place.block).toBe(11);
    expect(p.row.run[place.block]).toBeGreaterThanOrEqual(0);
    expect(p.runN).toBe(1);
    const pt = new PlacePoint();
    expect(p.place(place.unit, pt)).toBe(true);
    const k = pt.strip * p.strips.bands + 1;
    // its features layer laid out on its own strip, the feature that covers it among its items (a bar a device px)
    const lay = sec.lays[k]!.find((l) => l.layer.address?.endsWith('/f1'))!.lay as RangesLay;
    expect(Array.from({ length: lay.items }, (_, i) => String(lay.keys[i]))).toContain('f2');
    const h = hitView(f.layout, pt.x, p.strips.bandY[k] + 6);
    expect([h?.block, h?.system]).toEqual([11, '/data/s1']);
    expect(h!.position! >= 50 && h!.position! < 55).toBe(true);
    // a block in a run has no label of its own: the run's
    expect(p.shown[pt.strip]).toBe('');
  });

  it('across a block gap: a hit in the gap names no block’s data, the next block’s first position is its own', async () => {
    const s = await settle(mixedSpec());
    const f = await s.frame({ size: Size });
    const p = placementOf(f.layout);
    const a0 = p.strips.axes[0];
    const a1 = p.strips.axes[1];
    expect(a1.x0 - (a0.x0 + a0.width)).toBeGreaterThan(1);
    const y = p.strips.bandY[0] + 4;
    const gap = new PlacePoint();
    p.at([(a0.x0 + a0.width + a1.x0) / 2, y], gap);
    expect(gap.kind).not.toBe('data');
    expect(hitView(f.layout, a1.x0 + 0.2, y)?.position).toBe(399);
  });
});

describe('gaps at scale and runs (EC-102)', () => {
  const row = (layout: Layout) => placementOf(layout).row;

  it('the gaps in view share at most a quarter of the width, 12 px each at most; none inside a run', async () => {
    const s = await settle(genomeSpec(toyGenome(25, 3000)));
    const f = await s.frame({ size: { width: 1200, height: 300 } });
    const r = row(f.layout);
    let gaps = 0;
    for (let b = 1; b < r.n; b++) {
      const g = r.x[b] - (r.x[b - 1] + r.w[b - 1]);
      expect(g).toBeLessThanOrEqual(12 + 1);
      if (r.run[b] >= 0 && r.run[b] === r.run[b - 1]) expect(Math.abs(g)).toBeLessThan(1);
      gaps += Math.max(0, g);
    }
    expect(gaps).toBeLessThanOrEqual(0.25 * r.width + 2);
    const p = placementOf(f.layout);
    expect(p.runN).toBeGreaterThan(0);
    expect(p.runText.slice(0, p.runN).some((t) => /contigs$/.test(t))).toBe(true);
  });

  it('inside a run a block is read and drawn while a device pixel wide, else not (the joint rule)', async () => {
    const s = await settle(genomeSpec(toyGenome(3, 300)));
    const at = async (dpr: number) => {
      const f = await s.frame({ size: { width: 1000, height: 300 }, dpr });
      const p = placementOf(f.layout);
      let inRuns = 0;
      for (let i = 0; i < p.strips.n; i++) if (p.row.run[p.strips.block[i]] >= 0) inRuns++;
      return { inRuns, runs: p.runN };
    };
    const one = await at(1);
    const four = await at(4);
    expect(one.runs).toBeGreaterThan(0);
    expect(four.inRuns).toBeGreaterThan(one.inRuns);
  });

  it('runs part again as the zoom grows, and are found again on a zoom, never on a pan', async () => {
    const s = await settle(genomeSpec(toyGenome(3, 300)));
    const f = await s.frame({ size: { width: 1000, height: 300 } });
    const total = f.layout.sections[0].axis!.total;
    const runs = (l: Layout) => placementOf(l).runN;
    expect(runs(f.layout)).toBeGreaterThan(0);
    const zoomed = await s.frame({ size: { width: 1000, height: 300 }, viewport: { start: total - 2e6, end: total } });
    expect(
      placementOf(zoomed.layout)
        .row.run.slice(0, placementOf(zoomed.layout).row.n)
        .every((r) => r < 0),
    ).toBe(true);
    // a pan keeps the zoom, so the runs: the same blocks in each
    const a = await s.frame({ size: { width: 1000, height: 300 }, viewport: { start: total * 0.5, end: total * 0.9 } });
    const runsA = Array.from(placementOf(a.layout).row.run.slice(0, placementOf(a.layout).row.n));
    const b = await s.frame({
      size: { width: 1000, height: 300 },
      viewport: { start: total * 0.55, end: total * 0.95 },
    });
    expect(Array.from(placementOf(b.layout).row.run.slice(0, placementOf(b.layout).row.n))).toEqual(runsA);
  });

  it('a ring of 3,000 blocks keeps three quarters of its circle for positions, with runs or without (EC-102)', async () => {
    const spec = genomeSpec(toyGenome(24, 2976), { layout: 'rings' });
    for (const runPx of [2, 0]) {
      const s = await settle(spec, [], { blocks: { runPx } });
      const f = await s.frame({ size: { width: 800, height: 800 } });
      const p = f.layout.sections[0].placement as RingPlacement;
      let positions = 0;
      for (let b = 0; b < 3000; b++) positions += p.bw[b];
      expect(positions / p.circumference).toBeGreaterThanOrEqual(0.75 - 1e-9);
      // the old rule (1° a block) left no ring past 360 blocks: without runs the gaps are 0.03° each
      if (runPx === 0) {
        expect((p.gapAngle * 180) / Math.PI).toBeCloseTo(0.03, 3);
        expect(p.runN).toBe(0);
        expect(p.strips.n).toBe(3000);
      } else expect([p.runN > 0, p.strips.n >= 24]).toEqual([true, true]);
    }
  });

  it('a ring’s title is its assembly’s name and length in bases', async () => {
    const s = await settle(genomeSpec(toyGenome(5, 0), { layout: 'rings', name: 'GRCh38' }));
    const f = await s.frame({ size: { width: 600, height: 600 } });
    const p = f.layout.sections[0].placement as RingPlacement;
    expect(p.title).toBe('GRCh38');
    expect(p.subtitle).toMatch(/^\d+(\.\d+)? (Gb|Mb)$/);
  });
});

describe('two assemblies on one section (EC-210)', () => {
  it('a row: a group gap wider than a block gap, a bracket a group with its name', async () => {
    const s = await settle(twoGenomesSpec());
    const f = await s.frame({ size: { width: 1200, height: 300 } });
    const p = placementOf(f.layout);
    const r = p.row;
    const gapAt = (b: number) => r.x[b] - (r.x[b - 1] + r.w[b - 1]);
    // block 8 is grape's first
    expect(gapAt(8)).toBeGreaterThan(gapAt(4) * 2);
    expect(p.groupText.slice(0, p.groupN)).toEqual(['peach', 'grape']);
    expect(p.footRow).toBeGreaterThan(0);
  });

  it('a ring: a 4° group gap each side, the groups’ names placed outside', async () => {
    const s = await settle(twoGenomesSpec('rings'));
    const f = await s.frame({ size: { width: 800, height: 800 } });
    const p = f.layout.sections[0].placement as RingPlacement;
    const gap = (p.bx[8] - (p.bx[7] + p.bw[7])) / p.r;
    expect((gap * 180) / Math.PI).toBeCloseTo(4, 1);
    expect(p.groupN).toBe(2);
    expect(p.title).toBe('peach · grape');
  });
});

describe('a hue per block (EC-218)', () => {
  it('`block` takes the section’s blocks in order as its domain, never the data’s first occurrence', async () => {
    const b = createSeqSpec();
    const s0 = b.data.sequence('s0', { alphabet: 'dna', letters: dna(200, 1) });
    const s1 = b.data.sequence('s1', { alphabet: 'dna', letters: dna(200, 2) });
    const f0 = b.data.ranges('f0', { target: s0, rows: { start: [10], end: [50], id: ['a'] } });
    const f1 = b.data.ranges('f1', { target: s1, rows: { start: [10], end: [50], id: ['b'] } });
    b.views.view('v', (v) =>
      v.section('g', {}, (s) => {
        // s1's block first: its feature takes the first category
        s.block('x', { of: s1 });
        s.block('y', { of: s0 });
        const t = s.track('f');
        t.layer('f0', { type: 'ranges', data: f0, color: theme.field('block') });
        t.layer('f1', { type: 'ranges', data: f1, color: theme.field('block') });
      }),
    );
    const st = await settle(b.build());
    const f = await st.frame({ size: Size });
    const sec = f.layout.sections[0];
    const fills: Record<string, number> = {};
    for (const lays of sec.lays)
      for (const l of lays ?? []) {
        const lay = l.lay as RangesLay;
        for (let i = 0; i < lay.items; i++) fills[String(lay.keys[i])] = lay.fills[i];
      }
    expect(fills.b).toBe(LightTokens.colors['seq-cat-0']);
    expect(fills.a).toBe(LightTokens.colors['seq-cat-1']);
  });
});

describe('the export', () => {
  it('writes each block’s system, direction and scale, and an assembly’s address', async () => {
    const { svg } = await renderSvg(mixedSpec(), { size: Size });
    expect(svg).toContain('data-block="/views/v/g/b1" data-system="/data/s0" data-direction="reversed" data-scale="1"');
    expect(svg).toContain('data-block="/views/v/g/b2" data-system="/data/s1" data-direction="genomic" data-scale="3"');
    const two = await renderSvg(twoGenomesSpec(), { size: Size });
    expect(two.svg).toContain('data-assembly="/data/grape"');
  });

  it('takes a viewport by `{ system, position }` ends, as the same window in units', async () => {
    const at = '/views/v/g';
    const byPlace = await renderSvg(mixedSpec(), {
      size: Size,
      idPrefix: 'v',
      viewport: {
        positions: {
          [at]: [
            { system: '/data/s0', position: 300 },
            { system: '/data/s1', position: 50 },
          ],
        },
      },
    });
    // 300 on the reversed block [150, 400) is unit 150 + 99 = 249; 50 on the ×3 block is 400 + 150 = 550
    const byUnits = await renderSvg(mixedSpec(), {
      size: Size,
      idPrefix: 'v',
      viewports: { [at]: { start: 250, end: 550 } },
    });
    expect(byPlace.svg).toBe(byUnits.svg);
  });
});

describe('equal block lists follow', () => {
  it('a ring of the genome and its linear detail: the detail’s window is the ring’s, no issue', async () => {
    const b = createSeqSpec();
    const { assembly, genome } = await import('@molstar/seqx-core/lang');
    b.data.assembly('g', assembly('g', toyGenome(6, 0)));
    b.views.view('v', (v) => {
      genome(v, 'ring', '/data/g', { layout: 'rings' });
      genome(v, 'detail', '/data/g');
    });
    const spec = b.build();
    // the detail follows the ring through the section's `follow`
    const views = spec.children!.find((c) => c.kind === 'views')!.children![0];
    const detail = views.children!.find((c) => c.id === 'detail')!;
    (detail as { params?: Node['params'] }).params = { ...detail.params, follow: '/views/v/ring' };
    const st = await settle(spec);
    const f = await st.frame({
      size: { width: 800, height: 900 },
      viewports: { '/views/v/ring': { start: 1e6, end: 5e7 } },
    });
    expect(f.layout.issues.map((i) => i.code)).not.toContain('seq.render.follow.blocks');
    const d = f.layout.sections.find((x) => x.section.address === '/views/v/detail')!;
    expect([d.axis!.u0, d.axis!.u1]).toEqual([1e6, 5e7]);
  });
});

void bundledTextMetrics;
