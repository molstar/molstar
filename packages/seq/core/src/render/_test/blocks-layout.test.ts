/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Blocks: the resolved view's blocks (copies merged `×N`, or each; a layer placed on each block its data reaches),
// the tracks layout's row (one axis over the blocks end to end, gaps between) and stack (a block a row, aligned at
// their first position), placement as the inverse of hit-testing in both, hits naming their block, the rings layout's
// arcs, a follow across blocks of several systems refused, drawing against the data.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { type BlockPlacement, type Frame, hitView, PlacePoint, type RingPlacement } from '../index.js';
import { settle } from './helpers.js';
import { blocksSpec } from './rings-fixtures.js';

const Size = { width: 800, height: 600 };

describe('blocks in the resolved view', () => {
  it('copies merged into one block labeled ×N (the default), or a block each', () => {
    const merged = resolveViews(resolveSpec(blocksSpec()))[0].sections[0];
    expect(merged.blocks!.map((b) => b.label)).toEqual(['CDK2 ×2', 'cyclin A', 'ATP']);
    expect(merged.blocks![0].members).toEqual(['/views/v/complex/a', '/views/v/complex/c']);
    expect(merged.blocks!.map((b) => b.offset)).toEqual([0, 120, 200]);
    expect(merged.axis).toMatchObject({ pieces: 'blocks', ranges: [0, 201], cyclic: false });
    const each = resolveViews(resolveSpec(blocksSpec({ copies: 'each' })))[0].sections[0];
    expect(each.blocks!.map((b) => b.label)).toEqual(['CDK2', 'cyclin A', 'CDK2', 'ATP']);
    expect(each.axis.ranges).toEqual([0, 321]);
  });

  it('a path a block: a layer placed on each block its data reaches, no issue for those it does not', () => {
    const v = resolveViews(resolveSpec(blocksSpec()))[0];
    const [cdk2, cyclin, atp] = v.sections[0].blocks!;
    const placed = (b: typeof cdk2, layer: number) => b.layers[0][layer].placed;
    expect([placed(cdk2, 0), placed(cdk2, 1)]).toEqual([true, false]);
    expect([placed(cyclin, 0), placed(cyclin, 1)]).toEqual([false, true]);
    expect([placed(atp, 0), placed(atp, 1)]).toEqual([false, false]);
    expect(v.issues.filter((i) => i.code === 'seq.view.layer.unplaced')).toEqual([]);
  });

  it('a stack’s axis is its longest block', () => {
    const s = resolveViews(resolveSpec(blocksSpec({ arrange: 'stack' })))[0].sections[0];
    expect(s.axis.ranges).toEqual([0, 120]);
    expect(s.blocks!.map((b) => b.offset)).toEqual([0, 0, 0]);
  });
});

/** The feature keys under each position of each block, by hits across its features band. */
function keysAt(f: Frame, p: BlockPlacement, strip: number, unit: number): string[] {
  const s = p.strips;
  const a = s.axes[strip];
  const x = a.xOfUnit(unit - s.first[strip] + 0.5);
  const keys = new Set<string>();
  for (let y = s.bandY[strip * s.bands]; y < s.bandY[strip * s.bands] + s.bandH[strip * s.bands]; y += 2) {
    const h = hitView(f.layout, x, y);
    if (h?.element) keys.add(String(h.element.key));
  }
  return [...keys];
}

describe('blocks in the tracks layout', () => {
  it('a row: the blocks end to end with gaps, placement the inverse of hit-testing, hits naming their block', async () => {
    const s = await settle(blocksSpec());
    const f = await s.frame({ size: Size });
    const sec = f.layout.sections[0];
    const p = sec.placement as BlockPlacement;
    expect(p.strips.n).toBe(3);
    // gaps between the blocks, in order
    for (let i = 1; i < 3; i++)
      expect(p.strips.axes[i].x0).toBeGreaterThan(p.strips.axes[i - 1].x0 + p.strips.axes[i - 1].width);
    const pt = new PlacePoint();
    let last = -Infinity;
    for (let u = 0; u < 201; u++) {
      expect(p.place(u + 0.5, pt)).toBe(true);
      expect(p.at(pt.x, pt.y, pt)).toBe(true);
      expect(pt.kind).toBe('data');
      expect(Math.floor(pt.offset)).toBe(u);
      const o = p.offsetAt(pt.x, pt.y);
      expect(o).toBeGreaterThan(last);
      last = o;
    }
    // the features drawn where the data has them, each on its own block
    expect(keysAt(f, p, 0, 45)).toEqual(['psta']);
    expect(keysAt(f, p, 0, 10)).toEqual([]);
    expect(keysAt(f, p, 1, 120 + 30)).toEqual(['box']);
    const h = hitView(f.layout, p.strips.axes[1].xOfUnit(30.5), p.strips.bandY[p.strips.bands] + 2);
    expect(h).toMatchObject({ block: 1, system: '/data/cyclin' });
    // the block labels above their tracks
    const texts = f.lists.passes.labels.texts;
    const strings = Array.from({ length: texts.n }, (_, i) => f.lists.string(texts.str[i]));
    expect(strings).toEqual(expect.arrayContaining(['CDK2 ×2', 'cyclin A', 'ATP']));
  });

  it('a row zoomed into a block shows it alone; its lays kept across a pan', async () => {
    const s = await settle(blocksSpec());
    const f1 = await s.frame({ size: Size, viewport: { start: 125, end: 185 } });
    const p1 = f1.layout.sections[0].placement as BlockPlacement;
    expect(p1.strips.n).toBe(1);
    expect(p1.strips.axes[0].system).toBe('/data/cyclin');
    expect([p1.strips.axes[0].u0, p1.strips.axes[0].u1]).toEqual([5, 65]);
    const lay = f1.layout.sections[0].lays[0][0].lay;
    const f2 = await s.frame({ size: Size, viewport: { start: 130, end: 190 } });
    expect(f2.layout.sections[0].lays[0][0].lay).toBe(lay);
  });

  it('a row of many blocks: the gaps never take more than a quarter of the width (40 and 100 blocks)', async () => {
    for (const n of [40, 100]) {
      const b = createSeqSpec({ title: 'many blocks' });
      const systems = Array.from({ length: n }, (_, k) =>
        b.data.sequence(`s${k}`, { alphabet: 'protein', letters: 'ACDEFGHIKL' }),
      );
      b.views.view('v', (v) => {
        const sec = v.section('many', {});
        systems.forEach((sys, k) => sec.block(`b${k}`, { of: sys, label: `B${k}` }));
        sec.track('letters').layer('letters', { type: 'letters' });
      });
      const s = await settle(b.build());
      const f = await s.frame({ size: Size });
      const sec = f.layout.sections[0];
      const p = sec.placement as BlockPlacement;
      expect(p.strips.n).toBe(n);
      const axis = sec.axis!;
      let gaps = 0;
      for (let i = 1; i < n; i++) gaps += p.strips.axes[i].x0 - (p.strips.axes[i - 1].x0 + p.strips.axes[i - 1].width);
      expect(gaps).toBeGreaterThan(0);
      expect(gaps).toBeLessThanOrEqual(0.25 * axis.width + 1e-6);
      const last = p.strips.axes[n - 1];
      expect(last.x0 + last.width).toBeCloseTo(axis.x0 + axis.width, 6);
    }
  });

  it('a stack: a block a row of the tracks, aligned at their first position', async () => {
    const s = await settle(blocksSpec({ arrange: 'stack' }));
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as BlockPlacement;
    expect(p.stacked).toBe(true);
    expect(p.strips.n).toBe(3);
    // aligned: each block's first position at the same x, the rows one under another
    expect(p.strips.axes.map((a) => a.x0)).toEqual([p.strips.axes[0].x0, p.strips.axes[0].x0, p.strips.axes[0].x0]);
    expect(p.strips.top[1]).toBeGreaterThan(p.strips.top[0] + p.strips.height[0] - 1);
    expect(p.strips.axes[1].scale).toBeCloseTo(p.strips.axes[0].scale);
    const pt = new PlacePoint();
    p.place(30.5, pt);
    expect(pt.strip).toBe(0);
    // a point on the second row is the second block's
    const y = p.strips.bandY[p.strips.bands] + 2;
    expect(p.at(p.strips.axes[1].xOfUnit(30.5), y, pt) && pt.strip).toBe(1);
    expect(keysAt(f, p, 1, 30)).toEqual(['box']);
    // each row draws its own letters (a layer of the same key on every block: its lays kept a block each)
    for (let i = 0; i < p.strips.n; i++) {
      const k = i * p.strips.bands + 1;
      const h = hitView(f.layout, p.strips.axes[i].xOfUnit(0.5), p.strips.bandY[k] + 2);
      expect(h).toMatchObject({ block: i, position: 0, system: p.strips.axes[i].system });
    }
  });
});

describe('blocks on a ring', () => {
  it('a block an arc with a gap between, each on its system', async () => {
    const s = await settle(blocksSpec({ layout: 'rings' }));
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as RingPlacement;
    expect(p.strips.n).toBe(3);
    for (let i = 1; i < 3; i++)
      expect(p.strips.axes[i].x0).toBeGreaterThan(p.strips.axes[i - 1].x0 + p.strips.axes[i - 1].width);
    expect(p.strips.axes.map((a) => a.system)).toEqual(['/data/cdk2', '/data/cyclin', '/data/atp']);
    const pt = new PlacePoint();
    for (let u = 0; u < 201; u += 3) {
      p.place(u + 0.5, pt);
      expect(p.at(pt.x, pt.y, pt) && Math.floor(pt.offset)).toBe(u);
    }
    // the block labels lead the outside labels
    const placed = Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]);
    expect(placed).toEqual(expect.arrayContaining(['CDK2 ×2', 'cyclin A', 'ATP', 'PSTAIRE', 'cyclin box']));
  });

  it('a section of blocks on several systems neither follows nor is followed', async () => {
    const { createSeqSpec } = await import('@molstar/seqx-core/lang');
    const b = createSeqSpec({ title: 'follow blocks' });
    const a = b.data.sequence('a', { alphabet: 'protein', letters: 'ACDEFGHIKL' });
    const c = b.data.sequence('c', { alphabet: 'protein', letters: 'MNPQRSTVWY' });
    b.views.view('v', (v) => {
      const s = v.section('blocks', {});
      s.block('a', { of: a });
      s.block('c', { of: c });
      s.track('l').layer('l', { type: 'letters' });
      v.section('detail', { of: a, follow: s }).track('l').layer('l', { type: 'letters' });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: Size });
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.follow.blocks');
    expect(f.layout.sections[1].axis!.total).toBe(10);
  });
});
