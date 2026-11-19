/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved axis (_spec/README.md "Resolved view", Axis): a section's direction as its writer spells it, its pieces
// with their breaks, and offsets and positions through them both ways; property tests over random sections.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../../lang/builder.js';
import { resolveSpec } from '../../resolve/model.js';
import {
  type AxisShape,
  axisBreaks,
  axisKey,
  axisTotal,
  breakHolding,
  cellOffset,
  elementAt,
  monotonic,
  offsetOf,
  positionAt,
} from '../axis.js';
import { resolveView } from '../view.js';

// a protein of 4 residues on a 30-nucleotide gene, exons [2, 7) and [12, 19) (12 nucleotides), on a strand
function spec(
  gene: Record<string, unknown> = {},
  strand: '+' | '-' = '-',
  options: { layer?: boolean; second?: boolean } = {},
) {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGTAC'.repeat(5) });
  const other = b.data.sequence('other', { alphabet: 'dna', letters: 'ACGT' });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [2, 12], end: [7, 19] }, strand });
  if (options.second) b.data.map('m2', { from: p, to: g, method: 'codon', exons: { start: [0], end: [12] }, strand });
  b.data.ranges('exons', { target: g, rows: { start: [2, 12], end: [7, 19] } });
  b.data.ranges('far', { target: other, rows: { start: [0], end: [2] } });
  b.views.view('v', (v) => {
    const protein = v.section('protein', { of: p });
    protein.track('t').layer('l', { type: 'letters' });
    const s = v.section('gene', { of: g, ...(gene as object) });
    s.track('t').layer('l', { type: 'letters' });
    if (options.layer) s.track('d').layer('d', { type: 'ranges', data: '/data/exons' });
    v.section('lone', { of: other });
  });
  return b.build();
}

function gene(
  params: Record<string, unknown>,
  strand: '+' | '-' = '-',
  options: { layer?: boolean; second?: boolean } = {},
) {
  const view = resolveView(resolveSpec(spec(params, strand, options)), '/views/v')!;
  return { axis: view.sections[1].axis, issues: view.issues.map((i) => i.code) };
}

describe('a section’s direction', () => {
  it('resolves each value as written, defaults by the section’s own follow param', () => {
    expect(gene({}).axis.dir).toBe(1);
    expect(gene({ direction: 'reversed' }).axis.dir).toBe(-1);
    expect(gene({ follow: '/views/v/protein' }).axis.dir).toBe(-1);
    expect(gene({ follow: '/views/v/protein' }, '+').axis.dir).toBe(1);
    expect(gene({ follow: '/views/v/protein', direction: 'genomic' }).axis.dir).toBe(1);
    expect(gene({ direction: { along: '/data/p' } }).axis.dir).toBe(-1);
    expect(gene({ direction: { along: '/data/m' } }).axis.dir).toBe(-1);
    expect(gene({ direction: { along: '/data/p' } }, '+').axis.dir).toBe(1);
    for (const g of [gene({}), gene({ direction: 'reversed' }), gene({ direction: { along: '/data/p' } })])
      expect(g.issues).toEqual([]);
  });

  it('reports `along` a system no path reaches and `follow` with no leader, and is genomic', () => {
    expect(gene({ direction: { along: '/data/other' } })).toMatchObject({
      axis: { dir: 1 },
      issues: ['seq.view.direction.unresolved'],
    });
    expect(gene({ direction: 'follow' })).toMatchObject({
      axis: { dir: 1 },
      issues: ['seq.view.direction.unresolved'],
    });
  });

  it('is unchanged when a layer is added: never read off what the section draws', () => {
    for (const params of [{}, { follow: '/views/v/protein' }, { direction: { along: '/data/p' } }])
      expect(gene(params, '-', { layer: true }).axis).toEqual(gene(params).axis);
  });

  it('follows a follower’s leader: a chain multiplies the directions', () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
    const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGTAC'.repeat(5) });
    b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [2, 12], end: [7, 19] }, strand: '-' });
    b.views.view('v', (v) => {
      const protein = v.section('protein', { of: p, direction: 'reversed' });
      const gs = v.section('gene', { of: g, follow: protein });
      v.section('back', { of: p, follow: gs });
    });
    const view = resolveView(resolveSpec(b.build()), '/views/v')!;
    expect(view.sections.map((s) => s.axis.dir)).toEqual([-1, 1, -1]);
  });
});

describe('a section’s pieces', () => {
  it('whole by default, coding through the follow path or the path read along, refs and their union', () => {
    expect(gene({}).axis).toMatchObject({ pieces: 'whole', ranges: [0, 30] });
    // coding on a reversed axis: the higher exon first
    expect(gene({ follow: '/views/v/protein', pieces: 'coding' }).axis).toMatchObject({
      pieces: 'coding',
      dir: -1,
      ranges: [12, 19, 2, 7],
    });
    expect(gene({ direction: { along: '/data/p' }, pieces: 'coding' }).axis.ranges).toEqual([12, 19, 2, 7]);
    expect(gene({ pieces: '/data/m' }).axis.ranges).toEqual([2, 7, 12, 19]);
    expect(gene({ pieces: ['/data/exons', '/data/m'], range: [4, 30] }).axis.ranges).toEqual([4, 7, 12, 19]);
    expect(gene({ pieces: '/data/exons', pad: 3 }).axis.ranges).toEqual([0, 22]);
    expect(gene({ pieces: '/data/exons', pad: 2 }).axis.ranges).toEqual([0, 9, 10, 21]);
  });

  it('reports pieces that reach nothing, and shows whole', () => {
    expect(gene({ pieces: '/data/far' })).toMatchObject({
      axis: { pieces: 'whole', ranges: [0, 30] },
      issues: ['seq.view.section.pieces'],
    });
  });

  it('`coding` on a section that neither follows nor reads along: the one codon map onto its system', () => {
    expect(gene({ pieces: 'coding' })).toMatchObject({
      axis: { pieces: 'coding', ranges: [2, 7, 12, 19] },
      issues: [],
    });
    expect(gene({ pieces: 'coding', direction: 'reversed' }).axis.ranges).toEqual([12, 19, 2, 7]);
    // two codon maps onto it and none named: ambiguous, whole
    expect(gene({ pieces: 'coding' }, '-', { second: true })).toMatchObject({
      axis: { pieces: 'whole', ranges: [0, 30] },
      issues: ['seq.view.section.pieces'],
    });
    // the one named still decides
    expect(gene({ pieces: '/data/m' }, '-', { second: true }).axis.ranges).toEqual([2, 7, 12, 19]);
  });

  it('pieces all outside the range: reported, the range shown (never an axis waiting for pieces)', () => {
    expect(gene({ pieces: '/data/exons', range: [20, 30] })).toMatchObject({
      axis: { pieces: 'whole', ranges: [20, 30] },
      issues: ['seq.view.section.pieces'],
    });
    expect(gene({ pieces: '/data/exons', range: [4, 14] }).axis.ranges).toEqual([4, 7, 12, 14]);
  });

  it('has its breaks with their hidden ranges, in display order', () => {
    const fwd = gene({ pieces: '/data/m' }).axis;
    expect(axisBreaks(fwd)).toEqual([5, 7, 12]);
    const rev = gene({ follow: '/views/v/protein', pieces: 'coding' }).axis;
    expect(axisBreaks(rev)).toEqual([7, 7, 12]);
    expect(breakHolding(rev, 9)).toBe(0);
    expect(breakHolding(rev, 12)).toBe(-1);
    // contiguous across a cyclic origin: no break, either way; a gap across it hides both sides
    for (const [ranges, dir] of [
      [[8, 10, 0, 2], 1],
      [[0, 2, 8, 10], -1],
    ] as const)
      expect(axisBreaks({ ranges, dir, length: 10, cyclic: true })).toEqual([]);
    expect(axisBreaks({ ranges: [1, 3, 5, 9], dir: 1, length: 10, cyclic: true })).toEqual([2, 3, 5]);
    const wrap = { ranges: [5, 9, 1, 3], dir: 1 as const, length: 10, cyclic: true };
    expect(axisBreaks(wrap)).toEqual([4, 9, 1]);
    expect(breakHolding(wrap, 0)).toBe(0);
    expect(elementAt(wrap, 4.5)).toBe(1);
  });
});

describe('offsets and positions through the axis', () => {
  const rev = gene({ follow: '/views/v/protein', pieces: 'coding' }).axis;

  it('a reversed piece shows position `end − 1 − k` at its offset `k` (EC-32)', () => {
    expect(elementAt(rev, 0)).toBe(18);
    expect(elementAt(rev, 6.5)).toBe(12);
    expect(elementAt(rev, 7)).toBe(6);
    expect(cellOffset(rev, 18)).toBe(0);
    expect(cellOffset(rev, 2)).toBe(11);
    expect(positionAt(rev, 0.25)).toBe(18.75);
    expect(offsetOf(rev, 18.75)).toBe(0.25);
  });

  it('an offset on a break is the next piece’s first edge as a start, the previous one’s last as an end', () => {
    expect(positionAt(rev, 7, 'end')).toBe(12);
    expect(positionAt(rev, 7, 'start')).toBe(7);
    const fwd = gene({ pieces: '/data/m' }).axis;
    expect(positionAt(fwd, 5, 'end')).toBe(7);
    expect(positionAt(fwd, 5, 'start')).toBe(12);
  });

  it('round trips over random sections: direction, pieces, a cyclic origin', () => {
    const rng = new PCG(7);
    for (let trial = 0; trial < 200; trial++) {
      const length = 5 + Math.floor(rng.float() * 60);
      const cuts = new Set<number>();
      const n = 2 * (1 + Math.floor(rng.float() * 5));
      while (cuts.size < Math.min(n, length)) cuts.add(Math.floor(rng.float() * (length + 1)));
      const sorted = [...cuts].sort((a, b) => a - b);
      if (sorted.length % 2) sorted.pop();
      if (!sorted.length) continue;
      const dir = Math.floor(rng.float() * 2) ? 1 : -1;
      const pieces: number[] = [];
      for (let i = 0; i + 1 < sorted.length; i += 2)
        if (sorted[i + 1] > sorted[i]) pieces.push(sorted[i], sorted[i + 1]);
      if (!pieces.length) continue;
      const display = dir === 1 ? pieces : reverse(pieces);
      const axis: AxisShape = { ranges: display, dir, length, cyclic: false };
      const total = axisTotal(axis)!;
      // every shown position placed once
      const seen = new Map<number, number>();
      for (let u = 0; u < total; u++) {
        const p = elementAt(axis, u + 0.5)!;
        expect(seen.has(p)).toBe(false);
        seen.set(p, u);
        expect(cellOffset(axis, p)).toBe(u);
      }
      let shown = 0;
      for (let i = 0; i + 1 < pieces.length; i += 2) shown += pieces[i + 1] - pieces[i];
      expect(seen.size).toBe(shown);
      // fractional offsets and positions both ways, inside a piece
      for (let i = 0; i < 10; i++) {
        const u = rng.float() * total;
        const c = positionAt(axis, u, 'start')!;
        expect(offsetOf(axis, c)).toBeCloseTo(u, 9);
      }
      // the breaks hide exactly what no piece shows between the first and last shown
      for (let i = 0; i + 2 < axisBreaks(axis).length; i += 3) {
        const [o, s, e] = axisBreaks(axis).slice(i, i + 3);
        expect(positionAt(axis, o, 'end')).toBe(dir === 1 ? s : e);
        expect(positionAt(axis, o, 'start')).toBe(dir === 1 ? e : s);
        for (let p = s; p < e; p++) expect(seen.has(p)).toBe(false);
      }
    }
  });
});

function reverse(pieces: readonly number[]): number[] {
  const out: number[] = [];
  for (let i = pieces.length - 2; i >= 0; i -= 2) out.push(pieces[i], pieces[i + 1]);
  return out;
}

describe('an axis’s keys', () => {
  it('a history step’s key is short however many pieces, and tells pieces apart', () => {
    const many = Array.from({ length: 300_000 }, (_, i) => 2 * (i >> 1) + (i & 1) + (i >> 1));
    const key = axisKey({ system: '/s', ranges: many, dir: -1, length: 450_000, cyclic: false });
    expect(key.length).toBeLessThan(100);
    const moved = many.slice();
    moved[150_001] += 1;
    expect(axisKey({ system: '/s', ranges: moved, dir: -1, length: 450_000, cyclic: false })).not.toBe(key);
    expect(axisKey({ system: '/s', ranges: many.slice(), dir: -1, length: 450_000, cyclic: false })).toBe(key);
  });

  it('whether pieces run one way is kept per direction: one array asked both ways', () => {
    const r = [0, 2, 4, 6];
    expect(monotonic(r, 1)).toBe(true);
    expect(monotonic(r, -1)).toBe(false);
    const back = [4, 6, 0, 2];
    expect(monotonic(back, -1)).toBe(true);
    expect(monotonic(back, 1)).toBe(false);
  });
});
