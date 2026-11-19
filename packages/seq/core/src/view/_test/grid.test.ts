/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Pairs and grids resolved (_spec/README.md "Resolve", Pairs and grids; "Resolved view", Grids): held whole on their two
// systems, a grid of the wrong size dropped with the count it expected, a grid section's y and what it follows.

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../../lang/builder.js';
import { resolveSpec } from '../../resolve/model.js';
import { resolveView } from '../view.js';

function spec(followY?: string, ofY = true) {
  const b = createSeqSpec();
  const x = b.data.sequence('x', { alphabet: 'protein', letters: 'MEEPQ' });
  const y = b.data.sequence('y', { alphabet: 'protein', letters: 'MKL' });
  b.data.pairs('p', {
    a: x,
    b: y,
    rows: {
      a_start: [0, 4, 3],
      a_end: [1, 5, 3],
      b_start: [2, 0, 1],
      b_end: [3, 1, 2],
      type: ['contact', null, 'contact'],
      a_strand: ['+', '-', '+'],
    },
  });
  b.data.grid('g', { x, y, values: Float32Array.from({ length: 15 }, (_, i) => i) });
  b.views.view('v', (v) => {
    v.section('x', { of: x }).track('t').layer('l', { type: 'letters' });
    v.section('y', { of: y }).track('t').layer('l', { type: 'letters' });
    v.section('plot', {
      layout: 'grid',
      of: x,
      ...(ofY ? { of_y: y } : {}),
      ...(followY ? { follow_y: followY } : {}),
    })
      .track('t')
      .layer('g', { type: 'grid', data: '/data/g' });
  });
  return b.build();
}

describe('pairs and grids held whole', () => {
  it('reads the pairs’ ends on both systems, their kinds and keys; a grid x-major', () => {
    const m = resolveSpec(spec('none'));
    const p = m.annotations.get('/data/p')!;
    expect(p.shape).toBe('pairs');
    expect(p.target).toBe('/data/x');
    expect(p.ySystem).toBe('/data/y');
    expect(p.access).toBeUndefined();
    // the third row's first end is empty: left out (the checks say so)
    expect(p.pairs!.n).toBe(2);
    expect([...p.pairs!.aStart]).toEqual([0, 4]);
    expect([...p.pairs!.bEnd]).toEqual([3, 1]);
    expect(p.pairs!.kinds).toEqual(['contact']);
    expect([...p.pairs!.kind]).toEqual([0, -1]);
    expect([...p.pairs!.aStrand]).toEqual([1, -1]);
    expect(p.pairs!.keys).toEqual([0, 1]);
    const g = m.annotations.get('/data/g')!;
    expect(g.shape).toBe('grid');
    expect(g.grid).toMatchObject({ nx: 5, ny: 3 });
    expect(g.grid!.values[1 * 3 + 2]).toBe(5);
  });

  it('drops a grid of the wrong size, saying the count it expected', () => {
    const b = createSeqSpec();
    const x = b.data.sequence('x', { alphabet: 'protein', letters: 'MEE' });
    b.data.grid('g', { x, values: [1, 2, 3] });
    const m = resolveSpec(b.build());
    expect(m.annotations.get('/data/g')!.grid).toBeUndefined();
    const issue = m.issues.find((i) => i.code === 'seq.resolve.grid.size')!;
    expect(issue.message).toContain('not 9');
  });
});

describe('a grid section’s y', () => {
  it('follows x by default: the same system, its window', () => {
    const s = resolveView(resolveSpec(spec(undefined, false)), '/views/v')!.sections[2];
    expect(s.layout).toBe('grid');
    expect(s.y).toMatchObject({ followsX: true, axis: { system: '/data/x', length: 5, ranges: [0, 5] } });
    expect(s.y!.follow).toBeUndefined();
  });

  it('follows a section through the path from its system, or keeps a window of its own', () => {
    const m = resolveSpec(spec('/views/v/y'));
    const s = resolveView(m, '/views/v')!.sections[2];
    expect(s.y!.axis.system).toBe('/data/y');
    expect(s.y!.follow?.leader).toBe('/views/v/y');
    expect(s.y!.followsX).toBe(false);
    const own = resolveView(resolveSpec(spec('none')), '/views/v')!.sections[2];
    expect(own.y).toMatchObject({ followsX: false, axis: { system: '/data/y' } });
    expect(own.y!.follow).toBeUndefined();
  });

  it('says a y that cannot follow x (no path between the systems), and keeps its own window', () => {
    const v = resolveView(resolveSpec(spec()), '/views/v')!;
    expect(v.sections[2].y!.followsX).toBe(false);
    expect(v.issues.map((i) => i.code)).toContain('seq.view.follow.unplaced');
  });
});
