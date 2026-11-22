/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A view's bridges (_spec/README.md "Resolved view", Bridges): the stack, a bridge's sides checked (neighbours, rows),
// its layers (pairs, a segments map as pairs, another type drawn as its chain's `pairs`), its params' defaults, and a
// lock a follow makes redundant.

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../../lang/builder.js';
import { resolveSpec } from '../../resolve/model.js';
import { unreadParams } from '../../resolve/unread.js';
import { ViewIssues } from '../issues.js';
import { resolveView } from '../view.js';

type Order = readonly ('top' | 'bridge' | 'bottom' | 'third' | 'second')[];

/** Two (or three) sections of two sequences with a bridge, in the order given; the bridge's params and layer as given. */
function spec(
  order: Order = ['top', 'bridge', 'bottom'],
  bridge: Record<string, unknown> = {},
  layer: Record<string, unknown> = {},
  sides: { top?: Record<string, unknown>; bottom?: Record<string, unknown> } = {},
) {
  const b = createSeqSpec();
  const x = b.data.sequence('x', { alphabet: 'dna', letters: 'ACGT'.repeat(25) });
  const y = b.data.sequence('y', { alphabet: 'dna', letters: 'ACGT'.repeat(20) });
  const p = b.data.pairs('p', {
    a: x,
    b: y,
    rows: { a_start: [0, 40], a_end: [10, 60], b_start: [5, 30], b_end: [15, 50], b_strand: ['+', '-'] },
  });
  b.data.map('runs', {
    from: x,
    to: y,
    method: 'segments',
    segments: { from: [0, 50], to: [0, 79], length: [20, 10], strand: ['+', '-'] },
  });
  b.data.map('shift', { from: x, to: y, method: 'offset', by: 3 });
  b.views.view('v', (v) => {
    for (const o of order) {
      if (o === 'bridge' || o === 'second') {
        const id = o === 'bridge' ? 'br' : 'br2';
        v.bridge(id, { a: '/views/v/top', b: '/views/v/bottom', ...bridge } as never, (br) =>
          br.layer('ribbons', { type: 'pairs', data: p, ...layer } as never),
        );
        continue;
      }
      const s = v.section(o, { of: o === 'bottom' ? y : x, ...(o === 'top' ? sides.top : sides.bottom) } as never);
      s.track('t').layer('l', { type: 'letters' });
    }
  });
  return b.build();
}

const view = (s = spec()) => resolveView(resolveSpec(s), '/views/v')!;
const codes = (v: ReturnType<typeof view>) => v.issues.map((i) => i.code);

describe('bridges in the resolved view', () => {
  it('a view without a bridge has empty lists: its stack its sections in order', () => {
    const b = createSeqSpec();
    const x = b.data.sequence('x', { alphabet: 'dna', letters: 'ACGT' });
    b.views.view('v', (v) => {
      v.section('one', { of: x }).track('t').layer('l', { type: 'letters' });
      v.section('two', { of: x }).track('t').layer('l', { type: 'letters' });
    });
    const v = view(b.build());
    expect(v.bridges).toEqual([]);
    expect(v.stack).toEqual([
      { kind: 'section', index: 0 },
      { kind: 'section', index: 1 },
    ]);
  });

  it('stacks sections and bridges in order; a bridge its sides, its defaults and its pairs layer placed', () => {
    const v = view();
    expect(v.stack).toEqual([
      { kind: 'section', index: 0 },
      { kind: 'bridge', index: 0 },
      { kind: 'section', index: 1 },
    ]);
    const br = v.bridges[0];
    expect(br).toMatchObject({
      address: '/views/v/br',
      a: '/views/v/top',
      b: '/views/v/bottom',
      label: 'top to bottom',
      lock: 'none',
      hidden: false,
      above: 'a',
      placed: true,
      issues: [],
    });
    expect(br.height).toBeUndefined();
    expect(br.layers).toHaveLength(1);
    expect(br.layers[0]).toMatchObject({
      address: '/views/v/br/ribbons',
      type: 'pairs',
      shape: 'pairs',
      placed: true,
      data: { address: '/data/p', shape: 'pairs' },
    });
    expect(br.layers[0].path).toBeUndefined();
    expect(codes(v)).toEqual([]);
  });

  it('reads its params as written: label, height, lock, hidden; `b` above when it is the section above', () => {
    const v = view(
      spec(['top', 'bridge', 'bottom'], {
        a: '/views/v/bottom',
        b: '/views/v/top',
        label: 'P–G',
        height: 140,
        lock: 'px',
        hidden: true,
      }),
    );
    expect(v.bridges[0]).toMatchObject({
      label: 'P–G',
      height: 140,
      lock: 'px',
      hidden: true,
      above: 'b',
      placed: true,
    });
  });

  it('refuses sides that are not its neighbours (A against C across B): not drawn, said at the bridge', () => {
    const v = view(spec(['bridge', 'top', 'bottom']));
    expect(v.bridges[0].placed).toBe(false);
    expect(v.bridges[0].issues.map((i) => i.code)).toEqual(['seq.view.bridge.neighbours']);
    expect(v.issues.find((i) => i.code === 'seq.view.bridge.neighbours')?.address).toBe('/views/v/br');
    // two bridges in a row: the second has a bridge above it, never a section
    const two = view(spec(['top', 'bridge', 'second', 'bottom']));
    expect(two.bridges.map((b) => b.placed)).toEqual([false, false]);
  });

  it('refuses a side that is not a row of tracks (a ring, blocks in a stack)', () => {
    const ring = view(spec(undefined, {}, {}, { top: { layout: 'rings' } }));
    expect(ring.bridges[0].placed).toBe(false);
    expect(codes(ring)).toContain('seq.view.bridge.layout');
    const grid = view(spec(undefined, {}, {}, { bottom: { layout: 'text' } }));
    expect(grid.bridges[0].placed).toBe(false);
  });

  it('reads a segments map as pairs; a map of another method is the shape issue', () => {
    const runs = view(spec(undefined, {}, { data: '/data/runs' })).bridges[0].layers[0];
    expect(runs).toMatchObject({ placed: true, data: { address: '/data/runs', system: '/data/x', shape: 'pairs' } });
    const shift = view(spec(undefined, {}, { data: '/data/shift' }));
    expect(shift.bridges[0].layers[0].placed).toBe(false);
    expect(codes(shift)).toContain('seq.view.layer.shape');
  });

  it('draws a layer of another type as its chain’s first pairs, else not at all (said)', () => {
    const ext = view(spec(undefined, {}, { type: 'rondo:ribbons', fallback: 'pairs' }));
    expect(ext.bridges[0].layers[0]).toMatchObject({ placed: true, shape: 'pairs' });
    expect(codes(ext)).toContain('seq.view.bridge.layer');
    const ranges = view(spec(undefined, {}, { type: 'ranges' }));
    expect(ranges.bridges[0].layers[0].placed).toBe(false);
    expect(codes(ranges)).toContain('seq.view.bridge.layer');
  });

  it('says a lock between a leader and its follower moves nothing more', () => {
    const v = view(spec(undefined, { lock: 'px' }, {}, { bottom: { of: '/data/x', follow: '/views/v/top' } }));
    expect(codes(v)).toContain('seq.view.bridge.lock');
    expect(v.bridges[0].placed).toBe(true);
  });

  it('is read: no unread-param issue for a bridge any more', () => {
    const issues: never[] = [];
    unreadParams(spec(), issues);
    expect(issues).toEqual([]);
  });

  it('declares every bridge code in the catalog', () => {
    for (const c of [
      'seq.view.bridge.neighbours',
      'seq.view.bridge.layout',
      'seq.view.bridge.layer',
      'seq.view.bridge.lock',
    ])
      expect(ViewIssues.codes).toHaveProperty([c]);
  });
});
