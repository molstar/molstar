/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A plan of every channel kind, as large as asked: what the frame's time and allocation tests sample.

import type { Op } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '../../lang/index.js';
import { plan, type TransitionPlan } from '../index.js';
import { edit } from './fixture.js';

/** A spec of `n` tracks (a height and a hex color each) and a ring; a plan to a state changing all of them. */
export function manyChannels(n: number): TransitionPlan {
  const b = createSeqSpec({ title: 'many' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV'.repeat(10) });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    for (let i = 0; i < n; i++) s.track(`t${i}`, { height: 10 + i }).layer('l', { type: 'letters', color: '#336699' });
    v.section('ring', { of: c, layout: 'rings', rotation: 10 }).track('seq').layer('l', { type: 'letters' });
  });
  const from: Node = b.build();
  const ops: Op[] = [];
  const at = (i: number) => ['views', 'main', 'protein', `t${i}`];
  for (let i = 0; i < n; i++) {
    if (i % 10 === 9) ops.push({ op: 'remove', at: at(i) });
    else {
      ops.push({ op: 'set', at: at(i), key: 'height', value: 30 + i });
      ops.push({ op: 'set', at: [...at(i), 'l'], key: 'color', value: i % 3 ? '#cc3311' : 'token:cat-2' });
    }
  }
  ops.push({ op: 'set', at: ['views', 'main', 'ring'], key: 'rotation', value: 300 });
  ops.push({
    op: 'insert',
    parent: ['views', 'main', 'protein'],
    place: { last: true },
    node: { kind: 'track', id: 'new', params: { height: 12 }, children: [] },
  });
  const to = edit(from, ...ops);
  const bands = { '/views/main/protein/t0': { bandTop: 0, rowScale: 10 } };
  return plan(
    { root: from, view: { sections: { s: { start: 0, end: 50, top: 0, scale: 12, bands } } } },
    {
      root: to,
      view: {
        sections: {
          s: {
            start: 60,
            end: 80,
            top: 30,
            scale: 20,
            bands: { ...bands, '/views/main/protein/t0': { bandTop: 4, rowScale: 20 } },
          },
        },
      },
    },
  );
}
