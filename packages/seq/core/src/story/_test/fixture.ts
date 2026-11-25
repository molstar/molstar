/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small spec the stories' tests change: a protein with letters and domains, and a ring section.

import { applyPatch, type Op } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '../../lang/index.js';

export const Letters = ['views', 'main', 'protein', 'letters'];
export const Domains = ['views', 'main', 'protein', 'domains'];
export const Ring = ['views', 'main', 'ring'];

export function storySpec(): Node {
  const b = createSeqSpec({ title: 'story' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDD' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: { start: new Int32Array([0, 20]), end: new Int32Array([10, 30]), id: ['A', 'B'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters', { height: 20 }).layer('l', { type: 'letters', color: '#336699' });
    s.track('domains', { height: 16 }).layer('d', { type: 'ranges', data: domains });
    v.section('ring', { of: c, layout: 'rings', rotation: 0 }).track('seq').layer('l', { type: 'letters' });
  });
  return b.build();
}

/** The root with ops applied (the language's), throwing when they do not. */
export function edit(root: Node, ...ops: Op[]): Node {
  const r = applyPatch(root, { ops }, { maxOps: Number.POSITIVE_INFINITY });
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  return r.root;
}
