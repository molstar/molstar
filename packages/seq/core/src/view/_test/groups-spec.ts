/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A protein section in three groups (the track groups' Defaults, small): letters alone, Features open (domains,
// motifs), Variants folded with an authored density summary, Signals open with a hidden track. Shared by the core's
// and the viewer's group tests.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

export const GroupsLetters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWF';

const n = GroupsLetters.length;
const wave = (k: number) => Float32Array.from({ length: n }, (_, i) => 0.5 + 0.5 * Math.sin((i + k) / 4));

/** The spec: `/views/v/p53` with `letters`, `features/{domains,motifs}`, `variants/{clinvar,gnomad,density}`, `signals/{a,b,c}`. */
export function groupsSpec(
  options: { variants?: Record<string, unknown>; signals?: Record<string, unknown> } = {},
): Node {
  const b = createSeqSpec({ title: 'groups' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: GroupsLetters });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: { start: new Int32Array([0, 10, 30]), end: new Int32Array([8, 28, 50]), label: ['TAD', 'PRD', 'DBD'] },
  });
  const motifs = b.data.ranges('motifs', {
    target: p53,
    rows: { start: new Int32Array([4, 20]), end: new Int32Array([6, 22]), label: ['m1', 'm2'] },
  });
  const clinvar = b.data.ranges('clinvar', {
    target: p53,
    rows: { start: new Int32Array([3, 17, 40]), end: new Int32Array([4, 18, 41]), label: ['E4K', 'T18A', 'Q41R'] },
  });
  const gnomad = b.data.ranges('gnomad', {
    target: p53,
    rows: { start: new Int32Array([7, 44]), end: new Int32Array([8, 45]), label: ['S8F', 'M45I'] },
  });
  const density = b.data.values('density', { target: p53, values: wave(0) });
  const signals = ['a', 'b', 'c'].map((id, k) => b.data.values(`signal-${id}`, { target: p53, values: wave(k * 3) }));
  b.views.view('v', (v) => {
    const s = v.section('p53', { of: p53 });
    s.track('letters', {}, (t) => t.layer('l', { type: 'letters' }));
    const f = s.group('features', { label: 'Features' });
    f.track('domains', { label: 'Domains' }, (t) => t.layer('l', { type: 'ranges', data: domains }));
    f.track('motifs', { label: 'Motifs' }, (t) => t.layer('l', { type: 'ranges', data: motifs }));
    const va = s.group('variants', { label: 'Variants', collapsed: true, summary: 'density', ...options.variants });
    va.track('clinvar', { label: 'ClinVar' }, (t) => t.layer('l', { type: 'variants', data: clinvar }));
    va.track('gnomad', { label: 'gnomAD' }, (t) => t.layer('l', { type: 'variants', data: gnomad }));
    va.track('density', { label: 'Density' }, (t) => t.layer('l', { type: 'values', data: density }));
    const sg = s.group('signals', { label: 'Signals', ...options.signals });
    ['a', 'b', 'c'].forEach((id, k) =>
      sg.track(id, { label: `Signal ${id}`, ...(id === 'c' ? { hidden: true } : {}) }, (t) =>
        t.layer('l', { type: 'values', data: signals[k] }),
      ),
    );
  });
  return b.build();
}
