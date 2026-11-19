/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch, type Op } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { type StepInput, storeSteps, withSteps } from '@molstar/seqx-core/story';
import { addDomains, conservation, P53 } from './data.js';

// p53 as a story: four steps over its tracks, each a view of the protein and a small change of the spec (a track's
// height, the domains' color, the conservation track gone), stored by the story's encoder as patches over the spec.
// Played, the transitions pan or go out and back in, tween the heights and the color, and fade the track out.

export const title = 'p53, a story';

export const Section = '/views/main/protein';
const Domains = ['views', 'main', 'protein', 'domains'];
const DomainsLayer = [...Domains, 'ranges'];
const Conservation = ['views', 'main', 'protein', 'conservation'];

/** The spec the steps change: p53's letters, domains and conservation, at explicit heights (what a tween reads). */
function base(): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  const values = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains', height: 24 }).layer('ranges', {
      type: 'ranges',
      data: domains,
      color: '#4c78a8',
    });
    s.track('conservation', { label: 'Conservation', height: 32 }).layer('values', { type: 'values', data: values });
  });
  return b.build();
}

function edit(root: Node, ops: readonly Op[]): Node {
  const r = applyPatch(root, { ops }, { maxOps: Number.POSITIVE_INFINITY });
  if (!r.ok) throw new Error(`p53 story: ${r.issues.map((i) => i.message).join('; ')}`);
  return r.root;
}

const view = (start: number, end: number) => ({ sections: { [Section]: { start, end } } });

/** The steps' states, in order. */
export function steps(spec = base()): StepInput[] {
  const dnaBinding = edit(spec, [{ op: 'set', at: Conservation, key: 'height', value: 56 }]);
  const r175h = edit(dnaBinding, [{ op: 'set', at: DomainsLayer, key: 'color', value: '#e45756' }]);
  const tetramer = edit(spec, [
    { op: 'set', at: Domains, key: 'height', value: 40 },
    { op: 'remove', at: Conservation },
  ]);
  return [
    { id: 'whole', label: 'The whole protein', root: spec, view: view(0, P53.length), duration: 600 },
    {
      id: 'dna-binding',
      label: 'The DNA-binding domain',
      root: dnaBinding,
      view: view(93, 293),
      duration: 900,
    },
    { id: 'r175h', label: 'R175H, a hotspot', root: r175h, view: view(160, 190), duration: 900, easing: 'ease-out' },
    { id: 'tetramer', label: 'The tetramerization domain', root: tetramer, view: view(318, 363), duration: 1200 },
  ];
}

/** The view's spec: the tracks and their story. */
export function build(): Node {
  const spec = base();
  return withSteps(spec, storeSteps(spec, steps(spec)));
}
