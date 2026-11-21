/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, addSecondary, addVariants, conservation, P53 } from './data.js';

// p53 for the track headers' stories (track-header-parts.md): tracks with something to say in their header. One has a
// hidden layer (a second view of the same domains), one a filter (the one linker, PRR, thinned out of its data), and two tracks
// are hidden (the section header's "2 hidden" chip). The second view has every track of its section hidden; the third
// a short track beside tall ones.

type Builder = ReturnType<typeof createSeqSpec>;

function data(b: Builder) {
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  return {
    p53,
    domains: addDomains(b, p53),
    values: b.data.values('conservation', { target: p53, values: conservation() }),
  };
}

/** Domains with a hidden layer and a filter, conservation, and two hidden tracks (secondary structure, variants). */
export function build(): Node {
  const b = createSeqSpec({ title: 'p53, track headers' });
  const { p53, domains, values } = data(b);
  const ss = addSecondary(b, p53);
  const vs = addVariants(b, p53);
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    const d = s.track('domains', { label: 'Domains' });
    // the filter thins the layer's own data (what the "n filtered" chip counts); the second layer is hidden
    d.layer('ranges', { type: 'ranges', data: `${domains}/(except :of (item :category linker))` });
    d.layer('again', { type: 'ranges', data: domains, hidden: true });
    s.track('conservation', { label: 'Conservation' }).layer('values', { type: 'values', data: values });
    s.track('secondary', { label: 'Secondary', hidden: true }).layer('secondary', { type: 'secondary', data: ss });
    s.track('variants', { label: 'Variants', hidden: true }).layer('variants', { type: 'variants', data: vs });
  });
  return b.build();
}

/** Two sections: the first with every track hidden (its header keeps "2 hidden"), the second as it was. */
export function allHidden(): Node {
  const b = createSeqSpec({ title: 'p53, every track of a section hidden' });
  const { p53, domains, values } = data(b);
  b.views.view('main', (v) => {
    const hidden = v.section('hidden', { of: p53 });
    hidden.track('domains', { label: 'Domains', hidden: true }).layer('ranges', { type: 'ranges', data: domains });
    hidden
      .track('conservation', { label: 'Conservation', hidden: true })
      .layer('values', { type: 'values', data: values });
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}

/** Short tracks: a 14 px strip of domains and a 14 px conservation beside a tall one; their states are marks. */
export function short(): Node {
  const b = createSeqSpec({ title: 'p53, short tracks' });
  const { p53, domains, values } = data(b);
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    const d = s.track('domains', { label: 'Domains', height: 14 });
    d.layer('ranges', { type: 'ranges', data: `${domains}/(except :of (item :category linker))` });
    d.layer('again', { type: 'ranges', data: domains, hidden: true });
    s.track('conservation', { label: 'Conservation', height: 14 }).layer('values', { type: 'values', data: values });
    s.track('tall', { label: 'Conservation, tall', height: 48 }).layer('values', { type: 'values', data: values });
  });
  return b.build();
}
