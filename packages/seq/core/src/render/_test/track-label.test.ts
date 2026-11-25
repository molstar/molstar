/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track's label in the label column (_spec/README.md "Render", Track headers): an empty one (`label: ''`) draws no
// text at all, in the tracks layout and a grid's; a named one its name.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import type { Frame } from '../index.js';
import { settle } from './helpers.js';

/** Every text the frame draws, in every pass. */
function texts(f: Frame): string[] {
  const out: string[] = [];
  for (const p of Object.values(f.lists.passes))
    for (let i = 0; i < p.texts.n; i++) out.push(f.lists.string(p.texts.str[i]));
  return out;
}

function spec() {
  const b = createSeqSpec({ title: 'labels' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  const g = b.data.grid('g', { x: s, y: s, values: new Float32Array(1600).fill(1) });
  b.views.view('v', (v) => {
    const t = v.section('t', { of: s });
    t.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    t.track('letters', { label: 'Letters' }).layer('letters', { type: 'letters' });
    v.section('plot', { layout: 'grid', of: s, of_y: s, follow_y: 'none' })
      .track('plot', { label: '' })
      .layer('g', { type: 'grid', data: g });
  });
  return b.build();
}

describe('a track’s label', () => {
  it('an empty label draws no text; a named one its name', async () => {
    const f = await (await settle(spec())).frame({ size: { width: 800, height: 900 } });
    const t = texts(f);
    expect(t).toContain('Letters');
    expect(t).not.toContain('');
  });
});
