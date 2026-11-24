/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every built-in layer type declares its fallback (_spec/README.md "Resolved view", Layer kinds; extensibility, EC-235):
// the language's table names one for each (`none` ends the chain), the renderer's built-ins are the table's names, and a
// view resolves each type's chain as declared.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { LayerTypes, resolveView } from '@molstar/seqx-core/view';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { BuiltInLayerTypes } from '../built-ins.js';

const Declared: Readonly<Record<string, string>> = {
  letters: 'none',
  ruler: 'none',
  ranges: 'none',
  values: 'none',
  secondary: 'ranges',
  variants: 'ranges',
  frame: 'none',
  rows: 'none',
  logo: 'none',
  bands: 'ranges',
  pairs: 'ranges',
  grid: 'none',
  matrix: 'values',
};

describe('the built-in layer types’ declared fallbacks', () => {
  it('the renderer’s built-ins and the language’s table are the same types', () => {
    expect(Object.keys(BuiltInLayerTypes).sort()).toEqual(Object.keys(LayerTypes).sort());
    expect(Object.keys(Declared).sort()).toEqual(Object.keys(LayerTypes).sort());
  });

  for (const [name, fallback] of Object.entries(Declared)) {
    it(`EC-74: ${name} declares ${fallback}, a known type or none, and its chain ends`, () => {
      expect(LayerTypes[name].fallback).toBe(fallback);
      if (fallback !== 'none') expect(LayerTypes[fallback]).toBeDefined();
      // the chain, as the resolved view walks it
      const chain = [name];
      for (let t = LayerTypes[name].fallback; t !== 'none' && !chain.includes(t); t = LayerTypes[t].fallback)
        chain.push(t);
      expect(chain.length).toBeLessThanOrEqual(2);
      expect(chain.at(-1)).toBe(fallback === 'none' ? name : fallback);
    });
  }

  it('a view resolves a layer’s chain from the declarations', async () => {
    const b = createSeqSpec({ title: 'chains' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGTAC' });
    b.views.view('v', (v) => {
      const section = v.section('sec', { of: s });
      section.track('t').layer('letters', { type: 'letters' });
      section.track('r').layer('ruler', { type: 'ruler' });
    });
    const model = await resolveSpec(b.build());
    const view = resolveView(model, '/views/v');
    const chains = view.sections[0].tracks.flatMap((t) => t.layers.map((l) => l.fallbacks));
    expect(chains).toEqual([['letters'], ['ruler']]);
  });
});
