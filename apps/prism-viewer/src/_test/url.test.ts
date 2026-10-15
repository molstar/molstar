/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { unknownUrlParameters } from '../url.js';

describe('unknownUrlParameters', () => {
  const load = ['pdb', 'url', 'mvs-url'];

  it('names the parameters neither list knows, once each', () => {
    expect(unknownUrlParameters('?pdb=1cbs&hide-controls=1&foo=1&bar&foo=2&PIXEL-SCALE=2', load)).toEqual([
      'foo',
      'bar',
    ]);
  });

  it('knows the options, the dropped ones and the load parameters', () => {
    expect(unknownUrlParameters('?pick-padding=1&map-provider=rcsb&mvs-url=x&debug-mode=1', load)).toEqual([]);
    expect(unknownUrlParameters('', load)).toEqual([]);
  });
});
