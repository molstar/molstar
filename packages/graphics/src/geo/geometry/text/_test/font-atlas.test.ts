/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { getFontAtlas, setCanvasModule } from '../font-atlas.js';
import { Text } from '../text.js';

// Node has no canvas: a stand-in context that draws nothing.
beforeAll(() => {
  const context = {
    measureText: () => ({ width: 8 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});
afterAll(() => setCanvasModule(undefined));

describe('getFontAtlas', () => {
  // the smallest atlas (a few MB), so the test does not allocate 43 MB per atlas
  const font = { fontQuality: 0 };

  it('shares one atlas between props that differ outside the font fields (08/B6)', () => {
    const base = { ...PD.getDefaultValues(Text.Params), ...font };
    const a = getFontAtlas(base);
    expect(getFontAtlas({ ...base, sizeFactor: 2, backgroundMargin: 0.5, tetherLength: 3 })).toBe(a);
    expect(getFontAtlas({ ...base, items: [{ text: 'x', position: [0, 0, 0] }] } as any)).toBe(a);
    expect(getFontAtlas(font)).toBe(a);
  });

  it('keeps only the font fields', () => {
    const atlas = getFontAtlas({ ...PD.getDefaultValues(Text.Params), ...font, fontWeight: 'bold' });
    expect(Object.keys(atlas.props).sort()).toEqual([
      'fontFamily',
      'fontQuality',
      'fontStyle',
      'fontVariant',
      'fontWeight',
    ]);
  });

  it('builds a new atlas for each font and keeps the 4 most recently used', () => {
    const serif = getFontAtlas({ ...font, fontFamily: 'serif' });
    const mono = getFontAtlas({ ...font, fontFamily: 'monospace' });
    expect(mono).not.toBe(serif);
    getFontAtlas({ ...font, fontFamily: 'cursive' });
    getFontAtlas({ ...font, fontFamily: 'serif', fontWeight: 'bold' });
    // serif is the oldest entry; using it makes monospace the least recently used
    expect(getFontAtlas({ ...font, fontFamily: 'serif' })).toBe(serif);
    getFontAtlas({ ...font, fontFamily: 'serif', fontStyle: 'italic' });
    expect(getFontAtlas({ ...font, fontFamily: 'serif' })).toBe(serif);
    expect(getFontAtlas({ ...font, fontFamily: 'monospace' })).not.toBe(mono);
  });
});
