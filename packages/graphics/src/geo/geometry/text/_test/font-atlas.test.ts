/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterAll, afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  FontFamilyFonts,
  getFontAtlas,
  getFontAtlasFont,
  loadFontAtlasFonts,
  setCanvasModule,
  setFontFamilyFonts,
  withFontAtlasFonts,
} from '../font-atlas.js';
import { Text } from '../text.js';

// Node has no canvas: a stand-in context that draws nothing and, as a canvas does, ignores a font it cannot parse
// (here: one naming 'Invalid'). Atlases are compared by identity (`a === b`): `not.toBe` deep-compares the two on
// success, megabytes of texture.
const fontsSet: string[] = [];
function createContext() {
  let font = '10px sans-serif';
  return {
    get font() {
      return font;
    },
    set font(value: string) {
      fontsSet.push(value);
      if (!value.includes('Invalid')) font = value;
    },
    measureText: () => ({ width: 8 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
}
beforeAll(() => {
  setCanvasModule({ createCanvas: () => ({ getContext: createContext }) });
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
    expect(mono === serif).toBe(false);
    getFontAtlas({ ...font, fontFamily: 'cursive' });
    getFontAtlas({ ...font, fontFamily: 'serif', fontWeight: 'bold' });
    // serif is the oldest entry; using it makes monospace the least recently used
    expect(getFontAtlas({ ...font, fontFamily: 'serif' })).toBe(serif);
    getFontAtlas({ ...font, fontFamily: 'serif', fontStyle: 'italic' });
    expect(getFontAtlas({ ...font, fontFamily: 'serif' })).toBe(serif);
    expect(getFontAtlas({ ...font, fontFamily: 'monospace' }) === mono).toBe(false);
  });
});

describe('font family table and web fonts', () => {
  const font = { fontQuality: 0, fontFamily: 'serif' as const };
  const serif = FontFamilyFonts.serif;
  afterEach(() => {
    setFontFamilyFonts('serif', serif);
    vi.unstubAllGlobals();
  });

  /** A `document.fonts` whose faces load when `arrive()` is called, or fail to with `fail()`. */
  function stubFontFaceSet() {
    let loaded = false;
    let arrive = () => {};
    let fail = () => {};
    const arrived = new Promise<void>((resolve, reject) => {
      arrive = () => {
        loaded = true;
        resolve();
      };
      fail = () => reject(new Error('404'));
    });
    const fonts = { check: vi.fn(() => loaded), load: vi.fn(() => arrived.then(() => [])) };
    vi.stubGlobal('document', { fonts });
    const settled = async () => {
      await arrived.catch(() => {});
      for (let i = 0; i < 4; ++i) await Promise.resolve();
    };
    return { fonts, arrive: async () => (arrive(), await settled()), fail: async () => (fail(), await settled()) };
  }

  it('draws a family with the fonts of its table entry', () => {
    expect(getFontAtlasFont({ fontFamily: 'sans-serif', fontQuality: 0 })).toBe(
      `normal normal normal 64px ${FontFamilyFonts['sans-serif']}`,
    );
    setFontFamilyFonts('serif', '"Test Serif", serif');
    expect(getFontAtlasFont({ ...font, fontWeight: 'bold' })).toBe('normal normal bold 64px "Test Serif", serif');
    // a family the table does not have (an old state, a string from MVS) is used as is
    expect(getFontAtlasFont({ fontFamily: 'Inter' as any, fontQuality: 0 })).toBe('normal normal normal 64px Inter');
  });

  it('has nothing to wait for without document.fonts (Node)', () => {
    expect(loadFontAtlasFonts(font)).toBeUndefined();
    expect(withFontAtlasFonts(font, () => 1)).toBe(1);
  });

  it('waits for a web font, then builds', async () => {
    setFontFamilyFonts('serif', '"Lazy Serif A", serif');
    const { fonts, arrive } = stubFontFaceSet();
    const build = vi.fn(() => 'built');
    const built = withFontAtlasFonts(font, build);
    expect(built).toBeInstanceOf(Promise);
    await Promise.resolve();
    expect(build).not.toHaveBeenCalled();
    // one download per font, however many wait for it
    expect(loadFontAtlasFonts(font)).toBe(loadFontAtlasFonts(font));
    expect(fonts.load).toHaveBeenCalledTimes(1);
    expect(fonts.load).toHaveBeenCalledWith('normal normal normal 64px "Lazy Serif A", serif');
    await arrive();
    expect(await built).toBe('built');
    expect(withFontAtlasFonts(font, () => 'now')).toBe('now');
  });

  it('replaces an atlas drawn with fallback glyphs once the font arrives', async () => {
    setFontFamilyFonts('serif', '"Lazy Serif B", serif');
    const { arrive } = stubFontFaceSet();
    const fallback = getFontAtlas(font);
    expect(getFontAtlas(font)).toBe(fallback);
    await arrive();
    const loaded = getFontAtlas(font);
    expect(loaded === fallback).toBe(false);
    expect(getFontAtlas(font)).toBe(loaded);
  });

  it('builds a visual that waited with the font, though a text built meanwhile drew a fallback atlas', async () => {
    setFontFamilyFonts('serif', '"Lazy Serif C", serif');
    const { arrive } = stubFontFaceSet();
    const waited = withFontAtlasFonts(font, () => getFontAtlas(font));
    const fallback = getFontAtlas(font); // the camera helper, say: does not wait
    await arrive();
    expect((await waited) === fallback).toBe(false);
  });

  it('keeps the fallback of a font that fails to download, and stops waiting for it', async () => {
    setFontFamilyFonts('serif', '"Missing Serif", serif');
    const { fonts, fail } = stubFontFaceSet();
    const waited = withFontAtlasFonts(font, () => 'fallback');
    const fallback = getFontAtlas(font);
    await fail();
    expect(await waited).toBe('fallback');
    expect(getFontAtlas(font)).toBe(fallback);
    expect(withFontAtlasFonts(font, () => 'now')).toBe('now');
    expect(fonts.load).toHaveBeenCalledTimes(1);
  });

  it('stops waiting after 3 s; the font replaces the fallback when it arrives', async () => {
    vi.useFakeTimers();
    try {
      setFontFamilyFonts('serif', '"Slow Serif", serif');
      const { arrive } = stubFontFaceSet();
      const waited = withFontAtlasFonts(font, () => getFontAtlas(font));
      await vi.advanceTimersByTimeAsync(3000);
      const fallback = await waited;
      expect(getFontAtlas(font)).toBe(fallback);
      await arrive();
      expect(getFontAtlas(font) === fallback).toBe(false);
    } finally {
      vi.useRealTimers();
    }
  });

  it('draws the generic family when the table has a font list the canvas cannot parse', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    setFontFamilyFonts('serif', 'Invalid 4, serif');
    fontsSet.length = 0;
    getFontAtlas(font);
    expect(fontsSet).toEqual(['normal normal normal 64px Invalid 4, serif', 'normal normal normal 64px serif']);
    expect(warn).toHaveBeenCalledTimes(1);
  });
});
