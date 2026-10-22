/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { getFontAtlas, setCanvasModule } from '../font-atlas.js';
import { TextBuilder } from '../text-builder.js';

// A font of 40 px glyphs where 'A' then 'V' kern by -12 px and 'j' reaches 6 px left of its pen; it draws nothing.
beforeAll(() => {
  const context = {
    measureText: (text: string) => {
      const width = 40 * text.length - (text.includes('AV') ? 12 : 0);
      return { width, actualBoundingBoxLeft: text === 'j' ? 6 : 0, actualBoundingBoxRight: width };
    },
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});
afterAll(() => setCanvasModule(undefined));

describe('TextBuilder', () => {
  // a font no other test asks for: the atlas is drawn with this file's stand-in context
  const props = {
    fontFamily: 'monospace',
    fontStyle: 'oblique',
    fontVariant: 'small-caps',
    fontWeight: 'bold',
    fontQuality: 0,
    attachment: 'bottom-left',
  } as const;

  /** The left edge of each glyph quad, in lineheights. */
  function glyphLefts(str: string) {
    const builder = TextBuilder.create(props, 16, 16);
    builder.add(str, 0, 0, 0, 0, 1, 0);
    const text = builder.getText();
    const m = text.mappingBuffer.ref.value;
    const lefts: number[] = [];
    for (let i = 0; i < text.charCount; ++i) lefts.push(m[i * 8]);
    return lefts;
  }

  it('places glyphs by kerned advances, the cell starting its pen offset before them', () => {
    const atlas = getFontAtlas(props);
    const lh = atlas.lineHeight;
    const pen = atlas.get('A').left;
    expect(pen).toBeCloseTo(atlas.buffer / lh);

    const [a, v] = glyphLefts('AV');
    expect(a).toBeCloseTo(-pen);
    expect(v).toBeCloseTo(28 / lh - pen); // 40 less the pair's kerning

    const [, x] = glyphLefts('AX');
    expect(x).toBeCloseTo(40 / lh - pen);
  });

  it('fits ink left of the pen in the cell', () => {
    const atlas = getFontAtlas(props);
    const j = atlas.get('j');
    expect(j.left).toBeCloseTo((atlas.buffer + 6) / atlas.lineHeight);
    expect(j.w).toBe(Math.ceil(6 + 40 + 2 * atlas.buffer));
    const [, jLeft] = glyphLefts('Aj');
    expect(jLeft).toBeCloseTo(40 / atlas.lineHeight - j.left);
  });

  it('draws no quad for a space, which still advances', () => {
    const atlas = getFontAtlas(props);
    const lefts = glyphLefts('A X');
    expect(lefts.length).toBe(2);
    expect(lefts[1]).toBeCloseTo(80 / atlas.lineHeight - atlas.get('X').left);
  });
});
