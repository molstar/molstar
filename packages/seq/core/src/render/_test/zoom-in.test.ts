/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer whose wants together pass the store's bound on a layer's tiles or read cost (_spec/README.md "Render", Data):
// its band hatched with "zoom in to read this track", an info issue said, distinct from pending (which settles by
// itself); zoomed in, the layer reads and draws again. The drawn list is what the canvas and the SVG both draw.

import { describe, expect, it } from 'vitest';
import type { Frame } from '../index.js';
import { ZoomInText } from '../tracks-draw.js';
import { BigRuntime, bigSpec } from './fixtures.js';
import { settle } from './helpers.js';

function texts(f: Frame): string[] {
  const out: string[] = [];
  for (const p of Object.values(f.lists.passes))
    for (let i = 0; i < p.texts.n; i++) out.push(f.lists.string(p.texts.str[i]));
  return out;
}

const Size = { width: 720, height: 300 };

describe('the zoom-in state', () => {
  it('a layer past maxTilesPerLayer draws its band hatched with “zoom in to read this track”, said once', async () => {
    // 2 tiles a layer at most: 30 kb at its level (16 positions a bin, 16 kb a tile) with its margins are more
    const s = await settle(bigSpec(), [BigRuntime], {}, {}, { maxTilesPerLayer: 2, summaries: false });
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 30_000 } });
    expect(texts(f)).toContain(ZoomInText);
    expect(f.layout.issues.some((i) => i.code === 'seq.render.zoom_in')).toBe(true);
    expect(s.store.issues.value.some((i) => i.code === 'seq.store.too_many_tiles')).toBe(true);
    // zoomed in, the layers read again: no text, no issue
    const near = await s.frame({ size: Size, viewport: { start: 1000, end: 1400 } });
    expect(texts(near)).not.toContain(ZoomInText);
    expect(near.layout.issues.some((i) => i.code === 'seq.render.zoom_in')).toBe(false);
  });

  it('a pending band says nothing (it settles by itself)', async () => {
    const s = await settle(bigSpec(), [BigRuntime]);
    const f = s.renderer.frame({ size: Size, viewport: { start: 0, end: 3_000_000 } });
    expect(texts(f)).not.toContain(ZoomInText);
  });
});
