/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Muted text holds 4.5:1 (_spec/README.md "Themes", Colors): the ruler's tick labels and its "not numbered" note, as the
// ruler draws them, against the ruler's well and the page, in both themes, from the theme tokens.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { contrastRgba } from '@molstar/seqx-core/theme/index';
import { describe, expect, it } from 'vitest';
import { DarkTokens, LightTokens, type RulerLay } from '../index.js';
import { MutedTextContrast, mutedInk } from '../tokens.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

describe('muted text holds 4.5:1', () => {
  for (const [name, tokens] of [
    ['light', LightTokens],
    ['dark', DarkTokens],
  ] as const) {
    it(`${name}: the ruler's ink holds ${MutedTextContrast}:1 on its well and the page`, () => {
      const wells = [tokens.colors['seq-track'], tokens.colors['msui-bg']];
      const ink = mutedInk(tokens, tokens.colors['seq-ruler'], wells);
      for (const w of wells) expect(contrastRgba(ink, w)).toBeGreaterThanOrEqual(MutedTextContrast);
    });
  }

  it('a muted color that already holds is kept; one that falls short is the text color at 75%', () => {
    const t = LightTokens;
    const wells = [t.colors['seq-track']];
    expect(mutedInk(t, t.colors['msui-text'], wells)).toBe(t.colors['msui-text']);
    const pale = t.colors['msui-line'];
    expect(contrastRgba(pale, wells[0])).toBeLessThan(MutedTextContrast);
    expect(mutedInk(t, pale, wells)).not.toBe(pale);
  });

  it('the light theme’s seq-ruler alone falls short, which is why the ruler raises it', () => {
    const l = LightTokens;
    expect(contrastRgba(l.colors['seq-ruler'], l.colors['seq-track'])).toBeLessThan(MutedTextContrast);
  });

  it('the ruler draws its tick labels and note in an ink that holds it', async () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'M'.repeat(200) });
    b.views.view('v', (v) => {
      const s = v.section('p', { of: p });
      s.track('ruler').layer('ruler', { type: 'ruler' });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: { width: 600, height: 120 } });
    const lay = layerLay<RulerLay>(f, 0);
    for (const ink of lay.inkColors)
      for (const well of [LightTokens.colors['seq-track'], LightTokens.colors['msui-bg']])
        expect(contrastRgba(ink, well)).toBeGreaterThanOrEqual(MutedTextContrast);
    expect(lay.inkColors[0] & 0xff).toBe(0xff);
  });
});
