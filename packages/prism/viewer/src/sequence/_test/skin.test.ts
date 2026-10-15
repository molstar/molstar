/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { inkForColor } from '@molstar/prism-ui';
import { describe, expect, it } from 'vitest';
import { legendOf } from '../legend.js';
import { skinFor } from '../skin.js';

describe('skinFor', () => {
  it('draws on the ink: the surface is the canvas color, the theme sets the flag, ink and accent', () => {
    const light = skinFor(inkForColor('#f4f5f7'));
    expect(light.dark).toBe(false);
    expect(light.surface).toBe('rgb(244 245 247)');
    expect(light.ink).toBe('rgb(51 43 31)');
    const dark = skinFor(inkForColor('#111318'));
    expect(dark.dark).toBe(true);
    expect(dark.surface).toBe('rgb(17 19 24)');
    expect(dark.ink).toBe('rgb(204 212 224)');
    expect(dark.accent).toBe('rgb(81 162 251)');
  });

  it('writes the selection as the accent at 25% and takes fonts when given', () => {
    const s = skinFor(inkForColor('#ffffff'), { mono: 'Mono', ui: 'Sans' });
    expect(s.selection).toBe('rgb(174 93 4 / 0.25)');
    expect([s.fontMono, s.fontUi]).toEqual(['Mono', 'Sans']);
  });

  it('draws what the structure does not model in a quiet text tone, between the ink and the surface', () => {
    expect(skinFor(inkForColor('#ffffff')).unobserved).toBe('rgb(163 160 154)');
    expect(skinFor(inkForColor('#111318')).unobserved).toBe('rgb(101 106 114)');
  });

  it('is the same for the same ink', () => {
    expect(skinFor(inkForColor('#222222'))).toEqual(skinFor(inkForColor('#222222')));
  });
});

describe('the residue palette', () => {
  it('is the default scheme, with a legend of the classes for a light and a dark ink', () => {
    const light = legendOf('prism:residues', false);
    const dark = legendOf('prism:residues', true);
    expect(light.title).toBe('Prism residues');
    expect(light.items.map((i) => i.label)).toEqual([
      'Aliphatic',
      'Aromatic',
      'Positive',
      'Negative',
      'Polar',
      'Proline, glycine',
      'Cysteine',
    ]);
    expect(dark.items.map((i) => i.label)).toEqual(light.items.map((i) => i.label));
    expect(dark.items.map((i) => i.color)).not.toEqual(light.items.map((i) => i.color));
  });
});

describe('legendOf', () => {
  it("lists a letter scheme's groups, a ramp's ends and a categorical scheme's colors", () => {
    const zappo = legendOf('zappo', false);
    expect(zappo.title).toBe('Zappo');
    expect(zappo.items.map((i) => i.label)).toContain('aromatic');
    expect(zappo.ramps).toEqual([]);
    const viridis = legendOf('viridis', false);
    expect(viridis.ramps).toHaveLength(1);
    expect(viridis.ramps[0].stops.length).toBeGreaterThan(2);
    expect(legendOf('okabe_ito', true).items).toHaveLength(7);
  });

  it('is empty without a scheme', () => {
    expect(legendOf(null, false).items).toEqual([]);
    expect(legendOf('nope', false).items).toEqual([]);
  });
});
