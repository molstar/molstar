/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Color } from '../color.js';
import { checkColorVision, simulateColorVision } from '../cvd.js';
import { ColorLists } from '../lists.js';
import {
  AminoAcidNames,
  DefaultResidueHydrophobicity,
  ResidueHydrophobicity,
  ResidueSchemes,
} from '../residue-schemes.js';
import { Oklab } from '../spaces/oklab.js';

describe('oklab', () => {
  it('reads white, black and a primary as the reference does', () => {
    const w = Oklab.fromColor(Oklab.zero(), Color(0xffffff));
    expect(w[0]).toBeCloseTo(1, 4);
    expect(Math.abs(w[1]) + Math.abs(w[2])).toBeLessThan(1e-4);
    expect(Oklab.fromColor(Oklab.zero(), Color(0x000000))[0]).toBeCloseTo(0, 6);
    // Ottosson's published sRGB red: L 0.62796, a 0.22486, b 0.12585
    const r = Oklab.fromColor(Oklab.zero(), Color(0xff0000));
    expect(r[0]).toBeCloseTo(0.62796, 3);
    expect(r[1]).toBeCloseTo(0.22486, 3);
    expect(r[2]).toBeCloseTo(0.12585, 3);
  });

  it('round-trips every byte color it is given', () => {
    for (const c of [0x123456, 0xabcdef, 0x7a8fa6, 0xff00ff, 0x00ff00, 0x808080]) {
      expect(Oklab.toColor(Oklab.fromColor(Oklab.zero(), Color(c)))).toBe(c);
    }
  });

  it('interpolates with its ends exact and lightness monotone between', () => {
    const a = Color(0x0d0887);
    const b = Color(0xf0f921);
    expect(Color.interpolateOklab(a, b, 0)).toBe(a);
    expect(Color.interpolateOklab(a, b, 1)).toBe(b);
    let prev = -1;
    for (let i = 0; i <= 10; i++) {
      const l = Oklab.lightness(Color.interpolateOklab(a, b, i / 10));
      expect(l).toBeGreaterThan(prev);
      prev = l;
    }
  });

  it('maps lightness into a band, hue kept, in gamut', () => {
    const light = Color(0xcfe0f5);
    const dark = Oklab.mapLightness(light, [0.8, 1], [0.3, 0.45]);
    const l = Oklab.lightness(dark);
    expect(l).toBeGreaterThanOrEqual(0.29);
    expect(l).toBeLessThanOrEqual(0.46);
    const h = (c: Color) => {
      const o = Oklab.fromColor(Oklab.zero(), c);
      return Math.atan2(o[2], o[1]);
    };
    expect(Math.abs(h(dark) - h(light))).toBeLessThan(0.15);
    const saturated = Oklab.withLightness(Color(0x00ff00), 0.2);
    expect(Oklab.inGamut(Oklab.fromColor(Oklab.zero(), saturated))).toBe(true);
  });
});

describe('color vision', () => {
  it('keeps greys grey and confuses red and green under deuteranopia', () => {
    expect(simulateColorVision(Color(0x808080), 'protanopia')).toBe(0x808080);
    expect(simulateColorVision(Color(0x123456), 'normal')).toBe(0x123456);
    const red = Color(0xd62728);
    const green = Color(0x2ca02c);
    const normal = Oklab.colorDistance(red, green);
    const deut = Oklab.colorDistance(
      simulateColorVision(red, 'deuteranopia'),
      simulateColorVision(green, 'deuteranopia'),
    );
    expect(deut).toBeLessThan(normal / 2);
  });

  it('passes Okabe–Ito and fails a red-green pair', () => {
    // its seven hues; its grey sits 0.02 from its reddish purple under deuteranopia
    const okabe = ColorLists['okabe-ito'].list.map(Color.fromColorListEntry);
    expect(checkColorVision(okabe.slice(0, 7), 0.06).ok).toBe(true);
    expect(checkColorVision(okabe, 0.06).worst).toMatchObject({ vision: 'deuteranopia', i: 6, j: 7 });
    const check = checkColorVision([Color(0xd62728), Color(0x8c8c00)], 0.06);
    expect(check.ok).toBe(false);
    expect(check.worst.vision).not.toBe('normal');
    expect(checkColorVision([Color(0)]).worst.distance).toBe(Number.POSITIVE_INFINITY);
  });
});

describe('residue schemes', () => {
  it('colors every amino acid of each protein scheme with a hex color', () => {
    for (const s of Object.values(ResidueSchemes)) {
      expect(s.provenance.length).toBeGreaterThan(10);
      for (const c of Object.values(s.colors)) expect(c).toMatch(/^#[0-9a-f]{6}$/);
      if (s.alphabet === 'protein') for (const l of Object.keys(AminoAcidNames)) expect(s.colors).toHaveProperty(l);
      else for (const l of 'ACGT') expect(s.colors).toHaveProperty(l);
    }
  });

  it('holds the hydrophobicity scale for every amino acid (Mol*’s one table)', () => {
    for (const name of Object.values(AminoAcidNames)) expect(ResidueHydrophobicity).toHaveProperty(name);
    expect(ResidueHydrophobicity.TRP[0]).toBe(-1.85);
    expect(DefaultResidueHydrophobicity).toEqual([0, 0, 0]);
  });

  it('marks as safe only what passes the check', () => {
    for (const s of Object.values(ResidueSchemes)) {
      const colors = [...new Set(Object.values(s.colors))].map((c) => Color.fromHexStyle(c));
      if (s.cvdSafe) expect(checkColorVision(colors, 0.06).ok, s.name).toBe(true);
    }
  });
});
