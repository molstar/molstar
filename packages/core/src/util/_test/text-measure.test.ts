/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { TextMeasure, splitGraphemes } from '../text-measure.js';

/** A font where every grapheme is 10 wide and 'A' then 'V' kern by -3. */
function kerningFont() {
  return vi.fn((text: string) => {
    let w = 10 * splitGraphemes(text).length;
    for (let i = 0; i < text.length - 1; ++i) if (text[i] === 'A' && text[i + 1] === 'V') w -= 3;
    return w;
  });
}

describe('splitGraphemes', () => {
  it('splits ASCII by code unit', () => {
    expect(splitGraphemes('ALA 12')).toEqual(['A', 'L', 'A', ' ', '1', '2']);
    expect(splitGraphemes('')).toEqual([]);
  });

  it('keeps surrogate pairs, combining marks and emoji sequences whole', () => {
    expect(splitGraphemes('4.2 Å')).toEqual(['4', '.', '2', ' ', 'Å']);
    expect(splitGraphemes('ét')).toEqual(['é', 't']);
    expect(splitGraphemes('a𝛼b')).toEqual(['a', '𝛼', 'b']);
    expect(splitGraphemes('👍🏽!')).toEqual(['👍🏽', '!']);
  });
});

describe('TextMeasure', () => {
  it('places each grapheme by the pair it closes, kerning included', () => {
    const m = new TextMeasure(kerningFont());
    const { graphemes, offsets, width } = m.layout('AVA');
    expect(graphemes).toEqual(['A', 'V', 'A']);
    expect(Array.from(offsets)).toEqual([0, 7, 17]);
    expect(width).toBe(27);
    expect(m.advance('V', 'A')).toBe(7);
    expect(m.advance('V')).toBe(10);
  });

  it('measures each width once', () => {
    const measure = kerningFont();
    const m = new TextMeasure(measure);
    m.layout('AVAVAV');
    const calls = measure.mock.calls.length;
    m.layout('AVAVAV');
    m.layout('VAVA');
    expect(measure.mock.calls.length).toBe(calls);
    m.clear();
    m.layout('AV');
    expect(measure.mock.calls.length).toBeGreaterThan(calls);
  });
});
