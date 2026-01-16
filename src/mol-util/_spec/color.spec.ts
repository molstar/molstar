
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Zach Charlop-Powers <zach.charlop.powers@gmail.com>
 */

import { decodeColor } from '../color/utils';
import { Color } from '../color/color';

describe('decodeColor', () => {
  it('handles 3-character hex codes', () => {
    const color = decodeColor('#f0a');
    expect(color).toBe(Color.fromHexStyle('#ff00aa'));
  });

  it('handles 4-character hex codes (ignores alpha)', () => {
    const color = decodeColor('#f0ab');
    expect(color).toBe(Color.fromHexStyle('#ff00aa'));
  });

  it('handles 6-character hex codes', () => {
    const color = decodeColor('#ff00aa');
    expect(color).toBe(Color.fromHexStyle('#ff00aa'));
  });

  it('handles 8-character hex codes (strips alpha)', () => {
    const color = decodeColor('#ff00aacc');
    expect(color).toBe(Color.fromHexStyle('#ff00aa'));
  });

  it('returns undefined for invalid hex lengths', () => {
    // 2 characters (1 hex digit)
    expect(decodeColor('#f')).toBeUndefined();
    // 3 characters (2 hex digits)
    expect(decodeColor('#ff')).toBeUndefined();
    // 6 characters (5 hex digits) - invalid length
    expect(decodeColor('#ff00a')).toBeUndefined();
    // 8 characters (7 hex digits) - invalid length
    expect(decodeColor('#ff00aab')).toBeUndefined();
    // 10 characters (9 hex digits) - too long
    expect(decodeColor('#ff00aabbcc')).toBeUndefined();
  });

  it('returns undefined for invalid input', () => {
    expect(decodeColor(undefined)).toBeUndefined();
    expect(decodeColor(null)).toBeUndefined();
    expect(decodeColor('invalid')).toBeUndefined();
    expect(decodeColor('notahex')).toBeUndefined();
    expect(decodeColor('#gggggg')).toBeUndefined();
  });

  it('handles case insensitivity', () => {
    const lowerCase = decodeColor('#ff00aa');
    const upperCase = decodeColor('#FF00AA');
    const mixedCase = decodeColor('#Ff00Aa');
    expect(lowerCase).toBe(upperCase);
    expect(upperCase).toBe(mixedCase);
  });
});