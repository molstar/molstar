/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { parseNumber } from '../params-form.js';

// A typed number in a form: a point or a comma as the decimal point, whatever the browser's locale.

describe('a typed number', () => {
  it('reads a point or a comma as the decimal point, and nothing else as a number', () => {
    expect(parseNumber('0.5')).toBe(0.5);
    expect(parseNumber('0,5')).toBe(0.5);
    expect(parseNumber(' -1 ')).toBe(-1);
    expect(parseNumber('')).toBeNaN();
    expect(parseNumber('abc')).toBeNaN();
  });
});
