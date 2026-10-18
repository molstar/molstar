/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { UUID } from '../uuid.js';

describe('UUID.create22', () => {
  it('makes 22 URL-safe base64 characters', () => {
    for (let i = 0; i < 200; i++) expect(UUID.create22()).toMatch(/^[A-Za-z0-9_-]{22}$/);
  });

  it('is unique over many calls, across the counter carries', () => {
    const ids = new Set<string>();
    const n = 64 * 64 * 3 + 17;
    for (let i = 0; i < n; i++) ids.add(UUID.create22());
    expect(ids.size).toBe(n);
  });
});
