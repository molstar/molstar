/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { LanguageVersion } from '../index.js';

describe('the core entry', () => {
  it('names a semver language version', () => {
    expect(LanguageVersion).toMatch(/^\d+\.\d+\.\d+$/);
  });
});
