/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { LanguageVersion } from '@molstar/seqx-core';
import { ViewerLanguageVersion } from '../index.js';

describe('the viewer entry', () => {
  it('draws the core language version, in Node', () => {
    expect(ViewerLanguageVersion).toBe(LanguageVersion);
  });
});
