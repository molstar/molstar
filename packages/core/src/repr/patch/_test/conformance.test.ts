/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { applyPatch, type PatchResult } from '../../patch.js';
import { hashNode } from '../../tree/index.js';
import { checkPatches, PatchCases } from '../conformance.js';

describe('repr/patch: the conformance suite', () => {
  it('applyPatch holds every case, and every inverse gives the root back', () => {
    expect(checkPatches(applyPatch)).toEqual([]);
    expect(PatchCases.length).toBeGreaterThan(50);
  });

  it('finds an apply that refuses nothing, or loses the inverse', () => {
    const lenient = (...args: Parameters<typeof applyPatch>): PatchResult => {
      const r = applyPatch(...args);
      return r.ok
        ? r
        : {
            ok: true,
            changed: false,
            root: args[0],
            revision: hashNode(args[0]),
            changes: { made: [], params: [], children: [], removed: [], moved: [] },
            issues: [],
          };
    };
    expect(checkPatches(lenient).length).toBeGreaterThan(20);
    const noInverse = (...args: Parameters<typeof applyPatch>): PatchResult => {
      const r = applyPatch(...args);
      return r.ok ? { ...r, inverse: undefined } : r;
    };
    expect(checkPatches(noInverse).some((f) => /no inverse/.test(f))).toBe(true);
  });
});
