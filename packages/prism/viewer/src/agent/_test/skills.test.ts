/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { withChange } from '../skills.js';

// The declarations (agent.md §3): plain data a spec extends without the toolkit.

describe('the agent declarations', () => {
  it('import nothing, so a spec that extends them does not bring the toolkit or zod', () => {
    const source = readFileSync(new URL('../skills.ts', import.meta.url), 'utf8');
    expect(source).not.toMatch(/^\s*import\b/m);
    expect(source).not.toMatch(/^\s*export\b.*\bfrom\s+['"]/m);
    expect(source).not.toMatch(/\bimport\(/);
  });

  it('apply a change by id: removed, replaced in place, added at the end', () => {
    const defaults = [{ id: 'a' }, { id: 'b', v: 1 }, { id: 'c' }];
    expect(withChange(defaults, undefined)).toBe(defaults);
    expect(withChange(defaults, { add: [{ id: 'd' }, { id: 'b', v: 2 }], remove: ['a'] })).toEqual([
      { id: 'b', v: 2 },
      { id: 'c' },
      { id: 'd' },
    ]);
  });
});
