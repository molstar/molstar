/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createAppSpec } from '../spec.js';

// The app's spec from its URL (app.md §4): Home at start only when the URL loads nothing; the light restore only with
// `?restore` (journal.md §6).

describe('the app spec', () => {
  it('opens Home at start when the URL loads nothing, and keeps the restore only when asked', () => {
    expect(createAppSpec('').options?.ui?.homeWhenEmpty).toBe(true);
    expect(createAppSpec('?pdb=1cbs').options?.ui?.homeWhenEmpty).toBe(false);
    expect(createAppSpec('').persist?.restore).toBeUndefined();
    expect(createAppSpec('?restore').persist?.restore).toBeDefined();
  });

  it("reads parameter names case-insensitively, as v5's links are written", () => {
    expect(createAppSpec('?PDB=1cbs').options?.ui?.homeWhenEmpty).toBe(false);
    expect(createAppSpec('?Restore').persist?.restore).toBeDefined();
  });
});
