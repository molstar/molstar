/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, representationTypes } from '../../_test/harness.js';

// Redo (v5's): what Undo took back comes back, until a new step.

describe('history.redo', () => {
  it('puts back what Undo took back, and is blocked after a new step', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = representationTypes(v);
    expectOk(await v.api.styles.apply({ style: 'surface' }));
    const surface = representationTypes(v);
    expect((await v.api.history.redo({})).ok).toBe(false);
    expectOk(await v.api.history.undo({}));
    expect(representationTypes(v)).toEqual(before);
    expectOk(await v.api.history.redo({}));
    expect(representationTypes(v)).toEqual(surface);
    expect(v.managers.history.steps.value[0].label).toMatch(/Surface/);
    expectOk(await v.api.history.undo({}));
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    const r = await v.api.history.redo({});
    expect(r.ok).toBe(false);
  });
});
