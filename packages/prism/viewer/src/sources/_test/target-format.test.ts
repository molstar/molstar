/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../_test/harness.js';
import { createTargetFormat } from '../target-format.js';

// The command menu's `target` format (targets.md §1, parity.md §3).

describe('the target format', () => {
  it('checks targets as typed and completes keywords, chains and ligands', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const format = createTargetFormat(v.plugin);
    const ctx = { empty: false };
    expect(format.check('A:10-20', ctx)).toMatchObject({ ok: true });
    expect(format.check('resn (((', ctx)).toMatchObject({ ok: false });
    const all = format.complete('', ctx).map((c) => c.value);
    expect(all).toContain('ligand');
    expect(all).toContain('selection');
    expect(all).toContain('A');
    expect(all).toContain('resn REA');
    expect(format.complete('re', ctx).map((c) => c.value)).toContain('resn REA');
  });
});
