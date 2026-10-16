/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import { describe, expect, it } from 'vitest';
import { crambin, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// `load.files` judges each file by what it made: a file the plugin cannot parse fails, and leaves no node behind.

const broken = () => new File(['this is not a CIF file'], 'broken.cif');

function rootChildren(v: PrismViewer): number {
  return v.plugin.state.data.tree.children.get(StateTransform.RootRef).size;
}

describe('load.files', () => {
  it('fails for a file that does not parse', async () => {
    const v = await headlessViewer();
    const r = await v.api.load.files({ files: [broken()] });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    if (!r.ok) expect(r.message).toMatch(/broken\.cif/);
    expect(rootChildren(v)).toBe(0);
    expect(v.managers.history.steps.value).toHaveLength(0);
  });

  it('opens the files that parse, warns about the rest, and leaves no raw data behind', async () => {
    const v = await headlessViewer();
    const r = await v.api.load.files({ files: [crambin(), broken()] });
    expectOk(r);
    expect(r.did).toMatch(/Opened 1crn\.cif \(1 structure\)/);
    expect(r.did).not.toMatch(/broken/);
    expect(r.warnings).toEqual(['Could not open broken.cif.']);
    expect(v.api.structure.list.value).toHaveLength(1);
    expect(rootChildren(v)).toBe(1);
  });
});
