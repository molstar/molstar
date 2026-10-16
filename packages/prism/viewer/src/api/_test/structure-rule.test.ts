/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer } from '../../_test/harness.js';
import { findStructure } from '../../sources/naming.js';

// The structure an operation acts on when none is named (targets.md §1): the one being worked with (focused, else the
// one holding the selection), else the only one; a target falls back to every one; a destructive operation needs a
// name when several are loaded.

describe('the structure, when none is named', () => {
  it('is the one holding the selection, for an operation on one structure and for a target', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const move = { translation: [1, 0, 0] as never };
    expect(await v.api.structure.transform(move)).toMatchObject({ ok: false, code: 'ambiguous' });
    expectOk(await v.api.select.set({ target: 'A:10-12', structure: '1CRN' }));
    const moved = await v.api.structure.transform(move);
    expectOk(moved);
    const structures = v.plugin.managers.structure.hierarchy.current.structures;
    const crn = findStructure(structures, '1CRN')!;
    expect(moved.changed.nodes).toEqual([crn.cell.transform.ref]);
    // A target without a structure: the one holding the selection too.
    const colored = await v.api.representation.color({ target: 'polymer', color: 'red' });
    expectOk(colored);
    expect(colored.did).toMatch(/1CRN/i);
  });

  it('is never the selected one for a destructive operation', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    expectOk(await v.api.select.set({ target: 'A:10-12', structure: '1CRN' }));
    expect(await v.api.structure.remove({})).toMatchObject({ ok: false, code: 'ambiguous' });
    expectOk(await v.api.structure.remove({ structure: '1CRN' }));
    // One left: it is the only one.
    expectOk(await v.api.structure.remove({}));
  });
});
