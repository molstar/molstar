/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';

// Absent versus given (api.md §1): a value left out is kept as it is, a value given is set; no sentinel values
// (0, -1, '') stand for "keep". Switches take the state wanted; the commands flip them.

describe('absent values', () => {
  it('keeps what a unit-cell call leaves out, and shows the cell when given nothing', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const cell = () => v.plugin.managers.structure.hierarchy.current.structures[0].model?.unitcell?.cell;
    // A color alone restyles; it does not show the cell (its visibility is kept: hidden by the default preset).
    expectOk(await v.api.structure.unitCell({ color: 'red' }));
    expect(cell()?.state.isHidden ?? true).toBe(true);
    // Nothing given: shown (what the operation is for).
    expectOk(await v.api.structure.unitCell({}));
    expect(cell()?.state.isHidden).toBe(false);
    expect(await v.api.structure.unitCell({})).toMatchObject({ ok: true, unchanged: true });
  });

  it('asks how to move a structure given nothing, and an angle for an axis', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(await v.api.structure.transform({})).toMatchObject({ ok: false, code: 'invalid', message: 'Move how?' });
    expect(await v.api.structure.transform({ axis: [1, 0, 0] as never })).toMatchObject({
      ok: false,
      message: 'An axis needs an angle.',
    });
    expectOk(await v.api.structure.transform({ translation: [5, 0, 0] as never }));
    // A rotation given alone keeps no translation of its own: the transform is from where it was loaded.
    expectOk(await v.api.structure.transform({ axis: [0, 0, 1] as never, angle: 90 }));
  });

  it('flips a switch from its command, and takes the state wanted from the operation', async () => {
    const v = await headlessViewer();
    const nav = v.managers.navigation;
    expect(nav.expanded.value).toBe(false);
    await v.managers.commands.get('view.expand')!.run({}, { empty: false }, new AbortController().signal);
    expect(nav.expanded.value).toBe(true);
    await v.managers.commands.get('view.expand')!.run({}, { empty: false }, new AbortController().signal);
    expect(nav.expanded.value).toBe(false);
    // The operation itself never toggles: on, and on again, is already so.
    expectOk(await v.api.view.expand({ on: true }));
    expect(await v.api.view.expand({ on: true })).toMatchObject({ ok: true, unchanged: true });
  });

  it('refuses a switch value that is not one, naming the shape', async () => {
    const v = await headlessViewer();
    const r = await v.api.run('canvas.lighting', { shadow: 'on' });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    expect(!r.ok && r.message).toMatch(/takes \{ name, params \}/);
  });
});
