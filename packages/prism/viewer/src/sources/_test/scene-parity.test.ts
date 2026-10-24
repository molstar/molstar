/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, recordRuns, until } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';
import { otherOf } from '../../testing/scene-groups.js';

// parity.md §2 in the scene: the unit cell's row under its structure, and the structure inspector's
// Transform (the Conformation Transform), both acting through the API.

const unitCellRow = (v: PrismViewer) => otherOf(v.views.scene.roots.value[0]).find((c) => c.kind === 'unit-cell');
const structure = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures[0];
const params = (v: PrismViewer) => structure(v).model?.unitcell?.cell.transform.params as { cellColor: Color };

describe('the unit cell row', () => {
  it("is in a crystal structure's Other group, hidden; its eye shows it as a step", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const row = unitCellRow(v)!;
    expect(row).toBeDefined();
    expect(row.label.value).toBe('Unit cell');
    expect(row.description.value).toBe('P 21 21 21');
    expect(row.isHidden.value).toBe(true);
    expect(row.removable).toBe(false);
    // First in the structure's Other group (the volume streaming row, when the entry has density, follows it).
    expect(otherOf(v.views.scene.roots.value[0])[0]).toBe(row);

    const runs = recordRuns(v);
    row.toggleVisibility();
    await until(() => runs.some((r) => r.operation === 'structure.unit-cell' && !!r.result));
    await until(() => !row.isHidden.value);
    expect(v.managers.history.steps.value[0].label).toMatch(/unit cell/);
    // Turned on, it is picked: its controls are in the inspector at once.
    expect(v.views.scene.picked.value).toBe(row);
  });

  it('inspects the cell: space group, size and angles; its form sets the color', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const row = unitCellRow(v)!;
    const inspection = v.views.scene.inspect(row);
    expect(inspection.facts?.find(([k]) => k === 'Space group')?.[1]).toBe('P 21 21 21');
    expect(inspection.facts?.find(([k]) => k === 'Cell')?.[1]).toMatch(/ × .* Å$/);
    inspection.form!.set('color', Color(0x00ff00));
    await until(() => params(v).cellColor === Color(0x00ff00));
    // Still hidden: the form styles it, the eye shows it.
    expect(structure(v).model?.unitcell?.cell.state.isHidden).toBe(true);
    // The form follows Undo.
    await v.api.history.undo();
    await until(() => (inspection.form!.values.value as { color: Color }).color !== Color(0x00ff00));
  });
});

describe('the structure inspector: Transform', () => {
  it('moves the structure from its form, and shows a superposition-like matrix as its parts', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const node = v.views.scene.roots.value[0];
    const inspection = v.views.scene.inspect(node);
    const form = inspection.form!;
    const values = () =>
      (form.values.value as { transform: { translation: Vec3; axis: Vec3; angle: number } }).transform;
    expect(values()).toMatchObject({ angle: 0 });
    form.set('transform', { ...values(), translation: Vec3.create(5, 0, 0) });
    await until(() => !!structure(v).transform);

    // A rotation about the z axis by 90° and a shift, as a matrix: the form reads it back as parts.
    const m = Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3)),
      Mat4.fromRotation(Mat4(), Math.PI / 2, Vec3.create(0, 0, 1)),
    );
    expectOk(await v.api.structure.transform({ matrix: m }));
    await until(() => Math.abs(values().angle - 90) < 1e-3);
    expect(Array.from(values().axis)).toEqual([0, 0, 1]);

    // Its parts are the same transform.
    expect(await v.api.structure.transform({ ...values(), center: 'centroid' })).toMatchObject({
      ok: true,
      unchanged: true,
    });
  });

  it('resets from the row menu', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.structure.transform({ translation: Vec3.create(1, 1, 1) }));
    const node = v.views.scene.roots.value[0];
    const reset = node.menu?.find((m) => m.id === 'reset-transform');
    const runs = recordRuns(v);
    reset!.run!();
    await until(() => runs.some((r) => r.operation === 'structure.transform-reset' && !!r.result));
    expect(await v.api.structure.transformReset()).toMatchObject({ ok: true, unchanged: true });
  });
});
