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
import { cbs, exampleFile, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// parity.md §2, phase 6: the unit cell (`structure.unit-cell`) and the Conformation Transform (`structure.transform`,
// `structure.transform-reset`) on a headless viewer.

const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const structure = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures[0];
const unitCell = (v: PrismViewer) => structure(v).model?.unitcell?.cell;
const center = (v: PrismViewer) => {
  const s = structure(v);
  const leaf = s.transform?.cell.obj?.data ?? s.cell.obj!.data;
  return Vec3.clone(leaf.boundary.sphere.center);
};

describe('structure.unit-cell', () => {
  it('shows and hides the unit cell of a crystal structure, each a step, and undoes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = steps(v);
    const r = await v.api.structure.unitCell({ visible: true });
    expectOk(r);
    expect(r.did).toMatch(/shows its unit cell \(P 21 21 21\)/);
    expect(!!unitCell(v)?.state.isHidden).toBe(false);
    expect(steps(v)).toBe(before + 1);
    expect(await v.api.structure.unitCell({ visible: true })).toMatchObject({ ok: true, unchanged: true });
    expect(steps(v)).toBe(before + 1);

    expectOk(await v.api.structure.unitCell({ visible: false }));
    expect(unitCell(v)?.state.isHidden).toBe(true);
    await v.api.history.undo();
    expect(unitCell(v)).toBeDefined();
    expect(!!unitCell(v)?.state.isHidden).toBe(false);
  });

  it('sets its color and opacity, and keeps them when only the visibility changes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.structure.unitCell({ color: '#ff0000', opacity: 0.5 }));
    const params = () => unitCell(v)?.transform.params as { cellColor: Color; alpha: number };
    expect(params()).toMatchObject({ cellColor: Color(0xff0000), alpha: 0.5 });
    expectOk(await v.api.structure.unitCell({ visible: false }));
    expect(params()).toMatchObject({ cellColor: Color(0xff0000), alpha: 0.5 });
    expect(await v.api.structure.unitCell({ color: 'nonsense' })).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('makes the unit cell when the model has none', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const existing = unitCell(v);
    if (existing) await v.plugin.state.data.build().delete(existing.transform.ref).commit();
    expect(unitCell(v)).toBeUndefined();
    expectOk(await v.api.structure.unitCell({ visible: true }));
    expect(!!unitCell(v)?.state.isHidden).toBe(false);
  });

  it('is blocked without a crystal cell, and asks for a structure that is there', async () => {
    const v = await headlessViewer({ files: [exampleFile('adrenalin.mol2')] });
    const r = await v.api.structure.unitCell({ visible: true });
    expect(r).toMatchObject({ ok: false, code: 'blocked' });
    expect(await v.api.structure.unitCell({ structure: 'nothing' })).toMatchObject({ ok: false, code: 'not-found' });
  });
});

describe('structure.transform', () => {
  it('moves a structure by a translation as one step, says already so, and undoes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const c0 = center(v);
    const before = steps(v);
    const r = await v.api.structure.transform({ translation: Vec3.create(10, 0, 0) });
    expectOk(r);
    expect(center(v)[0]).toBeCloseTo(c0[0] + 10, 3);
    expect(steps(v)).toBe(before + 1);
    expect(await v.api.structure.transform({ translation: Vec3.create(10, 0, 0) })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    await v.api.history.undo();
    expect(center(v)[0]).toBeCloseTo(c0[0], 3);
  });

  it('rotates about the structure’s center, which stays put', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const c0 = center(v);
    expectOk(await v.api.structure.transform({ axis: Vec3.create(0, 0, 1), angle: 90 }));
    const c1 = center(v);
    expect(Vec3.distance(c0, c1)).toBeLessThan(1e-3);
    const t = structure(v).transform?.cell.transform.params?.transform;
    expect(t?.name).toBe('components');
  });

  it('turns a structure in place when composing a rotation on a moved one', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const c0 = center(v);
    expectOk(await v.api.structure.transform({ translation: Vec3.create(50, 0, 0) }));
    expectOk(await v.api.structure.transform({ axis: Vec3.create(0, 0, 1), angle: 90, compose: true }));
    // Its center stays where the move put it: the turn is about the structure, not about where it was loaded.
    const c1 = center(v);
    expect(Vec3.distance(c1, Vec3.create(c0[0] + 50, c0[1], c0[2]))).toBeLessThan(1e-2);
  });

  it('takes a matrix, composes on top, and refuses a matrix with parts', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const c0 = center(v);
    const m = Mat4.fromTranslation(Mat4(), Vec3.create(0, 5, 0));
    expectOk(await v.api.structure.transform({ matrix: m }));
    expect(center(v)[1]).toBeCloseTo(c0[1] + 5, 3);
    expectOk(await v.api.structure.transform({ translation: Vec3.create(0, 5, 0), compose: true }));
    expect(center(v)[1]).toBeCloseTo(c0[1] + 10, 3);
    const both = await v.api.structure.transform({ matrix: m, angle: 10 });
    expect(both).toMatchObject({ ok: false, code: 'invalid' });
    const zero = await v.api.structure.transform({ axis: Vec3.zero(), angle: 10 });
    expect(zero).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('resets to where it was loaded, keeping the components', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const c0 = center(v);
    const components = structure(v).components.length;
    expect(await v.api.structure.transformReset()).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.structure.transform({ translation: Vec3.create(3, 3, 3) }));
    expectOk(await v.api.structure.transformReset());
    expect(Vec3.distance(center(v), c0)).toBeLessThan(1e-3);
    expect(structure(v).components.length).toBe(components);
    expect(await v.api.structure.transformReset()).toMatchObject({ ok: true, unchanged: true });
  });
});
