/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, recordRuns, until } from '../../_test/harness.js';

// Measurements in the scene tree (parity.md §2, "measurement nodes"): a group after the data, a row each.

describe('the scene tree', () => {
  it('lists measurements in their group, and removes one through measure.remove', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const scene = v.views.scene;
    expect(scene.roots.value.some((n) => n.kind === 'measurements')).toBe(false);
    expectOk(await v.api.run('measure.orientation', { target: 'ligand' }));
    await until(() => scene.roots.value.some((n) => n.kind === 'measurements'));
    const group = scene.roots.value.find((n) => n.kind === 'measurements')!;
    expect(group.children.value).toHaveLength(1);
    const row = group.children.value[0];
    expect(row.kind).toBe('measurement');
    const runs = recordRuns(v);
    row.remove();
    await until(() => runs.some((r) => r.operation === 'measure.remove'));
    await until(() => !v.views.scene.roots.value.some((n) => n.kind === 'measurements'));
  });
});
