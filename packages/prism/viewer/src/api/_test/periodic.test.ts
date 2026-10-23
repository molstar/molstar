/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import type { Model } from '@molstar/model/model/structure';
import { expectOk, headlessViewer, until } from '../../_test/harness.js';
import { createTrajectoryControls } from '../../ui/trajectory.js';
import type { PrismViewer } from '../../viewer.js';
import { ensemble } from './trajectory-fixtures.js';

// `trajectory.periodic` (api.md §2 `trajectory`) on a headless viewer, over the repository's split peptide: a
// peptide and three waters in a dodecahedral box, split across a face in its first frame.

const gro = () =>
  new File(
    [
      fs.readFileSync(
        path.resolve(import.meta.dirname, '../../../../../../data/examples/trajectory/periodic/split-peptide.gro'),
      ),
    ],
    'split-peptide.gro',
  );

const model = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures[0].model!.cell.obj!.data as Model;
const params = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures[0].model!.cell.transform.params as { periodic?: string };

/** The CA–CA distance of residues 2 and 3: across the face in the first frame as read. */
function caca(m: Model) {
  const { x, y, z } = m.atomicConformation;
  return Math.hypot(x[6] - x[11], y[6] - y[11], z[6] - z[11]);
}

describe('trajectory.periodic', () => {
  it('shows the molecules whole by default, and as the file has them when asked, as a history step', async () => {
    const v = await headlessViewer({ files: [gro()] });
    expect(caca(model(v))).toBeLessThan(4);
    const steps = v.managers.history.steps.value.length;

    const off = await v.api.trajectory.periodic({ mode: 'off' });
    expectOk(off);
    expect(off.did).toMatch(/as in the file/);
    expect(caca(model(v))).toBeGreaterThan(10);
    expect(v.managers.history.steps.value.length).toBe(steps + 1);

    expect(await v.api.trajectory.periodic({ mode: 'off' })).toMatchObject({ ok: true, unchanged: true });

    expectOk(await v.api.history.undo());
    expect(params(v).periodic ?? 'whole').toBe('whole');
    expect(caca(model(v))).toBeLessThan(4);
  });

  it('centers the box on the polymers, and keeps the center given before', async () => {
    const v = await headlessViewer({ files: [gro()] });
    const r = await v.api.trajectory.periodic({ mode: 'box', center: 'polymers' });
    expectOk(r);
    expect(r.did).toMatch(/centered on the polymers/);
    expect(caca(model(v))).toBeLessThan(4);
    // nojump, then box again without a center: still centered on the polymers
    expectOk(await v.api.trajectory.periodic({ mode: 'nojump' }));
    expect(await v.api.trajectory.periodic({ mode: 'box' })).toMatchObject({ ok: true });
    expect(
      (v.plugin.managers.structure.hierarchy.current.structures[0].model!.cell.transform.params as any).periodicCenter
        .name,
    ).toBe('selection');
  });

  it('is greyed out, with the reason, when no structure is in a simulation box', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const r = await v.api.trajectory.periodic({ mode: 'box' });
    expect(r).toMatchObject({ ok: false });
    expect(r.ok ? '' : r.message).toMatch(/simulation box/);
  });

  it("drives the trajectory bar's box menu", async () => {
    const v = await headlessViewer({ files: [gro()] });
    const bar = createTrajectoryControls(v.plugin, v.api);
    const periodic = bar.model.periodic!;
    const off = periodic.mode.listen?.(() => {});
    await until(() => periodic.mode.value === 'whole');
    periodic.set('box-centered');
    await until(() => periodic.mode.value === 'box-centered');
    expect(params(v).periodic).toBe('box');
    periodic.set('nojump');
    await until(() => periodic.mode.value === 'nojump');
    off?.();
    bar.dispose();

    // an ensemble has no box: no menu
    const plain = await headlessViewer({ files: [ensemble()] });
    const none = createTrajectoryControls(plain.plugin, plain.api);
    expect(none.model.periodic!.mode.value).toBeUndefined();
  });
});
