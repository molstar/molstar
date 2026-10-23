/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { Model } from '@molstar/model/model/structure';
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

  it("reads a MolViewSpec structure's custom.molstar_periodic", async () => {
    const text = fs.readFileSync(
      path.resolve(import.meta.dirname, '../../../../../../data/examples/trajectory/periodic/split-peptide.gro'),
      'utf8',
    );
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: string | URL | Request) => {
        const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
        return url.includes('files.test/split-peptide.gro')
          ? new Response(text, { status: 200 })
          : new Response('not found', { status: 404 });
      }),
    );
    const mvs = JSON.stringify({
      metadata: { version: '1', timestamp: '2026-10-23T00:00:00Z' },
      root: {
        kind: 'root',
        children: [
          {
            kind: 'download',
            params: { url: 'https://files.test/split-peptide.gro' },
            children: [
              {
                kind: 'parse',
                params: { format: 'gro' },
                children: [
                  {
                    kind: 'structure',
                    params: { type: 'model' },
                    custom: { molstar_periodic: { mode: 'box', center: 'polymers' } },
                    children: [
                      {
                        kind: 'component',
                        params: { selector: 'all' },
                        children: [{ kind: 'representation', params: { type: 'ball_and_stick' } }],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    });
    const v = await headlessViewer({});
    expectOk(await v.api.load.mvs!({ data: mvs, keepCamera: true }));
    expect(params(v).periodic).toBe('box');
    expect(caca(model(v))).toBeLessThan(4);
    // the peptide's mean at the box's center
    const m = model(v);
    const box = Model.PeriodicBox.get(m)!;
    const { x, y, z } = m.atomicConformation;
    let cx = 0,
      cy = 0,
      cz = 0;
    for (let i = 0; i < 20; i++) {
      cx += x[i] / 20;
      cy += y[i] / 20;
      cz += z[i] / 20;
    }
    const b = box.basis;
    expect(cx).toBeCloseTo((b[0] + b[3] + b[6]) / 2, 3);
    expect(cy).toBeCloseTo((b[1] + b[4] + b[7]) / 2, 3);
    expect(cz).toBeCloseTo((b[2] + b[5] + b[8]) / 2, 3);
  });

  it("shows the periodic images from the trajectory bar's box menu", async () => {
    const v = await headlessViewer({ files: [gro()] });
    const bar = createTrajectoryControls(v.plugin, v.api);
    const images = bar.model.periodic!.images!;
    const off = images.shown.listen?.(() => {});
    const units = () =>
      (v.plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data as any).units.length;
    const alone = units();
    await until(() => images.shown.value === false);
    images.set(true);
    await until(() => images.shown.value === true);
    // the box and its 26 neighbours
    expect(units()).toBe(27 * alone);
    images.set(false);
    await until(() => images.shown.value === false);
    expect(units()).toBe(alone);
    off?.();
    bar.dispose();
  });

  it("takes a DCD's boxes with their origins from a NAMD extended system file", async () => {
    const data = (name: string) =>
      fs.readFileSync(path.resolve(import.meta.dirname, `../../../../../../data/examples/trajectory/${name}`));
    // a 40 Å cube centered on (5, -3, 2), for every step
    const xst = '#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n0 40 0 0 0 40 0 0 0 40 5 -3 2\n';
    const served: Record<string, BodyInit> = {
      'files.test/protein.pdb': data('protein.pdb'),
      'files.test/protein.dcd': data('protein.dcd'),
      'files.test/protein.xst': xst,
    };
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: string | URL | Request) => {
        const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
        const hit = Object.keys(served).find((k) => url.includes(k));
        return hit ? new Response(served[hit], { status: 200 }) : new Response('not found', { status: 404 });
      }),
    );
    const v = await headlessViewer({});
    expectOk(
      await v.api.load.trajectory!({
        model: 'https://files.test/protein.pdb',
        modelFormat: 'pdb',
        coordinates: 'https://files.test/protein.dcd',
        coordinatesFormat: 'dcd',
        boxes: 'https://files.test/protein.xst',
      }),
    );
    const m = model(v);
    const box = Model.PeriodicBox.get(m)!;
    expect(Array.from(box.origin)).toEqual([5 - 20, -3 - 20, 2 - 20]);
    // in the box: every molecule's atoms around the box, from its corner
    expectOk(await v.api.trajectory.periodic({ mode: 'box' }));
    const { x } = model(v).atomicConformation;
    const mean = Array.from(x).reduce((a, b) => a + b, 0) / x.length;
    expect(mean).toBeGreaterThan(-15);
    expect(mean).toBeLessThan(25);
  });

  it('gives a dropped NAMD box file to the trajectory loaded', async () => {
    const data = (name: string) =>
      fs.readFileSync(path.resolve(import.meta.dirname, `../../../../../../data/examples/trajectory/${name}`));
    const served: Record<string, BodyInit> = {
      'files.test/protein.pdb': data('protein.pdb'),
      'files.test/protein.dcd': data('protein.dcd'),
    };
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: string | URL | Request) => {
        const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
        const hit = Object.keys(served).find((k) => url.includes(k));
        return hit ? new Response(served[hit], { status: 200 }) : new Response('not found', { status: 404 });
      }),
    );
    const v = await headlessViewer({});
    expectOk(
      await v.api.load.trajectory!({
        model: 'https://files.test/protein.pdb',
        modelFormat: 'pdb',
        coordinates: 'https://files.test/protein.dcd',
        coordinatesFormat: 'dcd',
      }),
    );
    expectOk(
      await v.api.load.files({ files: [new File(['0 40 0 0 0 40 0 0 0 40 5 -3 2\n'], 'protein.xsc')] as never }),
    );
    expect(Array.from(Model.PeriodicBox.get(model(v))!.origin)).toEqual([5 - 20, -3 - 20, 2 - 20]);
  });

  it("reads a MolViewSpec coordinates node's custom.molstar_boxes, its address relative to the view's", async () => {
    const data = (name: string) =>
      fs.readFileSync(path.resolve(import.meta.dirname, `../../../../../../data/examples/trajectory/${name}`));
    const download = (url: string, format: string, child: object) => ({
      kind: 'download',
      params: { url },
      children: [{ kind: 'parse', params: { format }, children: [child] }],
    });
    const mvs = JSON.stringify({
      metadata: { version: '1', timestamp: '2026-10-23T00:00:00Z' },
      root: {
        kind: 'root',
        children: [
          download('protein.dcd', 'dcd', {
            kind: 'coordinates',
            ref: 'frames',
            params: {},
            custom: { molstar_boxes: { url: 'protein.xsc' } },
          }),
          download('protein.pdb', 'pdb', {
            kind: 'structure',
            params: { type: 'model', coordinates_ref: 'frames' },
            children: [
              {
                kind: 'component',
                params: { selector: 'all' },
                children: [{ kind: 'representation', params: { type: 'ball_and_stick' } }],
              },
            ],
          }),
        ],
      },
    });
    const served: Record<string, BodyInit> = {
      'files.test/views/view.mvsj': mvs,
      'files.test/views/protein.pdb': data('protein.pdb'),
      'files.test/views/protein.dcd': data('protein.dcd'),
      'files.test/views/protein.xsc': '0 40 0 0 0 40 0 0 0 40 5 -3 2\n',
    };
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: string | URL | Request) => {
        const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
        const hit = Object.keys(served).find((k) => url.includes(k));
        return hit ? new Response(served[hit], { status: 200 }) : new Response('not found', { status: 404 });
      }),
    );
    const v = await headlessViewer({});
    expectOk(await v.api.load.mvs!({ url: 'https://files.test/views/view.mvsj', keepCamera: true }));
    // MVS builds model, trajectory with coordinates, model, structure: the hierarchy names the model the structure is
    // built from (the frames'), so the trajectory and the periodic box are the viewer's
    const shown = (v.plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data as { model: Model })
      .model;
    expect(model(v)).toBe(shown);
    expect(Array.from(Model.PeriodicBox.get(model(v))!.origin)).toEqual([5 - 20, -3 - 20, 2 - 20]);
    const frames = await v.api.trajectory.get();
    expectOk(frames);
    expect(frames.value).toMatchObject({ structures: [{ count: 10, unit: 'frames' }] });
    expectOk(await v.api.trajectory.periodic({ mode: 'box' }));
  });
});
