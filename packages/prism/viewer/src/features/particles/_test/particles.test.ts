/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { StateSelection } from '@molstar/core/state';
import { describe, expect, it } from 'vitest';
import { expectOk, headlessViewer, until } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { type ParticleListInfo, particleUnitCellOf } from '../operations.js';
import { particlesFile } from '../../../api/_test/volume-fixtures.js';

// `particles` (api.md §2) on a headless viewer, over a small CryoET NDJSON file made in memory.

/** A fixture of `data/examples/particles` (real files, README.md there) */
const particleFixture = (name: string) =>
  new File(
    [fs.readFileSync(path.resolve(import.meta.dirname, '../../../../../../../data/examples/particles', name))],
    name,
  );

const infos = async (v: PrismViewer) => {
  const r = await v.api.particles.list();
  expectOk(r);
  return r.value as ParticleListInfo[];
};

const reps = (v: PrismViewer) =>
  v.plugin.managers.particles.hierarchy.current.lists.flatMap((l) =>
    l.representations.map((r) => (r.cell.transform.params as any).type.name as string),
  );

describe('particles', () => {
  it('has nothing to act on where no particles are loaded', async () => {
    const v = await headlessViewer();
    const none = await v.api.particles.list();
    expectOk(none);
    expect(none.did).toMatch(/No particles/);
    for (const r of [
      await v.api.particles.representation({ type: 'spacefill' }),
      await v.api.particles.visibility({ visible: false }),
      await v.api.particles.remove({}),
    ]) {
      expect(r).toMatchObject({ ok: false, code: 'not-found' });
    }
  });

  it('describes a particle list with its count and representations', async () => {
    const v = await headlessViewer();
    expectOk(await v.api.load.files({ files: [particlesFile()] }));
    const [list] = await infos(v);
    expect(list.count).toBe(5);
    expect(list.representations.map((r) => r.type)).toEqual(['spacefill']);
    const r = await v.api.particles.list();
    expectOk(r);
    expect(r.did).toMatch(/5 particles; spacefill/);
  });

  it('adds a representation or updates the one of that type, and leaves it when nothing changes', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [particlesFile()] });

    expectOk(await v.api.particles.representation({ type: 'orientation' }));
    expect(reps(v).sort()).toEqual(['orientation', 'spacefill']);
    const same = await v.api.particles.representation({ type: 'orientation' });
    expectOk(same);
    expect(same.unchanged).toBe(true);

    expectOk(await v.api.particles.representation({ type: 'spacefill', params: { sizeFactor: 2 } }));
    const spacefill = v.plugin.managers.particles.hierarchy.current.lists[0].representations.find(
      (r) => (r.cell.transform.params as any).type.name === 'spacefill',
    )!;
    expect((spacefill.cell.transform.params as any).type.params.sizeFactor).toBe(2);
    expect(
      ((await v.api.particles.representation({ type: 'spacefill', params: { sizeFactor: 2 } })) as any).unchanged,
    ).toBe(true);

    const bad = await v.api.particles.representation({ type: 'blob' });
    expect(bad).toMatchObject({ ok: false });
    expect((bad as any).hint).toMatch(/spacefill/);
    expect(await v.api.particles.representation({ type: 'spacefill', params: { nope: 1 } })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it('updates a representation from the inspector form', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [particlesFile()] });
    const rep = v.plugin.managers.particles.hierarchy.current.lists[0].representations[0];
    const current = rep.cell.transform.params as any;
    const params = { type: { name: 'spacefill', params: { ...current.type.params, alpha: 0.4 } } };
    expectOk(await v.api.particles.update({ node: rep.cell.transform.ref, params }));
    expect((rep.cell.transform.params as any).type.params.alpha).toBe(0.4);
    expect(((await v.api.particles.update({ node: rep.cell.transform.ref, params })) as any).unchanged).toBe(true);
    expect(await v.api.particles.update({ node: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('shows, hides and removes the list or a representation, and undoes', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [particlesFile()] });
    const list = () => v.plugin.managers.particles.hierarchy.current.lists[0];

    expectOk(await v.api.particles.visibility({ visible: false }));
    expect(list().cell.state.isHidden).toBe(true);
    expect(((await v.api.particles.visibility({ visible: false })) as any).unchanged).toBe(true);
    expectOk(await v.api.particles.visibility({ visible: true }));
    expectOk(await v.api.particles.visibility({ visible: false, representation: 'spacefill' }));
    expect((await infos(v))[0].representations[0].hidden).toBe(true);

    expectOk(await v.api.particles.remove({ representation: 'spacefill' }));
    expect(reps(v)).toEqual([]);
    expectOk(await v.api.history.undo());
    expect(reps(v)).toEqual(['spacefill']);

    expectOk(await v.api.particles.remove({}));
    expect(await infos(v)).toHaveLength(0);
    expectOk(await v.api.history.undo());
    expect(await infos(v)).toHaveLength(1);
  });

  it('colors the representations of a list by a theme, and says when it already is', async () => {
    const v = await headlessViewer({ files: [particleFixture('relion-subtomo.star')] });
    const theme = () =>
      (v.plugin.managers.particles.hierarchy.current.lists[0].representations[0].cell.transform.params as any)
        .colorTheme.name;
    expectOk(await v.api.particles.color({ theme: 'uniform', params: { value: 0xff0000 } }));
    expect(theme()).toBe('uniform');
    expect(((await v.api.particles.color({ theme: 'uniform', params: { value: 0xff0000 } })) as any).unchanged).toBe(
      true,
    );
    expectOk(await v.api.particles.color({ theme: 'particle-index', representation: 'spacefill' }));
    expect(theme()).toBe('particle-index');
    expectOk(await v.api.history.undo());
    expect(theme()).toBe('uniform');
    expect(await v.api.particles.color({ theme: 'nope' })).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('focuses a list, as a view change', async () => {
    const v = await headlessViewer({ files: [particleFixture('cryoet-ribosomes.ndjson')] });
    const r = await v.api.particles.focus({});
    expectOk(r);
    expect(r.did).toMatch(/Focused/);
  });

  it("shows, styles and hides a list's unit cell; a list without one says so", async () => {
    const v = await headlessViewer({ files: [particleFixture('cryoet-ribosomes.ndjson')] });
    expect(await v.api.particles.unitCell({})).toMatchObject({ ok: false, code: 'blocked' });

    // Simularium trajectories carry their box
    const w = await headlessViewer({ files: [particleFixture('small.simularium')] });
    const list = w.plugin.managers.particles.hierarchy.current.lists[0];
    expect(list.cell.obj?.data.cell).toBeDefined();
    const cell = () => particleUnitCellOf(w.plugin.state.data, list.cell.transform.ref);
    const shown = cell() && !cell()!.state.isHidden;
    if (shown) expect(((await w.api.particles.unitCell({ visible: true })) as any).unchanged).toBe(true);
    else expectOk(await w.api.particles.unitCell({}));
    expect(cell()!.state.isHidden).toBeFalsy();

    expectOk(await w.api.particles.unitCell({ color: '#ff8800', opacity: 0.5 }));
    expect((cell()!.transform.params as any).cellColor).toBe(0xff8800);
    expect((cell()!.transform.params as any).alpha).toBe(0.5);
    expectOk(await w.api.particles.unitCell({ visible: false }));
    expect(cell()!.state.isHidden).toBe(true);
    expect(await w.api.particles.unitCell({ color: 'nope' })).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('places a structure at the particles, as their target, and shows it', async () => {
    const v = await headlessViewer({ files: [particleFixture('cryoet-ribosomes.ndjson')] });
    const crambin = new File(
      [fs.readFileSync(path.resolve(import.meta.dirname, '../../../../../../../data/examples/1crn.cif'))],
      '1crn.cif',
    );
    expect(await v.api.particles.targets({})).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.particles.targets({ file: crambin, target: 'nope' })).toMatchObject({
      ok: false,
      code: 'not-found',
    });
    const r = await v.api.particles.targets({ file: crambin });
    expectOk(r);
    expect(r.did).toMatch(/stand for 1crn.cif/);
    const list = v.plugin.managers.particles.hierarchy.current.lists[0];
    const decorated = v.plugin.state.data.select(
      StateSelection.Generators.byRef(list.cell.transform.ref)
        .subtree()
        .filter((c) => c.transform.transformer.id === 'ms-plugin.particle-list-with-targets'),
    )[0];
    expect((decorated.obj as any).data.targetMapping.get(0)).toBeDefined();
    expect(reps(v).sort()).toEqual(['spacefill', 'target']);
    expectOk(await v.api.history.undo());
    expect(reps(v)).toEqual(['spacefill']);
  });

  it('reads a list again from its file with other params (a tomogram), and says what it does not take', async () => {
    const v = await headlessViewer({ files: [particleFixture('dynamo.tbl')] });
    const count = async () => (await infos(v))[0].count;
    const all = await count();
    const r = await v.api.particles.source({ params: { tomos: ['2'] } });
    expectOk(r);
    expect(await count()).toBe(150);
    expect(r.did).toMatch(/150 particles/);
    expect(((await v.api.particles.source({ params: { tomos: ['2'] } })) as any).unchanged).toBe(true);
    expect(await v.api.particles.source({ params: { nope: 1 } })).toMatchObject({ ok: false, code: 'invalid' });
    expectOk(await v.api.history.undo());
    expect(await count()).toBe(all);
  });

  it("moves through a particle trajectory's frames, as viewing", async () => {
    const v = await headlessViewer({ files: [particleFixture('small.simularium')] });
    const frame = () =>
      (v.plugin.managers.particles.hierarchy.current.lists[0].cell.transform.params as any).frameIndex;
    const steps = v.managers.history.steps.value.length;
    expectOk(await v.api.particles.frame({ index: 3 }));
    expect(frame()).toBe(3);
    expectOk(await v.api.particles.frame({ by: 2 }));
    expect(frame()).toBe(0);
    expect(((await v.api.particles.frame({ index: 0 })) as any).unchanged).toBe(true);
    expect(await v.api.particles.frame({ index: 9 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.particles.frame({})).toMatchObject({ ok: false, code: 'invalid' });
    // viewing: no history step
    expect(v.managers.history.steps.value.length).toBe(steps);

    const w = await headlessViewer({ files: [particleFixture('cryoet-ribosomes.ndjson')] });
    expect(await w.api.particles.frame({ index: 1 })).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('plays the frames with the animation, and stops', async () => {
    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await v.api.load.files({ files: [particleFixture('small.simularium')] }));
      const frame = () =>
        (v.plugin.managers.particles.hierarchy.current.lists[0].cell.transform.params as any).frameIndex;
      expectOk(await v.api.particles.play({ mode: 'loop', seconds: 2 }));
      expect(await v.api.particles.play({})).toMatchObject({ ok: true, unchanged: true });
      const animation = v.plugin.managers.animation;
      await animation.tick(0, true);
      await animation.tick(1000, true);
      await until(() => frame() !== 0);
      expectOk(await v.api.particles.stop());
      expect(await v.api.particles.stop()).toMatchObject({ ok: true, unchanged: true });
    } finally {
      await v.dispose();
    }
  });
});
