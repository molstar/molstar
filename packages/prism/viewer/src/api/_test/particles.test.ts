/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';
import type { ParticleListInfo } from '../particles.js';
import { particlesFile } from './volume-fixtures.js';

// `particles` (api.md §2) on a headless viewer, over a small CryoET NDJSON file made in memory.

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
});
