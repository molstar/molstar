/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { type Cut, findCuts } from '../cuts.js';
import { ClippingExtension } from '../prism-viewer.js';

// The `clipping` area over a headless viewer: cuts written into the representations' params and read back from them.

const cbs = () =>
  new File([readFileSync(new URL('../../../../data/examples/1cbs_updated.cif', import.meta.url))], '1cbs.cif');
const triangle = () => new File(['v 0 0 0\nv 10 0 0\nv 0 10 0\nf 1 2 3\n'], 'triangle.obj');

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

async function create(...files: File[]) {
  const base = DefaultPrismViewerSpec();
  viewer = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ClippingExtension] },
  });
  for (const f of files) expect((await viewer.api.load.files({ files: [f] })).ok).toBe(true);
  return viewer;
}

const ok = <T>(r: { ok: boolean; message?: string; value?: T }) => {
  if (!r.ok) throw new Error(r.message);
  return r.value as T;
};
const clipOf = (v: PrismViewer, ref: string) =>
  RepresentationClip.get(v.plugin.state.data.cells.get(ref)!.transform.params)!;

describe('the clipping area', () => {
  it('cuts the representations drawing the selection, and lists the cut', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping!;
    expect(ok<Cut[]>(await clipping.list())).toEqual([]);
    expect(await clipping.add({ shape: 'plane' })).toMatchObject({ ok: false, code: 'not-found' });

    ok(await v.api.select.set({ target: 'ligand' }));
    const cut = ok<Cut>(await clipping.add({ shape: 'plane' }));
    expect(cut.index).toBe(1);
    expect(cut.targets.map((t) => t.label)).toEqual(['Ball & Stick of Ligand']);
    expect(clipOf(v, cut.targets[0].ref).objects).toHaveLength(1);
    // the rest is not cut
    const all = findCuts(v.plugin);
    expect(all).toHaveLength(1);
    const listed = await clipping.list();
    expect(listed.ok && listed.did).toMatch(/Cut 1: plane at .* · Ball & Stick of Ligand/);
  });

  it('narrows to representation types, and says what it is drawn as when none matches', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping!;
    const r = await clipping.add({ shape: 'sphere', target: 'polymer', types: ['molecular-surface'] });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
    if (!r.ok) expect(r.hint).toMatch(/cartoon/);
    const cut = ok<Cut>(await clipping.add({ shape: 'sphere', target: 'polymer', types: ['cartoon'] }));
    expect(cut.targets.map((t) => t.type)).toEqual(['cartoon']);
  });

  it('changes, inverts and removes a cut, each a step undo takes back', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping!;
    const { targets } = ok<Cut>(await clipping.add({ shape: 'plane', target: 'polymer' }));
    const ref = targets[0].ref;

    ok(await clipping.update({ cut: 1, normal: [0, -1, 0], position: [1, 2, 3] }));
    const o = clipOf(v, ref).objects[0];
    expect([...o.position]).toEqual([1, 2, 3]);
    // a normal pointing straight down turns half a turn (not none: H005)
    expect(o.rotation.angle).toBeCloseTo(180);
    expect((await clipping.update({ cut: 1, position: [1, 2, 3] })) as any).toMatchObject({ unchanged: true });

    ok(await clipping.update({ cut: 1, invert: true }));
    expect(clipOf(v, ref).objects[0].invert).toBe(true);
    ok(await v.api.history.undo());
    expect(clipOf(v, ref).objects[0].invert).toBe(false);

    ok(await clipping.remove({ cut: 1 }));
    expect(clipOf(v, ref).objects).toHaveLength(0);
    expect(await clipping.remove({ cut: 1 })).toMatchObject({ ok: false, code: 'not-found' });
    ok(await v.api.history.undo());
    expect(findCuts(v.plugin)).toHaveLength(1);
    ok(await clipping.clear());
    expect(findCuts(v.plugin)).toHaveLength(0);
  });

  it('lists a clip set by code, and not the canvas-wide one (the component options)', async () => {
    const v = await create(cbs());
    const structure = v.plugin.managers.structure.hierarchy.current.structures[0];
    const cell = structure.components[0].representations[0].cell;
    const object = {
      type: 'sphere' as const,
      invert: false,
      position: Vec3.create(10, 20, 30),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
      scale: Vec3.create(8, 8, 8),
      transform: Mat4.identity(),
    };
    await v.plugin
      .build()
      .to(cell)
      .update((old: any) => {
        old.type.params.clip = { variant: 'pixel', objects: [object] };
      })
      .commit();
    const cuts = findCuts(v.plugin);
    expect(cuts).toHaveLength(1);
    expect(cuts[0].object.type).toBe('sphere');

    const manager = v.plugin.managers.structure.component;
    await manager.setOptions({
      ...manager.state.options,
      clipObjects: { variant: 'pixel', objects: [{ ...object, type: 'plane' }] },
    });
    expect(findCuts(v.plugin).map((c) => c.object.type)).toEqual(['sphere']);
  });

  it('cuts a shape by node selector', async () => {
    const v = await create(triangle());
    const cut = ok<Cut>(await v.api.clipping!.add({ shape: 'cube', node: 'shape' }));
    expect(cut.targets.map((t) => t.type)).toEqual(['shape']);
  });
});
