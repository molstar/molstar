/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { RepresentationClip } from '@molstar/plugin/state/helpers/representation-clip';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { cbs } from '../../../_test/harness.js';
import { mapFile } from '../../../api/_test/volume-fixtures.js';
import { afterEach, describe, expect, it } from 'vitest';
import { type Cut, findCuts } from '../../../sources/clipping/cuts.js';

// The `clipping` area over a headless viewer: cuts written into the representations' params and read back from them.

const triangle = () => new File(['v 0 0 0\nv 10 0 0\nv 0 10 0\nf 1 2 3\n'], 'triangle.obj');

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

async function create(...files: File[]) {
  viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
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
    const clipping = v.api.clipping;
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
    const clipping = v.api.clipping;
    const r = await clipping.add({ shape: 'sphere', target: 'polymer', types: ['molecular-surface'] });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
    if (!r.ok) expect(r.hint).toMatch(/cartoon/);
    const cut = ok<Cut>(await clipping.add({ shape: 'sphere', target: 'polymer', types: ['cartoon'] }));
    expect(cut.targets.map((t) => t.type)).toEqual(['cartoon']);
  });

  it('changes, inverts and removes a cut, each a step undo takes back', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
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
    const cut = ok<Cut>(await v.api.clipping.add({ shape: 'cube', node: 'shape' }));
    expect(cut.targets.map((t) => t.type)).toEqual(['shape']);
  });

  it('turns a cut on and off for one representation, never off the last', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    const polymer = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
      (c) => c.cell.obj?.label === 'Polymer',
    )!;
    const cartoon = polymer.representations[0].cell.transform.ref;
    const ligand = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
      (c) => c.cell.obj?.label === 'Ligand',
    )!.representations[0].cell.transform.ref;
    ok(await clipping.add({ shape: 'plane', target: 'polymer' }));
    const cut = findCuts(v.plugin)[0];
    expect(clipping.candidates(cut).map((t) => t.ref)).toContain(ligand);

    ok(await clipping.apply({ cut: 1, representation: ligand, on: true }));
    expect(findCuts(v.plugin)[0].targets.map((t) => t.ref)).toEqual([cartoon, ligand]);
    expect(await clipping.apply({ cut: 1, representation: ligand, on: true })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    ok(await clipping.apply({ cut: 1, representation: cartoon, on: false }));
    expect(clipOf(v, cartoon).objects).toHaveLength(0);
    // the last one: removing the cut is clipping.remove's
    expect(await clipping.apply({ cut: 1, representation: ligand, on: false })).toMatchObject({
      ok: false,
      code: 'blocked',
    });
  });

  it('refuses a direction of length 0, and keeps a rotation axis unit', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    ok(await clipping.add({ shape: 'cube', target: 'polymer' }));
    expect(await clipping.update({ cut: 1, normal: [0, 0, 0] })).toMatchObject({ ok: false, code: 'invalid' });
    ok(await clipping.update({ cut: 1, rotation: { axis: [0, 0, 2], angle: 30 } }));
    const { axis } = findCuts(v.plugin)[0].object.rotation;
    expect(Math.hypot(axis[0], axis[1], axis[2])).toBeCloseTo(1, 6);
  });

  it('Clip this… cuts what a row draws now, the selection left as it is', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    const polymer = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
      (c) => c.cell.obj?.label === 'Polymer',
    )!;
    const cartoon = polymer.representations[0].cell.transform.ref;
    expect(clipping.canClipThis(cartoon)).toBe(true);
    const done = await clipping.clipThis(v.api, cartoon);
    expect(done && done.cut?.targets.map((t) => t.ref)).toEqual([cartoon]);
    expect(v.plugin.managers.structure.selection.stats.elementCount).toBe(0);
    // a component's row: each of its representations; again: the same cut
    const all = await clipping.clipThis(v.api, polymer.cell.transform.ref);
    expect(all && all.cut?.targets.length).toBe(polymer.representations.length);
    const again = await clipping.clipThis(v.api, polymer.cell.transform.ref);
    expect(again && again.cut?.index).toBe(all && all.cut?.index);
  });

  it('turns a cut on for several representations as one step, and changes its shape in place', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    ok(await clipping.add({ shape: 'plane', target: 'polymer' }));
    const cut = findCuts(v.plugin)[0];
    const others = clipping.candidates(cut).filter((t) => !cut.targets.some((x) => x.ref === t.ref));
    expect(others.length).toBeGreaterThan(1);
    ok(await clipping.apply({ cut: 1, representations: others.map((t) => t.ref), on: true }));
    expect(findCuts(v.plugin)[0].targets).toHaveLength(cut.targets.length + others.length);
    await v.api.history.undo();
    expect(findCuts(v.plugin)[0].targets).toHaveLength(cut.targets.length);

    ok(await clipping.update({ cut: 1, shape: 'sphere' }));
    const sphere = findCuts(v.plugin)[0].object;
    expect(sphere.type).toBe('sphere');
    expect([...sphere.position]).toEqual([...cut.object.position]);
    expect(sphere.scale[0]).toBeGreaterThan(0);
  });

  it('adds a cut like one already there once: Add again gives that one', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    const first = ok<Cut>(await clipping.add({ shape: 'plane', target: 'polymer' }));
    const again = await clipping.add({ shape: 'plane', target: 'polymer' });
    expect(again).toMatchObject({ ok: true, unchanged: true, value: { index: first.index } });
    expect(findCuts(v.plugin)).toHaveLength(1);
    expect(clipOf(v, first.targets[0].ref).objects).toHaveLength(1);
    // another shape is another cut
    ok(await clipping.add({ shape: 'sphere', target: 'polymer' }));
    expect(findCuts(v.plugin)).toHaveLength(2);
  });

  it('cuts a map with nothing selected, and Clip this… on its isosurface makes a cut there', async () => {
    const v = await create(mapFile());
    const clipping = v.api.clipping;
    const iso = findCuts(v.plugin).length;
    expect(iso).toBe(0);
    // nothing selected: everything drawn
    const all = ok<Cut>(await clipping.add({ shape: 'plane', everything: true }));
    expect(all.targets.map((t) => t.type)).toContain('isosurface');
    ok(await clipping.remove({ cut: all.index }));
    // Clip this… on the representation: a cut on it, made now (no atoms to select)
    const ref = all.targets.find((t) => t.type === 'isosurface')!.ref;
    expect(clipping.canClipThis(ref)).toBe(true);
    const done = await clipping.clipThis(v.api, ref);
    expect(done && done.cut?.targets.map((t) => t.ref)).toEqual([ref]);
  });

  it('tags its cuts: an edit keeps the tag, and two cuts alike stay two', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    const first = ok<Cut>(await clipping.add({ shape: 'plane', target: 'polymer' }));
    expect(first.object.tag).toMatch(/^clip-/);
    ok(await clipping.update({ cut: 1, position: [1, 2, 3] }));
    expect(findCuts(v.plugin)[0].object.tag).toBe(first.object.tag);
    const second = ok<Cut>(await clipping.add({ shape: 'sphere', target: 'polymer' }));
    expect(second.object.tag).not.toBe(first.object.tag);
    // the second made a plane where the first is: still two cuts, each removed on its own
    ok(
      await clipping.update({
        cut: 2,
        shape: 'plane',
        position: [1, 2, 3],
        normal: [...findCuts(v.plugin)[0].object.rotation.axis],
      }),
    );
    ok(await clipping.update({ cut: 2, rotation: findCuts(v.plugin)[0].object.rotation }));
    expect(findCuts(v.plugin)).toHaveLength(2);
    ok(await clipping.remove({ cut: 2 }));
    expect(findCuts(v.plugin)).toHaveLength(1);
  });

  it('draws a cut in the scene once, moves the drawing with the cut, and takes it away with it', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    const state = v.plugin.state.data;
    const cut = ok<Cut>(await clipping.add({ shape: 'cube', target: 'polymer' }));
    const others = clipping.candidates(cut).filter((t) => !cut.targets.some((x) => x.ref === t.ref));
    ok(await clipping.apply({ cut: 1, representations: [others[0].ref], on: true }));

    ok(await clipping.show({ cut: 1, on: true }));
    const drawn = findCuts(v.plugin)[0];
    expect(drawn.drawing).toBeDefined();
    const drawing = state.cells.get(drawn.drawing!)!;
    expect(drawing.transform.tags).toContain(`clip-volume-${cut.object.tag}`);
    expect(drawing.transform.parent).toBe(cut.targets[0].ref);
    expect(drawing.obj).toBeDefined();
    // once a cut, and never itself something to cut
    expect(await clipping.show({ cut: 1, on: true })).toMatchObject({ ok: true, unchanged: true });
    expect(clipping.candidates(findCuts(v.plugin)[0]).some((t) => t.label.includes('Cut'))).toBe(false);
    // the cut moved: drawn again where it is (the shape rebuilt from its parent's clip)
    const before = state.cells.get(drawn.drawing!)!.obj;
    ok(await clipping.update({ cut: 1, position: [5, 6, 7] }));
    expect(state.cells.get(findCuts(v.plugin)[0].drawing!)!.obj).not.toBe(before);
    // its options
    ok(await clipping.show({ cut: 1, on: true, options: { wallOpacity: 0.4, edges: false } }));
    expect(state.cells.get(findCuts(v.plugin)[0].drawing!)!.transform.params).toMatchObject({
      wallOpacity: 0.4,
      edges: false,
    });

    // the cut off its first representation: the drawing goes under the other one
    ok(await clipping.apply({ cut: 1, representation: cut.targets[0].ref, on: false }));
    const moved = findCuts(v.plugin)[0];
    expect(state.cells.get(moved.drawing!)!.transform.parent).toBe(others[0].ref);
    // removed: so is its drawing
    ok(await clipping.remove({ cut: 1 }));
    expect([...state.cells.values()].some((c) => c.transform.tags?.includes('clip-volumes'))).toBe(false);
  });

  it('keeps atoms whole, releases them, and keeps them with their cut when an earlier one goes', async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    ok(await clipping.add({ shape: 'plane', target: 'polymer' }));
    ok(await clipping.add({ shape: 'sphere', target: 'polymer' }));
    expect(findCuts(v.plugin)[1].kept).toBeUndefined();
    ok(await clipping.keep({ cut: 2, target: 'A:10-20' }));
    const kept = findCuts(v.plugin)[1].kept;
    expect(kept?.count).toBeGreaterThan(0);
    expect(findCuts(v.plugin)[0].kept).toBeUndefined();
    // the first cut removed: the second (now 1) still keeps them (its bit moved down with it)
    ok(await clipping.remove({ cut: 1 }));
    expect(findCuts(v.plugin)[0].kept).toEqual(kept);
    ok(await clipping.release({ cut: 1 }));
    expect(findCuts(v.plugin)[0].kept).toBeUndefined();
    expect(await clipping.release({ cut: 1 })).toMatchObject({ ok: true, unchanged: true });
  });

  it('adds a cut that keeps its target whole', async () => {
    const v = await create(cbs());
    const cut = ok<Cut>(await v.api.clipping.add({ shape: 'sphere', target: 'A:10-20', keepWhole: true }));
    expect(cut.kept?.count).toBeGreaterThan(0);
  });

  it('finds a cut by its tag: a run after another that renumbered the cuts changes the cut it meant', async () => {
    const v = await create(cbs());
    const { clipping } = v.api;
    ok(await clipping.add({ shape: 'plane', target: 'A:10-20' }));
    const second = ok<Cut>(await clipping.add({ shape: 'sphere', target: 'A:10-20' }));
    const tag = second.object.tag!;
    // removed twice by its tag (a double click): the second finds it gone, the first cut stays
    ok(await clipping.remove({ tag }));
    expect(await clipping.remove({ tag })).toMatchObject({ ok: false, code: 'not-found' });
    expect(ok<Cut[]>(await clipping.list())).toHaveLength(1);
  });

  it('refuses numbers that are not a vector, a size not positive and drawing options it does not have', async () => {
    const v = await create(cbs());
    const { clipping } = v.api;
    ok(await clipping.add({ shape: 'cube', target: 'A:10-20' }));
    expect(await clipping.update({ cut: 1, position: [1, 2] })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await clipping.update({ cut: 1, scale: [0, 1, 1] })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await clipping.update({ cut: 1, rotation: { axis: [0, 1, 0], angle: Number.NaN } })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await clipping.show({ cut: 1, options: { wallOpacity: 'x' } as never })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await clipping.show({ cut: 1, options: { nope: 1 } as never })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it("keeps a box's twist when its direction changes", async () => {
    const v = await create(cbs());
    const { clipping } = v.api;
    ok(await clipping.add({ shape: 'cube', target: 'A:10-20' }));
    // twisted 30° about its up (+Y), then its up turned to +X
    ok(await clipping.update({ cut: 1, rotation: { axis: [0, 1, 0], angle: 30 } }));
    ok(await clipping.update({ cut: 1, normal: [1, 0, 0] }));
    const cut = ok<Cut[]>(await clipping.list())[0];
    const m = Mat4.fromRotation(Mat4(), (cut.object.rotation.angle * Math.PI) / 180, cut.object.rotation.axis);
    const up = Vec3.transformDirection(Vec3(), Vec3.create(0, 1, 0), m);
    const side = Vec3.transformDirection(Vec3(), Vec3.create(1, 0, 0), m);
    [1, 0, 0].forEach((x, i) => expect(up[i]).toBeCloseTo(x, 5));
    // the twist: its local X stays where the 30° turn put it, carried by the turn of the up (about Z)
    const twisted = Vec3.transformDirection(
      Vec3(),
      Vec3.create(1, 0, 0),
      Mat4.fromRotation(Mat4(), (30 * Math.PI) / 180, Vec3.create(0, 1, 0)),
    );
    const turn = Mat4.fromRotation(Mat4(), -Math.PI / 2, Vec3.create(0, 0, 1));
    Vec3.transformDirection(twisted, twisted, turn);
    [0, 1, 2].forEach((i) => expect(side[i]).toBeCloseTo(twisted[i], 5));
  });

  it('keeps nothing whole for a node selector, and says so', async () => {
    const v = await create(cbs());
    const surface = [...v.plugin.state.data.cells.values()].find((c) => RepresentationClip.isClippable(c))!.transform
      .ref;
    const r = await v.api.clipping.add({ shape: 'sphere', node: `#${surface}`, keepWhole: true });
    expect(r).toMatchObject({ ok: true });
    expect((r as { warnings?: string[] }).warnings?.[0]).toMatch(/Nothing is kept whole/);
    expect(ok<Cut[]>(await v.api.clipping.list())[0].kept).toBeUndefined();
  });

  it("draws a cut's new shape when its shape changes", async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    const state = v.plugin.state.data;
    ok(await clipping.add({ shape: 'plane', target: 'polymer' }));
    ok(await clipping.show({ cut: 1, on: true }));
    const vertices = () => {
      const drawing = findCuts(v.plugin)[0].drawing!;
      const [shape] = [...state.tree.children.get(drawing).values()];
      const repr = (state.cells.get(shape)!.obj!.data as { repr: { renderObjects: { values: any }[] } }).repr;
      return repr.renderObjects[0].values.drawCount.ref.value as number;
    };
    const plane = vertices();
    ok(await clipping.update({ cut: 1, shape: 'cube' }));
    expect(vertices()).not.toBe(plane);
  });

  it("redraws a drawing's geometry in place when its options change (its extent, relative, as drawn)", async () => {
    const v = await create(cbs());
    const clipping = v.api.clipping;
    const state = v.plugin.state.data;
    ok(await clipping.add({ shape: 'plane', target: 'polymer' }));
    ok(await clipping.show({ cut: 1, on: true }));
    const radius = () => {
      const [shape] = [...state.tree.children.get(findCuts(v.plugin)[0].drawing!).values()];
      const repr = (state.cells.get(shape)!.obj!.data as { repr: { renderObjects: { values: any }[] } }).repr;
      return repr.renderObjects[0].values.boundingSphere.ref.value.radius as number;
    };
    const auto = radius();
    ok(await clipping.show({ cut: 1, on: true, options: { extent: 2 } }));
    expect(radius()).toBeGreaterThan(1.9 * auto);
    // the render object's buffers are the mesh rebuilt in place: a smaller one is drawn smaller
    ok(await clipping.show({ cut: 1, on: true, options: { extent: 0.5 } }));
    expect(radius()).toBeLessThan(0.6 * auto);
  });
});

describe('clipping.settings', () => {
  it('makes every cut representation, and the cuts added later, cut by instance; solid surface is the scene default', async () => {
    const v = await create(cbs(), triangle());
    const clipping = v.api.clipping;
    ok(await v.api.select.set({ target: 'ligand' }));
    const cut = ok<Cut>(await clipping.add({ shape: 'plane' }));
    const ligand = cut.targets[0].ref;
    expect(clipOf(v, ligand).variant).toBe('pixel');

    const r = await clipping.settings({ variant: 'instance' });
    expect(r.ok && r.did).toBe('Cuts leave out whole copies.');
    expect(clipOf(v, ligand).variant).toBe('instance');
    expect(v.plugin.managers.structure.component.state.options.clipObjects.variant).toBe('instance');
    expect(clipping.cutSettings.value).toEqual({ variant: 'instance', solidSurface: true });
    expect(await clipping.settings({ variant: 'instance' })).toMatchObject({ ok: true, unchanged: true });

    // a cut added now, on the mesh (not the component manager's), cuts so too
    const mesh = ok<Cut>(await clipping.add({ shape: 'plane', node: 'shape' }));
    expect(clipOf(v, mesh.targets[0].ref).variant).toBe('instance');

    ok(await clipping.settings({ solidSurface: false }));
    expect(v.plugin.managers.structure.component.state.options.solidSurface).toBe(false);
    expect(clipping.cutSettings.value.solidSurface).toBe(false);
  });
});
