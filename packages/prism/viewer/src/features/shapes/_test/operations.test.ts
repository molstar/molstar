/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import type { SceneNode } from '../../../model/scene.js';
import { MeshFiles } from '../../../testing/meshes.js';
import type { PrismViewer } from '../../../viewer.js';
import type { ShapeInfo } from '../operations.js';

// The `shape` area and the Structure panel's Shapes group over a headless viewer: meshes made in code (an OBJ torus, a
// PLY octahedron), listed, shown and hidden, updated, given a look, removed.

async function withMeshes(files: File[] = [MeshFiles.torus(), MeshFiles.octahedron()]) {
  const v = await headlessViewer({ files: [cbs()] });
  expectOk(await v.api.load.files({ files }));
  return v;
}

const shapes = async (v: PrismViewer) => {
  const r = await v.api.shape.list();
  expectOk(r);
  return r.value as ShapeInfo[];
};
const params = (v: PrismViewer, ref: string) => v.plugin.state.data.cells.get(ref)!.transform.params as any;

/** Every node of the scene tree, depth first. */
function nodes(v: PrismViewer): SceneNode[] {
  const out: SceneNode[] = [];
  const walk = (list: readonly SceneNode[]) => {
    for (const n of list) {
      out.push(n);
      walk(n.children.value);
    }
  };
  walk(v.views.scene.roots.value);
  return out;
}

describe('the shape area', () => {
  it('lists the meshes loaded, by their files, with their representations; not the unit cell', async () => {
    const v = await withMeshes();
    const list = await shapes(v);
    expect(list.map((s) => s.label)).toEqual(['torus.obj', 'octahedron.ply']);
    expect(list.every((s) => s.representations.length === 1)).toBe(true);
    const r = await v.api.shape.list();
    expect(r.ok && r.did).toMatch(/torus\.obj \(.*\): Mesh/);
  });

  it('hides and shows a shape, one step; already so', async () => {
    const v = await withMeshes();
    const [torus] = await shapes(v);
    expectOk(await v.api.shape.visibility({ node: torus.ref, visible: false }));
    expect(v.plugin.state.data.cells.get(torus.ref)!.state.isHidden).toBe(true);
    expect(await v.api.shape.visibility({ node: torus.ref, visible: false })).toMatchObject({ unchanged: true });
    // empty: every shape
    expectOk(await v.api.shape.visibility({ visible: true }));
    expect(v.plugin.state.data.cells.get(torus.ref)!.state.isHidden).toBe(false);
  });

  it('updates a representation (the inspector form), and sets its look by node', async () => {
    const v = await withMeshes();
    const [torus] = await shapes(v);
    const mesh = torus.representations[0].ref;
    expectOk(await v.api.shape.update({ node: mesh, params: { alpha: 0.5 } }));
    expect(params(v, mesh).alpha).toBe(0.5);

    // the look: by the shape's node, written to its representations (their params, not a type's)
    const r = await v.api.look.set({
      node: torus.ref,
      ignoreLight: true,
      material: { metalness: 1, roughness: 0.6, bumpiness: 0 },
    });
    expectOk(r);
    expect(r.did).toMatch(/Set ignore light and material of Mesh/);
    expect(params(v, mesh)).toMatchObject({ ignoreLight: true, material: { metalness: 1 } });
    // hydrogens are a structure's
    expect(await v.api.look.set({ node: torus.ref, hydrogens: 'hide-all' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it('removes a shape with its data, and asks which one when several are loaded', async () => {
    const v = await withMeshes();
    expect(await v.api.shape.remove({})).toMatchObject({ ok: false, code: 'ambiguous' });
    const [torus] = await shapes(v);
    expectOk(await v.api.shape.remove({ node: torus.ref }));
    expect((await shapes(v)).map((s) => s.label)).toEqual(['octahedron.ply']);
    // its file's data went too
    let files = 0;
    v.plugin.state.data.cells.forEach((c) => {
      if (c.obj?.label === 'torus.obj') files++;
    });
    expect(files).toBe(0);
  });

  it('says when no shape is loaded', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(await v.api.shape.list()).toMatchObject({ ok: true, did: 'No shape is loaded.' });
    expect(await v.api.shape.visibility({ visible: false })).toMatchObject({ ok: false, code: 'not-found' });
  });
});

describe('the Shapes group', () => {
  it('shows each mesh under Shapes, its representation under it, and follows a removal', async () => {
    const v = await withMeshes();
    await until(() => nodes(v).some((n) => n.kind === 'shapes'));
    const group = nodes(v).find((n) => n.kind === 'shapes')!;
    expect(group.description.value).toBe('2 shapes');
    expect(group.children.value.map((n) => [n.kind, n.label.value])).toEqual([
      ['shape', 'torus.obj'],
      ['shape', 'octahedron.ply'],
    ]);
    const torus = group.children.value[0];
    expect(torus.children.value.map((n) => n.kind)).toEqual(['shape-representation']);

    // its eye, through shape.visibility
    torus.toggleVisibility();
    await until(() => torus.isHidden.value);

    // its inspector: a look (no hydrogens), its representation's form with its color at hand
    const look = v.views.scene.inspect(torus).look!;
    expect(Object.keys(look.form.params)).not.toContain('hydrogens');
    expect(Object.keys(look.form.params)).toContain('material');
    const repr = v.views.scene.inspect(torus.children.value[0]);
    expect(repr.form?.params.coloring?.isEssential).toBe(true);

    expectOk(await v.api.shape.remove({ node: torus.id }));
    await until(() => nodes(v).find((n) => n.kind === 'shapes')?.children.value.length === 1);
    expectOk(await v.api.shape.remove({}));
    await until(() => !nodes(v).some((n) => n.kind === 'shapes'));
  });
});
