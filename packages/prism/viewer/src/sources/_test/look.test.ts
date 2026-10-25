/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';
import type { SceneNode } from '../../model/scene.js';
import type { PrismViewer } from '../../viewer.js';

// The look's view models (sources/look.ts): the scene defaults' form (each change a `look.set` of the scene), the
// overrides it lists, a node's look in its inspector, and the tree rows' marks.

const Metallic = { metalness: 1, roughness: 0.6, bumpiness: 0 };

/** Every node of the scene tree, depth first (expanded or not). */
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
const node = (v: PrismViewer, id: string) => nodes(v).find((n) => n.id === id)!;

async function withLigand() {
  const v = await headlessViewer({ files: [cbs()] });
  expectOk(await v.api.component.add({ target: 'ligand', label: 'RA', representation: 'ball-and-stick' }));
  const component = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
    (c) => (c.cell.params?.values as any).label === 'RA',
  )!;
  return {
    v,
    structure: component.structure.cell.transform.ref,
    component: component.cell.transform.ref,
    repr: component.representations[0].cell.transform.ref,
  };
}

describe('the scene defaults form', () => {
  it('sets each changed default by its own name (quality, material: not the manager’s names)', async () => {
    const { v } = await withLigand();
    const options = () => v.plugin.managers.structure.component.state.options;
    v.views.look.form.set('quality', 'high');
    v.views.look.form.set('material', Metallic);
    await until(() => options().visualQuality === 'high' && options().materialStyle.metalness === 1);
    // changed elsewhere: the form follows
    expectOk(await v.api.look.set({ ignoreLight: true }));
    await until(() => v.views.look.form.values.value.ignoreLight === true);
  });

  it('lists the representations that keep their own, and Apply to all takes them back', async () => {
    const { v, component } = await withLigand();
    expect(v.views.look.overrides.value).toEqual([]);
    expectOk(await v.api.look.set({ node: component, material: Metallic }));
    await until(() => v.views.look.overrides.value.length === 1);
    expect(v.views.look.overrides.value[0]).toMatchObject({ key: 'material', label: 'Material' });
    expect(v.views.look.overrides.value[0].names[0]).toMatch(/Ball & Stick of RA/);

    v.views.look.applyToAll('material');
    await until(() => v.views.look.overrides.value.length === 0);
  });
});

describe('a node’s look and the rows’ marks', () => {
  it('marks the representation with its own look, and the rows above it', async () => {
    const { v, structure, component, repr } = await withLigand();
    expect(node(v, repr).look?.value).toBeUndefined();
    expectOk(await v.api.look.set({ node: repr, material: Metallic, ignoreLight: true }));
    await until(() => node(v, repr).look?.value !== undefined);
    expect(node(v, repr).look!.value).toEqual({
      kind: 'own',
      text: 'Its own ignore light and material, not the scene default',
    });
    expect(node(v, component).look!.value).toMatchObject({ kind: 'inside' });
    expect(node(v, structure).look!.value).toMatchObject({ kind: 'inside' });

    // the default becomes the override's value: no longer an override
    expectOk(await v.api.look.set({ material: Metallic, ignoreLight: true }));
    await until(() => node(v, repr).look!.value === undefined);
    expect(node(v, structure).look!.value).toBeUndefined();
  });

  it('a component’s inspector has its look: a change sets its representations, Use default resets', async () => {
    const { v, component, repr } = await withLigand();
    const inspection = v.views.scene.inspect(node(v, component));
    const look = inspection.look!;
    expect(look).toBeDefined();
    expect(look.own.value).toEqual([]);

    look.form.set('material', Metallic);
    await until(
      () => (v.plugin.state.data.cells.get(repr)!.transform.params as any).type.params.material.metalness === 1,
    );
    look.refresh();
    expect(look.own.value.map((o) => o.key)).toEqual(['material']);

    look.reset('material');
    await until(
      () => (v.plugin.state.data.cells.get(repr)!.transform.params as any).type.params.material.metalness === 0,
    );
    look.refresh();
    expect(look.own.value).toEqual([]);
  });
});
