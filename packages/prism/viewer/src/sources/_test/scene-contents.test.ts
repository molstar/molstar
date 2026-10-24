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
import { componentsOf } from '../../testing/scene-groups.js';

// The contents rows (README decision 18): an entity's or a chain's eye hides it and shows it again, the row says which,
// and hovering it highlights its atoms.

describe('the contents rows', () => {
  it('show what is hidden, toggle back, and highlight on hover', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const contents = v.views.scene.roots.value[0].children.value.find((c) => c.kind === 'contents')!;
    contents.load?.();
    const entity = () =>
      contents.children.value.find((c) => c.kind === 'entity' && /RETINOIC/i.test(c.label.value)) as SceneNode;
    expect(entity().isHidden.value).toBe(false);

    entity().toggleVisibility();
    await until(() => entity().isHidden.value);
    expect(entity().children.value.every((c) => c.isHidden.value)).toBe(true);
    const components = () => v.plugin.managers.structure.hierarchy.current.structures[0].components.length;
    expect(components()).toBe(2);

    entity().toggleVisibility();
    await until(() => !entity().isHidden.value);
    await until(() => components() === 3);

    entity().highlight(true);
    await until(() => v.sources.pointer.hoverLines.value.length > 0);
    entity().highlight(false);
    await until(() => v.sources.pointer.hoverLines.value.length === 0);
  });

  it('hides a chain of its entity only (a ligand and waters of the same author chain stay); emptied rows stay', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const structure = v.views.scene.roots.value[0];
    const contents = structure.children.value.find((c) => c.kind === 'contents')!;
    contents.load?.();
    const polymer = contents.children.value.find(
      (c) => c.kind === 'entity' && /polymer/.test(c.description.value ?? ''),
    )!;
    const chain = polymer.children.value[0] as SceneNode;
    const components = () => v.plugin.managers.structure.hierarchy.current.structures[0].components;
    const before = components().map((c) => c.cell.obj?.label);
    expect(before).toEqual(['Polymer', 'Ligand', 'Water']);

    chain.toggleVisibility();
    await until(() => chain.isHidden.value);
    // The ligand and the waters (author chain A too) are still there; the polymer is empty.
    expect(components().map((c) => c.cell.obj?.label)).toEqual(['Ligand', 'Water']);
    const rows = () => componentsOf(structure).filter((c) => c.kind === 'component');
    expect(rows().map((r) => r.label.value)).toEqual(['Polymer', 'Ligand', 'Water']);
    expect(rows()[0].isHidden.value).toBe(true);
    expect(rows()[0].description.value).toBe('hidden in Entities');
    // Its row kept, but it shows or hides nothing meanwhile: no eye, and the Components group goes around it.
    expect(rows()[0].hasEye).toBe(false);
    const group = structure.children.value.find((c) => c.kind === 'components')!;
    const steps = v.managers.history.steps.value.length;
    group.toggleVisibility();
    await until(() => rows()[1].isHidden.value && rows()[2].isHidden.value && group.isHidden.value);
    group.toggleVisibility();
    await until(() => !rows()[1].isHidden.value && !rows()[2].isHidden.value && !group.isHidden.value);
    expect(v.managers.history.steps.value).toHaveLength(steps + 2);

    chain.toggleVisibility();
    await until(() => components().length === 3);
    expect(rows()[0].isHidden.value).toBe(false);
    expect(rows()[0].hasEye).toBe(true);
  });

  it('gives an emptied component row made without a row before it its actions once it is shown again', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const structure = () => v.views.scene.roots.value[0];
    const contents = () => structure().children.value.find((c) => c.kind === 'contents')!;
    contents().load?.();
    const polymer = contents().children.value.find(
      (c) => c.kind === 'entity' && /polymer/.test(c.description.value ?? ''),
    )!;
    const chain = polymer.children.value[0] as SceneNode;
    chain.toggleVisibility();
    await until(() => chain.isHidden.value);

    // The structure removed and brought back by Undo: its rows are made anew, the polymer's as an emptied one.
    expectOk(await v.api.structure.remove({ structure: structure().id }));
    await until(() => v.views.scene.roots.value[0]?.kind !== 'structure');
    expectOk(await v.api.history.undo());
    await until(() => structure()?.kind === 'structure');
    const rows = () => componentsOf(structure()).filter((c) => c.kind === 'component');
    await until(() => rows()[0]?.description.value === 'hidden in Entities');

    contents().load?.();
    const again = contents().children.value.find(
      (c) => c.kind === 'entity' && /polymer/.test(c.description.value ?? ''),
    )!;
    (again.children.value[0] as SceneNode).toggleVisibility();
    await until(() => rows()[0]?.description.value !== 'hidden in Entities');
    const row = rows()[0];
    expect(row.removable).toBe(true);
    expect(row.menu?.length).toBeGreaterThan(0);
    row.toggleVisibility();
    await until(() => row.isHidden.value);
  });
});
