/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, recordRuns, until } from '../../_test/harness.js';
import { mapFile } from '../../api/_test/volume-fixtures.js';
import type { SceneNode } from '../../model/scene.js';
import type { PrismViewer } from '../../viewer.js';
import { componentsOf } from '../../testing/scene-groups.js';

// The scene's rows and inspections as the scene changes under them: a representation row removes itself whatever its
// type is now; an inspection follows what it was built from (an assembly, a volume's representations); only the picked
// node's is kept.

const typeOf = (v: PrismViewer, ref: string) =>
  (v.plugin.state.data.cells.get(ref)?.transform.params as { type?: { name: string } } | undefined)?.type?.name;
const polymer = (v: PrismViewer) =>
  componentsOf(v.views.scene.roots.value[0]).find((c) => c.kind === 'component' && c.label.value === 'Polymer')!;
const facts = (v: PrismViewer, n: SceneNode) => new Map(v.views.scene.inspect(n).facts ?? []);

describe('a representation row', () => {
  it('removes itself after its type changed, and only itself beside another of its type', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const component = polymer(v);
    const row = component.children.value[0];
    const quick = v.views.scene.inspect(row).quick!;
    quick.change({ type: 'spacefill' });
    quick.settle();
    await until(() => typeOf(v, row.id) === 'spacefill' && v.api.running.value.length === 0);
    // Another spacefill on the same component (a session's, a state tree's).
    const c = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
      (x) => x.cell.transform.ref === component.id,
    )!;
    await v.plugin.builders.structure.representation.addRepresentation(c.cell, { type: 'spacefill' });
    await until(() => component.children.value.length === 2);
    const other = component.children.value.find((n) => n !== row)!;

    const runs = recordRuns(v);
    row.remove();
    await until(() => runs.some((r) => r.operation === 'representation.remove'));
    expect(runs.find((r) => r.operation === 'representation.remove')?.result?.ok).toBe(true);
    await until(() => component.children.value.length === 1);
    expect(component.children.value[0]).toBe(other);
  });
});

describe('an inspection', () => {
  it("follows a structure's assembly: its facts and its contents", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const structure = v.views.scene.roots.value[0];
    const contents = structure.children.value.find((c) => c.kind === 'contents')!;
    const atoms = facts(v, structure).get('Atoms');
    const entityAtoms = () => v.views.scene.inspect(contents).contents!.entities.reduce((n, e) => n + e.atoms, 0);
    const before = entityAtoms();
    v.views.scene.pick(structure);
    const shown = v.views.scene.inspection!.value;
    expect(shown?.node).toBe(structure);

    expectOk(await v.api.structure.symmetry({ kind: 'mates', radius: 5 }));
    // The picked one built anew, so the inspector shows it.
    expect(v.views.scene.inspection!.value).not.toBe(shown);
    expect(v.views.scene.inspection!.value?.node).toBe(structure);
    expect(facts(v, structure).get('Atoms')).not.toBe(atoms);
    expect(entityAtoms()).toBeGreaterThan(before);
  });

  // (A direct volume needs WebGL, which a headless viewer has not: its transfer function follows the same key.)
  it("follows a volume's representations: the levels of its isosurfaces come and go with them", async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const volume = v.views.scene.roots.value[0];
    v.views.scene.pick(volume);
    expect(v.views.scene.inspection!.value?.levels).toHaveLength(1);

    expectOk(await v.api.volume.remove({ node: volume.children.value[0].id }));
    expect(v.views.scene.inspection!.value?.levels).toBeUndefined();
    expectOk(await v.api.volume.representation({ type: 'isosurface' }));
    expect(v.views.scene.inspection!.value?.node).toBe(volume);
    expect(v.views.scene.inspection!.value?.levels).toHaveLength(1);
  });

  it("follows a volume representation's type: built anew for its new type", async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const rep = v.views.scene.roots.value[0].children.value[0];
    v.views.scene.pick(rep);
    const shown = v.views.scene.inspection!.value;
    const uniform = { colorTheme: { name: 'uniform', params: { value: 0xff0000 } } };
    expectOk(await v.api.volume.update({ node: rep.id, params: uniform }));
    // Its params changed, its type not: the same inspection, its form following.
    expect(v.views.scene.inspection!.value).toBe(shown);

    expectOk(await v.api.volume.update({ node: rep.id, params: { type: { name: 'slice', params: {} } } }));
    expect(typeOf(v, rep.id)).toBe('slice');
    expect(v.views.scene.inspection!.value).not.toBe(shown);
    expect(v.views.scene.inspection!.value?.node).toBe(rep);
  });

  it("is kept for the picked node only: another's is built anew when looked at again", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const structure = v.views.scene.roots.value[0];
    const row = polymer(v).children.value[0];
    v.views.scene.pick(row);
    const quick = v.views.scene.inspection!.value?.quick;
    expect(v.views.scene.inspect(row).quick).toBe(quick);
    v.views.scene.pick(structure);
    expect(v.views.scene.inspect(row).quick).not.toBe(quick);
  });
});
