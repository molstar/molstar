/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// `contents` (api.md §2, README decision 18): entities and chains left out of the structure's components.

/** Atoms the components of the first structure hold, together. */
const componentAtoms = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures[0].components.reduce(
    (n, c) => n + (c.cell.obj?.data.elementCount ?? 0),
    0,
  );

describe('contents', () => {
  it('hides an entity from every component, and shows it again', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = componentAtoms(v);
    const listed = await v.api.contents.get();
    expectOk(listed);
    expect(listed.did).toMatch(/entities 1/);

    expectOk(await v.api.contents.visibility({ entity: '1', visible: false }));
    const hidden = componentAtoms(v);
    expect(hidden).toBeLessThan(before);
    const again = await v.api.contents.visibility({ entity: '1', visible: false });
    expect(again).toMatchObject({ ok: true, unchanged: true });

    expectOk(await v.api.contents.visibility({ entity: '1', visible: true }));
    expect(componentAtoms(v)).toBe(before);
  });

  it('keeps what is hidden out of a preset applied later, and undoes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = componentAtoms(v);
    expectOk(await v.api.contents.visibility({ chain: 'A', visible: false }));
    expectOk(await v.api.structure.preset({ preset: 'atomic-detail' }));
    // Chain A is the whole protein and its ligand: what stays is outside it (waters of other chains, if any).
    expect(componentAtoms(v)).toBeLessThan(before);
    await v.api.history.undo();
    await v.api.history.undo();
    expect(componentAtoms(v)).toBe(before);
  });

  it('keeps what is hidden out of the components a symmetry change rebuilds', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.contents.visibility({ chain: 'A', visible: false }));
    const hidden = componentAtoms(v);
    expectOk(await v.api.structure.symmetry({ kind: 'mates', radius: 5 }));
    const s = v.plugin.managers.structure.hierarchy.current.structures[0];
    // Chain A's copies are left out too: the components hold less than the structure they were made from.
    expect(componentAtoms(v)).toBeLessThan(s.cell.obj!.data.elementCount);
    expect(hidden).toBeLessThan(s.cell.obj!.data.elementCount);
  });

  it('refuses unknown entities and chains with the ones it has', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.contents.visibility({ chain: 'Z', visible: false });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
    if (!r.ok) expect(r.hint).toMatch(/Chains: A/);
    expect(await v.api.contents.visibility({ visible: false })).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('representation rows', () => {
  it('hide and show a representation, and update its params as one step', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const repr = v.plugin.managers.structure.hierarchy.current.structures[0].components[0].representations[0];
    const ref = repr.cell.transform.ref;
    expectOk(await v.api.representation.visibility({ node: ref, visible: false }));
    expect(v.plugin.state.data.cells.get(ref)!.state.isHidden).toBe(true);
    expect(await v.api.representation.visibility({ node: ref, visible: false })).toMatchObject({
      unchanged: true,
    });

    const params = repr.cell.transform.params as { colorTheme: { name: string } };
    expectOk(
      await v.api.representation.update({
        node: ref,
        params: { colorTheme: { name: 'uniform', params: {} } },
      }),
    );
    expect((v.plugin.state.data.cells.get(ref)!.transform.params as typeof params).colorTheme.name).toBe('uniform');
    await v.api.history.undo();
    expect((v.plugin.state.data.cells.get(ref)!.transform.params as typeof params).colorTheme.name).toBe(
      params.colorTheme.name,
    );
  });

  it('keeps a rename and a change of a component made while contents are hidden', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const ligand = () =>
      v.plugin.managers.structure.hierarchy.current.structures[0].components.find((c) =>
        ['Ligand', 'Retinoic'].includes(c.cell.obj?.label ?? ''),
      );
    const atoms = ligand()?.cell.obj?.data.elementCount ?? 0;
    expectOk(await v.api.contents.visibility({ entity: '1', visible: false }));
    const renamed = await v.api.component.rename({ component: 'Ligand', label: 'Retinoic' });
    expectOk(renamed);
    expect(renamed.did).toBe("Renamed 'Ligand' to 'Retinoic'.");
    // Residues of the hidden protein: added to its own params, not shown yet.
    expectOk(await v.api.component.modify({ component: 'Retinoic', target: 'resi 10-20', mode: 'add' }));
    expect(ligand()?.cell.obj?.data.elementCount).toBe(atoms);
    expectOk(await v.api.contents.visibility({ entity: '1', visible: true }));
    expect(ligand()?.cell.obj?.label).toBe('Retinoic');
    expect(ligand()?.cell.obj?.data.elementCount).toBeGreaterThan(atoms);
  });

  it('hides a chain while something is focused (the focus cells are not components)', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.focus.set({ target: 'ligand' }));
    expectOk(await v.api.contents.visibility({ entity: '1', visible: false }));
  });
});
