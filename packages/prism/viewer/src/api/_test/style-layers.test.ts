/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection } from '@molstar/core/state';
import { OverpaintStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';
import type { StyleLayerInfo } from '../../sources/style-layers.js';
import { createStyleLayersModel } from '../../ui/style-layers.js';
import type { PrismViewer } from '../../viewer.js';

// The applied style layers (`sources/style-layers.ts`) and the operations over them: `representation.layers`,
// `layer-update`, `layer-remove`, `reset` with every layer; and the Selection tab's model of them (rows that keep
// their identity through an edit). Headless, 1CBS.

const overpaints = (v: PrismViewer) =>
  v.plugin.state.data.select(StateSelection.Generators.ofTransformer(OverpaintStructureRepresentation3DFromBundle));

async function styled() {
  const v = await headlessViewer({ files: [cbs()] });
  expectOk(await v.api.representation.color({ target: 'A:10-17', color: '#ff0000' }));
  expectOk(await v.api.representation.color({ target: 'ligand', color: 'blue' }));
  expectOk(await v.api.representation.transparency({ target: 'ligand', alpha: 0.4 }));
  expectOk(await v.api.representation.emissive({ target: 'A:30-40', strength: 0.6 }));
  return v;
}

const layers = async (v: PrismViewer): Promise<readonly StyleLayerInfo[]> => {
  const r = await v.api.representation.layers();
  expectOk(r);
  return r.value as readonly StyleLayerInfo[];
};

const steps = async (v: PrismViewer) => {
  const r = await v.api.history.list();
  expectOk(r);
  return (r.value as unknown[]).length;
};

describe('representation.layers', () => {
  it('lists each value of a kind on a structure, with where it is, by kind', async () => {
    const v = await styled();
    const list = await layers(v);
    expect(list.map((l) => [l.kind, l.value, l.where])).toEqual([
      ['color', '#ff0000', 'A:10-17'],
      ['color', '#0000ff', 'A:200'],
      ['transparency', '0.4', 'A:200'],
      ['emissive', '0.6', 'A:30-40'],
    ]);
    expect(list[0]).toMatchObject({ color: '#ff0000', structure: { label: '1CBS' } });
    expect(list[2].amount).toBe(0.4);
    expect(list[0].line).toMatch(/^8 residues \(\d+ atoms\) in chain A$/);
    expect(list[0].id).toBe(`color:${list[0].structure.ref}:#ff0000`);
    // A read: no history step.
    const before = await steps(v);
    await layers(v);
    expect(await steps(v)).toBe(before);
  });

  it('is empty without layers, and the source follows the state', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(await layers(v)).toEqual([]);
    const seen: number[] = [];
    const sub = v.sources.layers.subscribe((l) => seen.push(l.length));
    try {
      expectOk(await v.api.representation.color({ target: 'ligand', color: 'red' }));
      await until(() => seen[seen.length - 1] === 1);
      expectOk(await v.api.history.undo());
      await until(() => seen[seen.length - 1] === 0);
    } finally {
      sub.unsubscribe();
    }
  });
});

describe('representation.layer-update and layer-remove', () => {
  it('changes a color where it is, one step, and gives the new id', async () => {
    const v = await styled();
    const [red] = await layers(v);
    const before = await steps(v);
    const r = await v.api.representation.layerUpdate({ layer: red.id, value: 'lime' });
    expectOk(r);
    expect(r.value).toEqual({ id: `color:${red.structure.ref}:#00ff00` });
    const after = await layers(v);
    expect(after.find((l) => l.id === red.id)).toBeUndefined();
    expect(after.find((l) => l.id === (r.value as { id: string }).id)?.where).toBe('A:10-17');
    expect(await steps(v)).toBe(before + 1);
    // The same value again: already so.
    expect(
      await v.api.representation.layerUpdate({ layer: (r.value as { id: string }).id, value: '#00ff00' }),
    ).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.history.undo());
    expect((await layers(v))[0].id).toBe(red.id);
  });

  it('changes a transparency, and refuses a bad value or an unknown layer', async () => {
    const v = await styled();
    const t = (await layers(v)).find((l) => l.kind === 'transparency')!;
    expectOk(await v.api.representation.layerUpdate({ layer: t.id, value: '0.8' }));
    expect((await layers(v)).find((l) => l.kind === 'transparency')?.amount).toBe(0.8);
    const now = (await layers(v)).find((l) => l.kind === 'transparency')!;
    expect(await v.api.representation.layerUpdate({ layer: now.id, value: '2' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await v.api.representation.layerUpdate({ layer: now.id, value: 'red' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await v.api.representation.layerRemove({ layer: 'color:nope:#000000' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it('removes one layer and keeps the rest of its kind; the last one takes its cells with it', async () => {
    const v = await styled();
    const [red, blue] = await layers(v);
    const cells = overpaints(v).length;
    expectOk(await v.api.representation.layerRemove({ layer: blue.id }));
    expect((await layers(v)).map((l) => l.id)).toEqual(expect.not.arrayContaining([blue.id]));
    expect((await layers(v)).some((l) => l.id === red.id)).toBe(true);
    expectOk(await v.api.representation.layerRemove({ layer: red.id }));
    expect(overpaints(v).length).toBeLessThan(cells);
    expect(overpaints(v)).toHaveLength(0);
    expectOk(await v.api.history.undo());
    expect((await layers(v)).some((l) => l.id === red.id)).toBe(true);
  });

  it('are not available without layers', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expect(v.api.describe('representation.layer-remove')?.available).toMatch(/No style layers/);
  });
});

describe('representation.reset with every layer', () => {
  it('clears every kind in one step', async () => {
    const v = await styled();
    const before = await steps(v);
    expectOk(await v.api.representation.reset({ layer: 'all' }));
    expect(await layers(v)).toEqual([]);
    expect(await steps(v)).toBe(before + 1);
    expect(await v.api.representation.reset({ layer: 'all' })).toMatchObject({ ok: true, unchanged: true });
  });
});

describe('the style layers model', () => {
  it('keeps a row through an edit of its value, removes and clears through the operations', async () => {
    const v = await styled();
    const failures: string[] = [];
    const host = createStyleLayersModel(v.api, v.sources.layers, (r) => {
      if (!r.ok) failures.push(r.message);
    });
    const { model } = host;
    const seen: number[] = [];
    const sub = model.items.subscribe((l) => seen.push(l.length));
    try {
      await until(() => seen[seen.length - 1] === 4);
      const row = model.items.value.find((i) => i.view.value.value.endsWith('% transparent'))!;
      expect(row.editable).toBe(true);
      expect(row.view.value).toMatchObject({ where: 'A:200', value: '40% transparent', icon: 'eye-off' });
      expect(model.items.value[0].view.value).toMatchObject({ swatch: '#ff0000', value: '' });

      const form = row.editor()!;
      form.set('value', 0.7);
      await until(() => row.view.value.value === '70% transparent');
      // The same row (and its open form): the list's row for the new value is it.
      expect(model.items.value.includes(row)).toBe(true);
      expect((await layers(v)).find((l) => l.kind === 'transparency')?.amount).toBe(0.7);

      // A second edit names the layer by the id the first gave.
      form.set('value', 0.2);
      await until(() => row.view.value.value === '20% transparent');
      expect(model.items.value.includes(row)).toBe(true);

      row.remove();
      await until(() => seen[seen.length - 1] === 3);
      model.clear();
      await until(() => seen[seen.length - 1] === 0);
      expect(failures).toEqual([]);
    } finally {
      sub.unsubscribe();
      host.dispose();
    }
  });
});
