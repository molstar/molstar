/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, recordRuns, until } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// parity.md §2, the structure's source and crystal symmetry in the structure inspector and the structure row's menu:
// fields and entries over `structure.model`, `.assembly`, `.symmetry` and `.preset`, acting as the person.

const structure = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures[0];
const typeName = (v: PrismViewer) => (structure(v).cell.transform.params as { type: { name: string } }).type.name;
const typeParams = (v: PrismViewer) =>
  (structure(v).cell.transform.params as { type: { params: Record<string, unknown> } }).type.params;

/** Two models of three atoms each, as PDB text. */
function twoModels(): File {
  const atom = (n: number, name: string, x: number) =>
    `ATOM  ${String(n).padStart(5)}  ${name.padEnd(3)} ALA A   1    ${x.toFixed(3).padStart(8)}   0.000   0.000  1.00  0.00           ${name[0]}`;
  const model = (i: number, dx: number) => [
    `MODEL     ${String(i).padStart(4)}`,
    atom(1, 'N', dx),
    atom(2, 'CA', dx + 1.5),
    atom(3, 'C', dx + 3),
    'ENDMDL',
  ];
  return new File([[...model(1, 0), ...model(2, 1), 'END'].join('\n')], 'two.pdb');
}

describe('the structure inspector: the source', () => {
  it("offers a crystal's types: the model, its assemblies, symmetry mates and unit cells; one model, no model field", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const form = v.views.scene.inspect(v.views.scene.roots.value[0]).form!;
    const type = form.params.type as { options: readonly (readonly [string, string])[]; isHidden?: boolean };
    expect(type.options.map(([id]) => id)).toEqual(expect.arrayContaining(['model', 'assembly:1', 'mates', 'cells']));
    expect(type.isHidden).toBeFalsy();
    expect(form.params.model.isHidden).toBe(true);
    expect(form.layout().essentials).toEqual(expect.arrayContaining(['type', 'transform']));
  });

  it('sets the type through the operations, and follows the structure', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const form = v.views.scene.inspect(v.views.scene.roots.value[0]).form!;
    const values = () => form.values.value as { type: string; radius: number; cells: number };
    const runs = recordRuns(v);

    form.set('type', 'mates');
    await until(() => runs.some((r) => r.operation === 'structure.symmetry' && !!r.result?.ok), 10000);
    expect(typeName(v)).toBe('symmetry-mates');
    expect(runs.find((r) => r.operation === 'structure.symmetry')?.values).toMatchObject({ kind: 'mates' });
    expect(values().type).toBe('mates');

    form.set('radius', 8);
    await until(() => typeParams(v).radius === 8, 10000);

    form.set('type', 'cells');
    await until(() => typeName(v) === 'symmetry', 10000);

    form.set('type', 'assembly:1');
    await until(() => typeName(v) === 'assembly', 10000);
    expect(typeParams(v).id).toBe('1');

    const before = runs.length;
    form.set('type', 'model');
    await until(() => runs.length > before && !!runs.at(-1)?.result, 10000);
    expect(typeName(v)).toBe('model');
    expect(runs.at(-1)?.operation).toBe('structure.assembly');

    // A change made elsewhere (a command, the agent): the form follows.
    expectOk(await v.api.structure.symmetry({ kind: 'cells', cells: 2 }));
    await until(() => values().type === 'cells' && values().cells === 2, 10000);
  });

  it('shows another model of a trajectory', async () => {
    const v = await headlessViewer({ files: [twoModels()] });
    const form = v.views.scene.inspect(v.views.scene.roots.value[0]).form!;
    expect(form.params.model.isHidden).toBeFalsy();
    expect(form.values.value.model).toBe(1);
    form.set('model', 2);
    await until(() => (structure(v).model?.cell.transform.params as { modelIndex: number }).modelIndex === 1, 10000);
    expect(form.values.value.model).toBe(2);
    // A crystal-less file: no symmetry among its types.
    const type = form.params.type as { options: readonly (readonly [string, string])[] };
    expect(type.options.map(([id]) => id)).not.toContain('mates');
  });

  it('puts the values back when a change is refused', async () => {
    const v = await headlessViewer({ files: [twoModels()] });
    const form = v.views.scene.inspect(v.views.scene.roots.value[0]).form!;
    // No crystal: `structure.symmetry` refuses mates, and the field returns to the model.
    form.set('type', 'mates');
    await until(() => form.values.value.type === 'model', 10000);
  });
});

describe('the structure row menu', () => {
  it('has the presets, the assemblies, crystal symmetry and the transform reset', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const menu = v.views.scene.roots.value[0].menu!;
    expect(menu.map((m) => m.id)).toEqual(
      expect.arrayContaining(['preset', 'assembly', 'symmetry', 'reset-transform']),
    );
    const assembly = menu.find((m) => m.id === 'assembly')!;
    expect(assembly.items!.map((i) => i.id)).toEqual(['model', 'assembly:1']);

    const runs = recordRuns(v);
    menu.find((m) => m.id === 'symmetry')!.items!.find((i) => i.id === 'mates')!.run!();
    await until(() => runs.some((r) => r.operation === 'structure.symmetry' && !!r.result?.ok), 10000);
    expect(typeName(v)).toBe('symmetry-mates');

    const preset = menu.find((m) => m.id === 'preset')!;
    expect(preset.items!.length).toBeGreaterThan(2);
    preset.items!.find((i) => i.id === 'preset-structure-representation-polymer-cartoon')!.run!();
    await until(() => runs.some((r) => r.operation === 'structure.preset' && !!r.result), 10000);
  });

  it('offers streaming density while it is off', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.structure.assembly({ id: '' }));
    const menu = v.views.scene.roots.value[0].menu!;
    expect(menu.some((m) => m.id === 'stream-density')).toBe(!!v.api.describe('volume.streaming'));
  });
});
