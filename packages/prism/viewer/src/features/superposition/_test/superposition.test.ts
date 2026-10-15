/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { TransformStructureConformation } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, recordRuns, until } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { SuperpositionFeature } from '../index.js';
import { createSuperpositionModel } from '../model.js';

// The `superposition` feature (operations.md §2, `superpose`): each operation superposes structures onto the first and is
// one history step. The fixtures are copies of one file, one of them moved by a known transform, so the answer is known:
// the transform undoes the move and the RMSD is zero. No canvas here (the camera is not touched).

/** The rotation and translation the second copy is moved by. */
const Move = (() => {
  const m = Mat4.fromRotation(Mat4(), 0.7, Vec3.normalize(Vec3(), Vec3.create(1, 2, 3)));
  Mat4.setTranslation(m, Vec3.create(12, -7, 5));
  return m;
})();

const structures = (plugin: PluginContext) => plugin.managers.structure.hierarchy.current.structures;

/** The matrix a structure's conformation transform holds now (identity without one). */
function matrixOf(s: StructureRef): Mat4 {
  const p = s.transform?.cell.transform.params?.transform as { params: { data: Mat4 } } | undefined;
  return p ? Mat4.clone(p.params.data) : Mat4.identity();
}

async function move(v: PrismViewer, s: StructureRef, matrix: Mat4) {
  const b = v.plugin.state.data
    .build()
    .to(s.cell)
    .insert(TransformStructureConformation, {
      transform: { name: 'matrix', params: { data: matrix, transpose: false } },
    });
  await v.plugin.runTask(v.plugin.state.data.updateTree(b));
}

/** Two copies of `file`, the second moved: `[reference, mobile]`. */
async function twoCopies(file: () => File = cbs, extra: { files?: File[] } = {}) {
  const v = await headlessViewer({ extensions: [SuperpositionFeature] });
  const a = extra.files?.[0] ?? file();
  expectOk(await v.api.load.files({ files: [a] }));
  expectOk(await v.api.load.files({ files: [extra.files?.[1] ?? file()] }, { copy: true }));
  const [ref, mobile] = structures(v.plugin);
  expect(structures(v.plugin)).toHaveLength(2);
  await move(v, mobile, Move);
  return { v, ref: ref.cell.transform.ref, mobile: mobile.cell.transform.ref, api: v.api.superpose! };
}

const mobileOf = (v: PrismViewer, ref: string) => structures(v.plugin).find((s) => s.cell.transform.ref === ref)!;
const steps = async (v: PrismViewer) => {
  const r = await v.api.history.list();
  expectOk(r);
  return (r.value as unknown[]).length;
};

/** Where the superposition puts the mobile copy: the copies have the same coordinates, so no transform at all. */
const Undone = Mat4.identity();

describe('superpose.chains', () => {
  it('superposes a moved copy onto the reference and undoes the move', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    const before = await steps(v);
    const r = await api.chains({
      targets: [
        { target: 'polymer', structure: ref },
        { target: 'polymer', structure: mobile },
      ],
    });
    expectOk(r);
    expect(r.changed.nodes).toEqual([mobile]);
    expect(r.value?.results[0]).toMatchObject({ moved: true, structure: { ref: mobile } });
    expect(r.value!.results[0].rmsd).toBeLessThan(0.05);
    expect(r.did).toMatch(/RMSD 0\.0/);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
    expect(await steps(v)).toBe(before + 1);
  });

  it('is already so when the transforms would not change, without a history step', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    const targets = [
      { target: 'polymer', structure: ref },
      { target: 'polymer', structure: mobile },
    ];
    expectOk(await api.chains({ targets }));
    const before = await steps(v);
    const again = await api.chains({ targets });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(await steps(v)).toBe(before);
  });

  it('puts the transform back with Undo', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    expectOk(
      await api.chains({
        targets: [
          { target: 'polymer', structure: ref },
          { target: 'polymer', structure: mobile },
        ],
      }),
    );
    expectOk(await v.api.history.undo());
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Move, 1e-6)).toBe(true);
  });

  it('positions bare targets by the structures at their place, and works without the sequence alignment', async () => {
    const { v, mobile, api } = await twoCopies();
    const r = await api.chains({ targets: ['polymer', 'polymer'], alignSequences: false, traceOnly: false });
    expectOk(r);
    expect(r.value?.results[0].structure.ref).toBe(mobile);
    expect(r.value!.results[0].rmsd).toBeLessThan(0.05);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
  });

  it('uses the selection of each structure when no targets are given', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    const none = await api.chains({});
    expect(none).toMatchObject({ ok: false, code: 'invalid' });
    expectOk(await v.api.select.set({ target: 'polymer', structure: ref }));
    expectOk(await v.api.select.set({ target: 'polymer', structure: mobile, mode: 'add' }));
    const r = await api.chains({});
    expectOk(r);
    expect(r.value!.results[0].rmsd).toBeLessThan(0.05);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
  });

  it('moves onto a reference that was itself moved', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    const shift = Mat4.fromTranslation(Mat4(), Vec3.create(3, 0, 0));
    await move(v, mobileOf(v, ref), shift);
    expectOk(
      await api.chains({
        targets: [
          { target: 'polymer', structure: ref },
          { target: 'polymer', structure: mobile },
        ],
      }),
    );
    // Both copies now sit where the reference is shown: the reference's shift leads.
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), shift, 0.02)).toBe(true);
  });

  it('fails when a target selects nothing, is not a chain, or two share a structure', async () => {
    const { ref, mobile, api } = await twoCopies();
    const empty = await api.chains({
      targets: [
        { target: 'polymer', structure: ref },
        { target: 'resn XYZ', structure: mobile },
      ],
    });
    expect(empty.ok).toBe(false);
    expect(empty).toMatchObject({ ok: false, hint: expect.any(String) });
    const ligand = await api.chains({
      targets: [
        { target: 'polymer', structure: ref },
        { target: 'ligand', structure: mobile },
      ],
    });
    expect(ligand).toMatchObject({ ok: false, code: 'invalid' });
    const twice = await api.chains({
      targets: [
        { target: 'polymer', structure: ref },
        { target: 'polymer', structure: ref },
      ],
    });
    expect(twice).toMatchObject({ ok: false, code: 'invalid', message: expect.stringContaining('both') });
    const few = await api.chains({ targets: ['polymer'] });
    expect(few).toMatchObject({ ok: false, code: 'invalid' });
    const many = await api.chains({ targets: ['polymer', 'polymer', 'polymer'] });
    expect(many).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('fails invalid for too few atoms', async () => {
    const { ref, mobile, api } = await twoCopies();
    const r = await api.chains({
      targets: [
        { target: { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 2 }, structure: ref },
        { target: { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 2 }, structure: mobile },
      ],
      alignSequences: false,
    });
    expect(r).toMatchObject({ ok: false, code: 'invalid', message: expect.stringContaining('too few') });
  });

  it('is blocked with one structure', async () => {
    const v = await headlessViewer({ extensions: [SuperpositionFeature], files: [cbs()] });
    const r = await v.api.superpose!.chains({ targets: ['polymer', 'polymer'] });
    expect(r).toMatchObject({ ok: false, code: 'blocked' });
    expect(v.api.describe('superpose.chains')?.available).toMatch(/two or more/);
  });
});

describe('superpose.tm-align', () => {
  it('superposes the copies with TM-score 1', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    const r = await api.tmAlign({
      targets: [
        { target: 'polymer', structure: ref },
        { target: 'polymer', structure: mobile },
      ],
    });
    expectOk(r);
    const res = r.value!.results[0];
    expect(res.tmScore).toBeGreaterThan(0.99);
    expect(res.rmsd).toBeLessThan(0.05);
    expect(r.did).toMatch(/TM-score 1\.0/);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.05)).toBe(true);
  });

  it('is already so the second time', async () => {
    const { ref, mobile, api } = await twoCopies();
    const targets = [
      { target: 'polymer', structure: ref },
      { target: 'polymer', structure: mobile },
    ];
    expectOk(await api.tmAlign({ targets }));
    expect(await api.tmAlign({ targets })).toMatchObject({ ok: true, unchanged: true });
  });
});

describe('superpose.atoms', () => {
  const CA = (from: number, to: number) => ({ label_atom_id: 'CA', beg_label_seq_id: from, end_label_seq_id: to });

  it('pairs the atoms in order and superposes', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    const r = await api.atoms({
      targets: [
        { target: CA(1, 20), structure: ref },
        { target: CA(1, 20), structure: mobile },
      ],
    });
    expectOk(r);
    expect(r.value!.results[0].paired).toBe(20);
    expect(r.value!.results[0].rmsd).toBeLessThan(0.05);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
    expect(
      await api.atoms({
        targets: [
          { target: CA(1, 20), structure: ref },
          { target: CA(1, 20), structure: mobile },
        ],
      }),
    ).toMatchObject({
      ok: true,
      unchanged: true,
    });
  });

  it('fails invalid for different counts and for too few atoms', async () => {
    const { ref, mobile, api } = await twoCopies();
    const counts = await api.atoms({
      targets: [
        { target: CA(1, 20), structure: ref },
        { target: CA(1, 10), structure: mobile },
      ],
    });
    expect(counts).toMatchObject({ ok: false, code: 'invalid', message: expect.stringContaining('different numbers') });
    const few = await api.atoms({
      targets: [
        { target: CA(1, 2), structure: ref },
        { target: CA(1, 2), structure: mobile },
      ],
    });
    expect(few).toMatchObject({ ok: false, code: 'invalid', message: expect.stringContaining('too few') });
  });
});

describe('superpose.ligands', () => {
  it('superposes the same ligand of the copies by atom names', async () => {
    const { v, ref, mobile, api } = await twoCopies();
    const r = await api.ligands({
      targets: [
        { target: 'ligand', structure: ref },
        { target: 'ligand', structure: mobile },
      ],
    });
    expectOk(r);
    expect(r.value!.results[0].matched).toBe('atom-name');
    expect(r.value!.results[0].rmsd).toBeLessThan(0.05);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
  });

  it('fails invalid for a target that is not one ligand', async () => {
    const { ref, mobile, api } = await twoCopies();
    const r = await api.ligands({
      targets: [
        { target: 'ligand', structure: ref },
        { target: 'polymer', structure: mobile },
      ],
    });
    expect(r).toMatchObject({ ok: false, code: 'invalid', message: expect.stringContaining('one ligand') });
  });
});

describe('superpose.db', () => {
  /** 1CRN with a SIFTS UniProt mapping added to its atoms (the columns the updated PDB files carry). */
  function withSifts(name: string): File {
    const lines = fs
      .readFileSync(path.resolve(import.meta.dirname, '../../../../../../../data/examples/1crn.cif'), 'utf8')
      .split('\n');
    const out: string[] = [];
    for (const line of lines) {
      if (line.startsWith('_atom_site.pdbx_PDB_model_num')) {
        out.push(
          line,
          '_atom_site.pdbx_sifts_xref_db_name',
          '_atom_site.pdbx_sifts_xref_db_acc',
          '_atom_site.pdbx_sifts_xref_db_num',
          '_atom_site.pdbx_sifts_xref_db_res',
        );
      } else if (line.startsWith('ATOM ') || line.startsWith('HETATM ')) {
        const seq = line.trim().split(/\s+/)[8];
        out.push(`${line.trimEnd()} UNP P01542 ${seq} T`);
      } else out.push(line);
    }
    return new File([out.join('\n')], name);
  }

  it('fails invalid, with a hint, for structures without the mapping', async () => {
    const { api } = await twoCopies(crambin);
    const r = await api.db({});
    expect(r).toMatchObject({
      ok: false,
      code: 'invalid',
      message: expect.stringContaining('no SIFTS'),
      hint: expect.any(String),
    });
  });

  it('superposes structures by the residues they share in UniProt', async () => {
    const { v, ref, mobile, api } = await twoCopies(crambin, { files: [withSifts('a.cif'), withSifts('b.cif')] });
    const r = await api.db({ structures: [ref, mobile] });
    expectOk(r);
    expect(r.value!.results[0].rmsd).toBeLessThan(0.05);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
    expect(await api.db({})).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.history.undo());
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Move, 1e-6)).toBe(true);
  });

  it('names a structure that is not there', async () => {
    const { api } = await twoCopies(crambin, { files: [withSifts('a.cif'), withSifts('b.cif')] });
    expect(await api.db({ structures: ['nope', 'a'] })).toMatchObject({ ok: false, code: 'not-found' });
  });
});

describe('the feature', () => {
  const context = { empty: false };

  it('brings its commands', async () => {
    const v = await headlessViewer({ extensions: [SuperpositionFeature] });
    for (const id of [
      'superpose.open',
      'superpose.chains',
      'superpose.tm-align',
      'superpose.ligands',
      'superpose.atoms',
      'superpose.db',
    ]) {
      expect(v.managers.commands.get(id), id).toBeDefined();
    }
  });

  it('runs a command with its targets as text, one per structure', async () => {
    const { v, mobile } = await twoCopies();
    const command = v.managers.commands.get('superpose.chains')!;
    expect(Object.keys(command.params!(context, {}))).toContain('targets');
    const r = await command.run(
      { targets: 'polymer; polymer', alignSequences: true, traceOnly: true },
      context,
      new AbortController().signal,
    );
    expect(r.undoable).toBe(true);
    expect(r.message).toMatch(/Superposed/);
    expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
  });

  it('opens the Selection tab at the section through view.open', async () => {
    const v = await headlessViewer({ extensions: [SuperpositionFeature] });
    const runs = recordRuns(v);
    // Headless has no Selection tab to open, so the command reports that; it asked for the right one.
    await expect(
      v.managers.commands.get('superpose.open')!.run({}, context, new AbortController().signal),
    ).rejects.toThrow();
    const open = runs.find((r) => r.operation === 'view.open');
    expect(open?.values).toMatchObject({ tab: 'core.tab.selection', section: 'superposition.section.superpose' });
  });
});

describe('the Superpose view model', () => {
  it('follows the structures and the selection and superposes by the form', async () => {
    const { v, ref, mobile } = await twoCopies();
    const messages: string[] = [];
    const model = createSuperpositionModel({
      api: v.api,
      sources: v.sources,
      notify: (text, kind) => messages.push(`${kind}: ${text}`),
    });
    try {
      await until(() => model.available.value && model.hint.value !== '');
      expect(model.ready.value).toBe(false);
      expect(model.hint.value).toMatch(/Select one chain/);
      expectOk(await v.api.select.set({ target: 'polymer', structure: ref }));
      expectOk(await v.api.select.set({ target: 'polymer', structure: mobile, mode: 'add' }));
      await until(() => model.ready.value);
      expect(model.entries.value.map((e) => e.id)).toEqual([ref, mobile]);
      model.superpose();
      await until(() => messages.length > 0);
      expect(messages[0]).toMatch(/^success: Superposed/);
      expect(Mat4.areEqual(matrixOf(mobileOf(v, mobile)), Undone, 0.02)).toBe(true);
      // UniProt uses every structure, selected or not.
      model.form.load({ method: 'db', alignSequences: true, traceOnly: true });
      await until(() => model.entries.value.length === 2 && /all chains/.test(model.entries.value[0].label));
      model.superpose();
      await until(() => messages.length > 1);
      expect(messages[1]).toMatch(/^error: .*no SIFTS/);
    } finally {
      model.dispose();
    }
  });
});
