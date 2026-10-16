/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { ModelExportExtension } from '../prism-viewer.js';

// Model export in the Prism viewer on a headless viewer (parity.md §6): 1CBS (a protein and its retinoic acid, REA) written
// as mmCIF, BinaryCIF and PDB, each loaded back. No network.

const cbs = () =>
  new File([readFileSync(new URL('../../../../data/examples/1cbs_updated.cif', import.meta.url))], '1cbs.cif');

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(files: File[], withExtension = true) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ...(withExtension ? [ModelExportExtension] : [])] },
  });
  viewers.push(v);
  for (const f of files) {
    const r = await v.api.load.files({ files: [f] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

const atoms = (v: PrismViewer) => v.sources.structures.value.map((s) => s.atoms);

async function exported(v: PrismViewer, values: Record<string, unknown>) {
  const r = await v.api['model-export']!.export(values);
  if (!r.ok) throw new Error(r.message);
  return r;
}

describe('the model-export extension in the Prism viewer', () => {
  it('exports 1CBS as mmCIF, which loads back with every atom, and is no history step', async () => {
    const v = await viewer([cbs()]);
    const steps = v.managers.history.steps.value.length;
    const r = await exported(v, { format: 'cif' });
    expect(r.value!.filename).toBe('1CBS.cif');
    expect(r.did).toMatch(/^Exported 1CBS as 1CBS\.cif \(\d+ KB\)\.$/);
    const text = await r.value!.blob.text();
    expect(text.startsWith('data_1CBS')).toBe(true);
    expect(text).toContain('_atom_site.Cartn_x');
    // Every category of the file is kept for a whole structure.
    expect(text).toContain('_struct_conf.');
    expect(v.managers.history.steps.value.length).toBe(steps);

    const back = await viewer([new File([text], 'again.cif')]);
    expect(atoms(back)).toEqual(atoms(v));
  });

  it('exports PDB records in their columns, ligands after a TER', async () => {
    const v = await viewer([cbs()]);
    const r = await exported(v, { format: 'pdb', name: 'my model.pdb' });
    expect(r.value!.filename).toBe('my_model.pdb');
    const lines = (await r.value!.blob.text()).split('\n');
    const first = lines[0];
    expect(first.startsWith('ATOM      1  N   PRO A   1')).toBe(true);
    expect(first.length).toBe(80);
    expect(first.slice(76, 78)).toBe(' N');
    const ter = lines.findIndex((l) => l.startsWith('TER'));
    expect(lines[ter - 1].startsWith('ATOM')).toBe(true);
    expect(lines.slice(ter + 1).some((l) => l.startsWith('HETATM') && l.slice(17, 20) === 'REA')).toBe(true);
    expect(lines.at(-2)).toBe('END');
    expect(r.warnings ?? []).toEqual([]);

    const back = await viewer([new File([lines.join('\n')], '1cbs.pdb')]);
    expect(atoms(back)).toEqual(atoms(v));
  });

  it('keeps a modified residue (a HETATM) inside its chain, and says coordinates the columns cannot hold', async () => {
    const atom = (rec: string, i: number, name: string, res: string, seq: number, x: number, el: string) =>
      `${rec.padEnd(6)}${String(i).padStart(5)}  ${name.padEnd(3)} ${res} A${String(seq).padStart(4)}    ${x.toFixed(3).padStart(8)}${'0.000'.padStart(8)}${'0.000'.padStart(8)}  1.00 10.00           ${el}`;
    const residues = [
      ['ATOM', 'ALA'],
      ['ATOM', 'ALA'],
      ['HETATM', 'MSE'],
      ['ATOM', 'ALA'],
      ['ATOM', 'ALA'],
    ];
    const lines: string[] = [];
    residues.forEach(([rec, res], r) => {
      for (const [k, [name, el]] of [
        ['N', 'N'],
        ['CA', 'C'],
        ['C', 'C'],
      ].entries()) {
        lines.push(atom(rec, r * 3 + k + 1, name, res, r + 1, r * 3.8 + k * 1.3, el));
      }
    });
    const v = await viewer([new File([[...lines, 'END', ''].join('\n')], 'mse.pdb')]);
    const r = await exported(v, { format: 'pdb' });
    const out = (await r.value!.blob.text()).split('\n');
    // One TER, after the last residue: not before the MSE.
    expect(out.filter((l) => l.startsWith('TER'))).toHaveLength(1);
    expect(out.findIndex((l) => l.startsWith('TER'))).toBe(out.findIndex((l) => l === 'END') - 1);
    expect(out.some((l) => l.startsWith('HETATM') && l.slice(17, 20) === 'MSE')).toBe(true);
    expect(r.warnings ?? []).toEqual([]);

    // A coordinate beyond the column is cut, which is said.
    const far = lines.map((l, i) => (i === 0 ? atom('ATOM', 1, 'N', 'ALA', 1, -1234.567, 'N') : l));
    const wide = await viewer([new File([[...far, 'END', ''].join('\n')], 'far.pdb')]);
    const w = await exported(wide, { format: 'pdb' });
    expect((w.warnings ?? []).join(' ')).toMatch(/too wide.*mmCIF/);
  });

  it('exports BinaryCIF that loads back', async () => {
    const v = await viewer([cbs()]);
    const r = await exported(v, { format: 'bcif' });
    expect(r.value!.filename).toBe('1CBS.bcif');
    const bytes = new Uint8Array(await r.value!.blob.arrayBuffer());
    const back = await viewer([new File([bytes], '1cbs.bcif')]);
    expect(atoms(back)).toEqual(atoms(v));
  });

  it('exports only the selection, and fails when nothing is selected', async () => {
    const v = await viewer([cbs()]);
    const none = await v.api['model-export']!.export({ format: 'pdb', selection: true });
    expect(none.ok).toBe(false);
    if (!none.ok) {
      expect(none.code).toBe('not-found');
      expect(none.message).toBe('Nothing is selected in 1CBS.');
    }
    expect((await v.api.select.set({ target: 'resn REA' })).ok).toBe(true);
    const r = await exported(v, { format: 'pdb', selection: true });
    expect(r.did).toMatch(/^Exported the selection of 1CBS/);
    const records = (await r.value!.blob.text()).split('\n').filter((l) => /^(ATOM|HETATM)/.test(l));
    expect(records.length).toBe(22);
    expect(records.every((l) => l.startsWith('HETATM') && l.slice(17, 20) === 'REA')).toBe(true);
    const cif = await (await exported(v, { format: 'cif', selection: true })).value!.blob.text();
    expect(cif).toContain('REA');
    expect(cif).not.toMatch(/\bPRO\b/);
  });

  it('puts several structures into a ZIP, or one by name; says what is wrong', async () => {
    const v = await viewer([cbs(), cbs()]);
    const r = await exported(v, { format: 'cif', name: 'both' });
    expect(r.value!.filename).toBe('both.zip');
    expect(r.value!.files).toHaveLength(2);
    expect(new Set(r.value!.files).size).toBe(2);
    expect(r.value!.blob.type).toBe('application/zip');
    const ref = v.plugin.managers.structure.hierarchy.current.structures[1].cell.transform.ref;
    expect((await exported(v, { structure: ref })).value!.files).toHaveLength(1);
    const missing = await v.api['model-export']!.export({ structure: '9xyz' });
    expect(!missing.ok && missing.code).toBe('not-found');
    const empty = await viewer([]);
    const blocked = await empty.api['model-export']!.export({});
    expect(!blocked.ok && blocked.code).toBe('blocked');
  });

  it('brings its commands and its skill, and asks the agent every time (no undo restores a file)', async () => {
    const v = await viewer([]);
    expect(v.managers.commands.get('model-export.export')).toBeDefined();
    expect(v.managers.commands.get('model-export.open')).toBeDefined();
    const def = v.api.describe('model-export.export');
    expect(def).toBeDefined();
    const op = v.api['model-export']!.export.definition;
    expect(op.undoable).toBe(false);
    expect(op.history).toBeUndefined();
    expect(ModelExportExtension.skills?.[0].tools.map((t) => t.id)).toEqual(['export_models']);
    const bare = await viewer([], false);
    expect(bare.api['model-export']).toBeUndefined();
  });
});
