/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { exampleFile, expectOk, headlessViewer } from '../../_test/harness.js';
import { describe, expect, it } from 'vitest';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { evaluateExpression, evaluateTarget } from '../evaluate.js';

// Coarse models (integrative, PDB-IHM): spheres and gaussians, no atoms. What a target selects in them is said by their
// label chains and the sequence ranges each element covers (targets.md §2).

describe('targets on a coarse model', () => {
  it('summarizes coarse elements by label chain and the residues they cover', async () => {
    const v = await headlessViewer({ files: [exampleFile('ihm-coarse-8zzh.cif')] });
    const structure = v.plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data;

    const all = evaluateTarget(structure, 'all');
    expect(all.ok).toBe(true);
    if (!all.ok) return;
    expect(all.selected.atoms).toBe(all.selected.residues);
    expect([...all.selected.chains].sort()).toEqual(['A', 'B', 'C', 'D']);

    const coarse = evaluateTarget(structure, 'coarse');
    expect(coarse.ok && coarse.selected.atoms).toBe(all.selected.atoms);

    // Chain A by its label (coarse elements have no author chain).
    const a = evaluateExpression(
      structure,
      MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.ammp('label_asym_id'), 'A']) }),
    );
    expect(a.ok).toBe(true);
    if (!a.ok) return;
    expect(a.selected.chains).toEqual(['A']);
    expect(a.selected.byChain[0].ranges.length).toBeGreaterThan(0);

    expectOk(await v.api.select.set({ target: 'coarse' }));
  });
});
