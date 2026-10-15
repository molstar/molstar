/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import { beforeAll, describe, expect, it } from 'vitest';
import {
  type SelectorT,
  checkMolQL,
  compactResidues,
  compileQuery,
  formatMolQL,
  joinSelectors,
  project,
  selectLoci,
  selectionItems,
} from '../index.js';
import { tinyStructure } from './fixture.js';

// Stories Untold's `any-molql` tests, on a structure built in `fixture.ts` (the heme A:500 beside CYS 104).

let structure: Structure;
beforeAll(async () => {
  structure = await tinyStructure();
});
const pymol = (text: string): SelectorT => {
  const r = compileQuery(text);
  if (!r.ok) throw new Error(r.message);
  return { molql: r.molql };
};
const residuesOf = (s: SelectorT) =>
  project(selectLoci(structure, s), 'residues', 1000).rows.map((r) => (r.kind === 'residue' ? r.address : ''));

describe('compiling selection text to MolQL', () => {
  it('reads PyMOL, VMD, Jmol and MolScript', () => {
    for (const [text, language] of [
      ['byres polymer within 5 of resn HEM', 'pymol'],
      ['protein and within 5 of resname HEM', 'vmd'],
      ['[HEM] or resno=104', 'jmol'],
      [
        '(structure-query.generator.atom-groups :residue-test (core.rel.eq (structure-query.atom-property.macromolecular.label_comp_id) HEM))',
        'molscript',
      ],
      ['(sel.atom.atom-groups :residue-test (= atom.label_comp_id HEM))', 'molscript'],
    ] as const) {
      const r = compileQuery(text, language);
      expect(r.ok, `${language}: ${r.ok ? '' : r.message}`).toBe(true);
    }
  });

  it('says why it cannot read a query, with a hint', () => {
    const bad = compileQuery('byres polymer within five of resn HEM');
    expect(bad.ok).toBe(false);
    if (!bad.ok) {
      expect(bad.message).toMatch(/^PyMOL could not read it/);
      expect(bad.hint).toMatch(/PyMOL selection syntax/);
    }
    expect(compileQuery('   ')).toMatchObject({ ok: false, message: 'The query is empty.' });
    const lopsided = compileQuery('byres (polymer and within 5 of (chain A and resi 104))');
    expect(!lopsided.ok && lopsided.hint).toMatch(/give "within" and "around" a left side/i);
    const doubled = compileQuery('byres ((chain A and polymer) within 4 of (chain A and resn HEM))');
    expect(!doubled.ok && doubled.hint).toMatch(/open no group with two parentheses/);
    expect(checkMolQL({ head: { name: 'not.a.symbol' } })).toMatch(/Invalid MolQL|Mol\* cannot run it/);
  });

  it('formats MolQL as MolScript text, which reads back to the same selection', () => {
    const heme = pymol('resn HEM');
    const text = formatMolQL((heme as { molql: unknown }).molql);
    expect(text).toMatch(/^\(structure-query\./);
    const again = compileQuery(text, 'molscript');
    expect(again.ok && residuesOf({ molql: again.molql })).toEqual(['A:500']);
  });
});

describe('evaluating selections (a tiny structure with a heme beside CYS 104)', () => {
  it('evaluates every selector form alike', () => {
    expect(residuesOf({ auth_asym_id: 'A', auth_seq_id: 104 })).toEqual(['A:104']);
    expect(residuesOf([{ auth_asym_id: 'A', beg_auth_seq_id: 102, end_auth_seq_id: 105 }])).toEqual([
      'A:102',
      'A:103',
      'A:104',
      'A:105',
    ]);
    expect(residuesOf('ligand')).toEqual(['A:500']);
    expect(residuesOf(pymol('resn HEM'))).toEqual(['A:500']);
  });

  it("selects what the old forms cannot say: the heme's surroundings, the cysteines", () => {
    const around = residuesOf(pymol('byres polymer within 5 of resn HEM'));
    expect(around).toContain('A:104');
    expect(around.length).toBeGreaterThan(2);
    expect(around).not.toContain('A:500');
    expect(around).not.toContain('A:119');
    const cys = project(selectLoci(structure, pymol('resn CYS and chain A')), 'residues');
    expect(cys.rows.every((r) => r.kind === 'residue' && r.comp === 'CYS')).toBe(true);
    expect(cys.residues).toBe(3);
  });

  it('joins selections: items stay a list, anything else becomes one MolQL', () => {
    const items = joinSelectors([{ auth_asym_id: 'A', auth_seq_id: 104 }, [{ auth_asym_id: 'A', auth_seq_id: 105 }]]);
    expect(Array.isArray(items)).toBe(true);
    const mixed = joinSelectors([{ auth_asym_id: 'A', auth_seq_id: 104 }, pymol('resn HEM')]);
    expect(mixed && typeof mixed === 'object' && 'molql' in mixed).toBe(true);
    expect(residuesOf(mixed as SelectorT)).toEqual(['A:104', 'A:500']);
    expect(joinSelectors(['ligand'])).toBe('ligand');
    expect(joinSelectors([])).toBeUndefined();
  });

  it('projects a selection to atoms, residues, chains and molecules', () => {
    const loci = selectLoci(structure, pymol('resn HEM or (chain A and resi 104)'));
    const atoms = project(loci, 'atoms', 5);
    expect(atoms.rows[0]).toMatchObject({ kind: 'atom', element: expect.any(String) });
    expect(atoms.more).toBe(atoms.atoms - 5);
    const iron = project(selectLoci(structure, pymol('resn HEM and name FE')), 'atoms');
    expect(iron.rows).toEqual([
      { kind: 'atom', address: 'A:500@FE', element: 'FE', comp: 'HEM', xyz: expect.any(Array) },
    ]);
    const chains = project(loci, 'chains');
    // The heme sits on author chain A but in its own label chain.
    expect(chains.rows.map((r) => (r.kind === 'chain' ? `${r.auth}/${r.label}` : ''))).toEqual(['A/A', 'A/D']);
    const molecules = project(loci, 'molecules');
    expect(molecules.rows.map((r) => (r.kind === 'molecule' ? r.type : ''))).toEqual(['polymer', 'non-polymer']);
    expect(project(loci).line).toBe(`2 residues (${atoms.atoms} atoms) in chain A`);
  });

  it('says when it selects nothing', () => {
    expect(project(selectLoci(structure, pymol('chain Z'))).line).toBe('selects nothing');
  });
});

describe('residue addresses', () => {
  it('fold consecutive residues of a chain into ranges; insertion codes stay single', () => {
    expect(compactResidues(['A:1', 'A:2', 'A:3', 'A:5', 'B:6', 'B:7', 'A:100A', 'A:101'])).toEqual([
      'A:1-3',
      'A:5',
      'B:6-7',
      'A:100A',
      'A:101',
    ]);
  });
});

describe('what a selection holds, as schema items', () => {
  it('keys residues by label chain and author number, apart from the author chain they share', () => {
    const heme = selectionItems(selectLoci(structure, pymol('resn HEM')));
    expect(heme.label_comp_id).toEqual(['HEM']);
    // The heme is in author chain A with the protein, but in its own label chain.
    expect(heme.auth_asym_id).toEqual(['A']);
    expect(heme.label_asym_id?.[0]).not.toBe('A');
    expect(heme.atom_id).toBeUndefined();
    const pocket = selectionItems(selectLoci(structure, pymol('byres polymer within 4 of resn HEM')));
    expect(pocket.label_asym_id?.length).toBe(pocket.auth_seq_id?.length);
    expect(new Set(pocket.label_asym_id)).toEqual(new Set(['A']));
  });

  it('keeps insertion codes when the selection has one', () => {
    const items = selectionItems(
      selectLoci(structure, { auth_asym_id: 'A', beg_auth_seq_id: 109, end_auth_seq_id: 111 }),
    );
    expect(items.auth_seq_id).toContain(110);
    expect(items.pdbx_PDB_ins_code).toContain('A');
  });

  it('lists atoms when the selection cuts through residues', () => {
    const iron = selectionItems(selectLoci(structure, pymol('resn HEM and name FE')));
    expect(iron.atom_id?.length).toBe(1);
    expect(iron.auth_atom_id).toEqual(['FE']);
  });
});
