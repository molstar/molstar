/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, StructureElement, StructureProperties as P } from '@molstar/model/model/structure';
import { beforeAll, describe, expect, it } from 'vitest';
import {
  type SelectorT,
  checkMolQL,
  compactResidues,
  compileQuery,
  formatMolQL,
  joinSelectors,
  selectLoci,
} from '../index.js';
import { tinyStructure } from './fixture.js';
import { registerQueryLanguages } from '../languages.js';
import { AllQueryLanguages } from '../languages/all.js';

// Every query language, as the default spec's viewer has them (targets/languages.ts).
registerQueryLanguages(AllQueryLanguages);

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
/** The residues a selection holds, in order: their addresses and names. */
function residues(loci: StructureElement.Loci): { address: string; comp: string }[] {
  const out = new Map<string, string>();
  StructureElement.Loci.forEachLocation(loci, (l) => {
    const address = `${P.chain.auth_asym_id(l)}:${P.residue.auth_seq_id(l)}${P.residue.pdbx_PDB_ins_code(l) || ''}`;
    if (!out.has(address)) out.set(address, P.atom.label_comp_id(l));
  });
  return [...out].map(([address, comp]) => ({ address, comp }));
}
const residuesOf = (s: SelectorT) => residues(selectLoci(structure, s)).map((r) => r.address);

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
    // A group holding a group and an operator reads (PyMOL's grammar: the operator is tried first in a group).
    expect(compileQuery('byres ((chain A and polymer) within 4 of (chain A and resn HEM))').ok).toBe(true);
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
    const cys = residues(selectLoci(structure, pymol('resn CYS and chain A')));
    expect(cys.every((r) => r.comp === 'CYS')).toBe(true);
    expect(cys).toHaveLength(3);
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
