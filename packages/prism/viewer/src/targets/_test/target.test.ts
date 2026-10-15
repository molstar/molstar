/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { beforeAll, describe, expect, it } from 'vitest';
import {
  type ComponentExpression,
  KEYWORDS,
  KEYWORD_EXPRESSIONS,
  compileTarget,
  formatMolQL,
  formatTarget,
  fromMolqlJson,
  parseAddressList,
  parseTarget,
  toMolqlJson,
} from '../index.js';
import { tinyStructure } from './fixture.js';
import { evaluateTarget } from '../index.js';
import { registerQueryLanguages } from '../languages.js';
import { AllQueryLanguages } from '../languages/all.js';

// Every query language, as the default spec's viewer has them (targets/languages.ts).
registerQueryLanguages(AllQueryLanguages);

function items(text: string): ComponentExpression[] {
  const r = parseAddressList(text);
  if (!r.ok) throw new Error(`${r.message} ${r.hint}`);
  return r.items;
}

describe('keywords', () => {
  it('are read in any case, with spaces and a plural', () => {
    for (const keyword of KEYWORDS) {
      expect(parseTarget(keyword)).toEqual({ ok: true, parsed: { kind: 'keyword', keyword } });
      expect(parseTarget(`  ${keyword.toUpperCase()} `)).toEqual({ ok: true, parsed: { kind: 'keyword', keyword } });
    }
    expect(parseTarget('Ligands')).toEqual({ ok: true, parsed: { kind: 'keyword', keyword: 'ligand' } });
    expect(parseTarget('waters')).toEqual({ ok: true, parsed: { kind: 'keyword', keyword: 'water' } });
  });

  it('compile to their fixed expressions', () => {
    for (const keyword of KEYWORDS) {
      const c = compileTarget(keyword);
      expect(c.ok && c.molql).toBe(KEYWORD_EXPRESSIONS[keyword]);
    }
  });

  it('selection and focus are the person: parsed, resolved by the viewer', () => {
    expect(parseTarget('selection')).toEqual({ ok: true, parsed: { kind: 'person', which: 'selection' } });
    expect(parseTarget('Focus')).toEqual({ ok: true, parsed: { kind: 'person', which: 'focus' } });
    const without = compileTarget('selection');
    expect(without).toMatchObject({ ok: false, message: expect.stringMatching(/only the viewer knows/) });
    expect(!without.ok && without.hint).toMatch(/Name what you mean/);

    const asked: string[] = [];
    const ligand = KEYWORD_EXPRESSIONS.ligand;
    const resolved = compileTarget('selection', {
      resolvePerson: (which) => {
        asked.push(which);
        return which === 'selection' ? ligand : undefined;
      },
    });
    expect(resolved.ok && resolved.molql).toBe(ligand);
    expect(asked).toEqual(['selection']);
    const nothing = compileTarget('focus', { resolvePerson: () => undefined });
    expect(nothing).toMatchObject({ ok: false, message: 'Nothing is focused.' });
  });
});

describe('address lists', () => {
  it('read the permissive forms', () => {
    expect(items('A')).toEqual([{ auth_asym_id: 'A' }]);
    expect(items('B:120')).toEqual([{ auth_asym_id: 'B', auth_seq_id: 120 }]);
    expect(items('B:200-210')).toEqual([{ auth_asym_id: 'B', beg_auth_seq_id: 200, end_auth_seq_id: 210 }]);
    expect(items('B:120A')).toEqual([{ auth_asym_id: 'B', auth_seq_id: 120, pdbx_PDB_ins_code: 'A' }]);
    expect(items('C:122@OG1')).toEqual([{ auth_asym_id: 'C', auth_seq_id: 122, auth_atom_id: 'OG1' }]);
    expect(items('A/120')).toEqual([{ auth_asym_id: 'A', auth_seq_id: 120 }]);
    expect(items('A 120-130')).toEqual([{ auth_asym_id: 'A', beg_auth_seq_id: 120, end_auth_seq_id: 130 }]);
    expect(items('A120')).toEqual([{ auth_asym_id: 'A', auth_seq_id: 120 }]);
    expect(items('A:120..130')).toEqual([{ auth_asym_id: 'A', beg_auth_seq_id: 120, end_auth_seq_id: 130 }]);
    expect(items('A:120 to 130')).toEqual([{ auth_asym_id: 'A', beg_auth_seq_id: 120, end_auth_seq_id: 130 }]);
    expect(items('A:-5')).toEqual([{ auth_asym_id: 'A', auth_seq_id: -5 }]);
    expect(items('A:-5-3')).toEqual([{ auth_asym_id: 'A', beg_auth_seq_id: -5, end_auth_seq_id: 3 }]);
    expect(items('120-130')).toEqual([{ beg_auth_seq_id: 120, end_auth_seq_id: 130 }]);
  });

  it('read lists, trailing commas, case and spaces', () => {
    expect(items('A, B:120,')).toEqual([{ auth_asym_id: 'A' }, { auth_asym_id: 'B', auth_seq_id: 120 }]);
    expect(items('  A : 120 - 130 ;  b:7a ')).toEqual([
      { auth_asym_id: 'A', beg_auth_seq_id: 120, end_auth_seq_id: 130 },
      { auth_asym_id: 'b', auth_seq_id: 7, pdbx_PDB_ins_code: 'A' },
    ]);
    // A range written backwards is the same range.
    expect(items('A:130-120')).toEqual([{ auth_asym_id: 'A', beg_auth_seq_id: 120, end_auth_seq_id: 130 }]);
    // Residues after an item with a chain belong to that chain.
    expect(items('A:120, 130-140, B, 7')).toEqual([
      { auth_asym_id: 'A', auth_seq_id: 120 },
      { auth_asym_id: 'A', beg_auth_seq_id: 130, end_auth_seq_id: 140 },
      { auth_asym_id: 'B' },
      { auth_asym_id: 'B', auth_seq_id: 7 },
    ]);
  });

  it('read ranges with insertion codes as ranges in the order of the file', () => {
    expect(items('A:100A-102')).toEqual([
      { auth_asym_id: 'A', beg_auth_seq_id: 100, beg_auth_ins_code: 'A', end_auth_seq_id: 102 },
    ]);
    expect(items('A:100-102b')).toEqual([
      { auth_asym_id: 'A', beg_auth_seq_id: 100, end_auth_seq_id: 102, end_auth_ins_code: 'B' },
    ]);
  });

  it('read the atom forms', () => {
    expect(items('id:1234')).toEqual([{ atom_id: 1234 }]);
    expect(items('@#1234')).toEqual([{ atom_id: 1234 }]);
    expect(items('index:1233')).toEqual([{ atom_index: 1233 }]);
    expect(items('element:fe')).toEqual([{ type_symbol: 'FE' }]);
    expect(items('id: 12')).toEqual([{ atom_id: 12 }]);
    expect(items('C:122 element:S')).toEqual([{ auth_asym_id: 'C', auth_seq_id: 122, type_symbol: 'S' }]);
    expect(items('A element:C')).toEqual([{ auth_asym_id: 'A', type_symbol: 'C' }]);
    expect(items('A@CA')).toEqual([{ auth_asym_id: 'A', auth_atom_id: 'CA' }]);
    expect(items("A:5@o5'")).toEqual([{ auth_asym_id: 'A', auth_seq_id: 5, auth_atom_id: "O5'" }]);
    expect(items('A[ASM-X0-1]:120')).toEqual([{ auth_asym_id: 'A', instance_id: 'ASM-X0-1', auth_seq_id: 120 }]);
    expect(items('A[1_555]')).toEqual([{ auth_asym_id: 'A', instance_id: '1_555' }]);
    expect(items('id:1, id:2')).toEqual([{ atom_id: 1 }, { atom_id: 2 }]);
  });

  it('are not mistaken for query text, nor query text for addresses', () => {
    for (const text of ['chain A', 'resi 100-120', 'resn HIS', 'name CA', 'chain A and resi 100']) {
      const r = parseTarget(text);
      expect(r.ok && r.parsed.kind, text).toBe('query');
    }
    for (const text of ['A', 'B:120', 'A 120-130', 'C:122@OG1']) {
      const r = parseTarget(text);
      expect(r.ok && r.parsed.kind, text).toBe('addresses');
    }
  });

  it('say which part did not read, with a hint', () => {
    for (const [text, part] of [
      ['B:12x3', '12x3'],
      ['A:100-', '100-'],
      ['id:abc', 'abc'],
      ['index:1.5', '1.5'],
      ['element:F3e9', 'F3e9'],
      ['A:10, B:x', 'x'],
      ['C:122@', '@'],
    ] as const) {
      const r = parseTarget(text);
      expect(r.ok, text).toBe(false);
      if (!r.ok) {
        expect(r.message, text).toContain(part);
        expect(r.hint, text).toMatch(/\S/);
      }
    }
    const twice = parseAddressList('id:1 id:2');
    expect(twice.ok).toBe(false);
    expect(!twice.ok && twice.message).toMatch(/atom_id twice/);
    expect(parseAddressList(' , ')).toMatchObject({ ok: false, message: 'The address list is empty.' });
  });
});

describe('query text', () => {
  it('is read in PyMOL by default, or the language asked for', () => {
    const pymol = compileTarget('chain A and resn HIS');
    expect(pymol.ok && pymol.parsed).toEqual({ kind: 'query', language: 'pymol', text: 'chain A and resn HIS' });
    const vmd = compileTarget('protein and within 5 of resname HEM', { language: 'vmd' });
    expect(vmd.ok && vmd.parsed.kind).toBe('query');
    const jmol = compileTarget('[HEM] or resno=104', { language: 'jmol' });
    expect(jmol.ok).toBe(true);
    const molscript = compileTarget('(sel.atom.atom-groups :residue-test (= atom.label_comp_id HEM))', {
      language: 'molscript',
    });
    expect(molscript.ok).toBe(true);
    // PyMOL text is not Jmol.
    expect(compileTarget('byres polymer within 5 of resn HEM', { language: 'jmol' }).ok).toBe(false);
  });

  it('fails with the language, the reason and a hint', () => {
    const r = compileTarget('byres polymer within five of resn HEM');
    expect(r.ok).toBe(false);
    if (!r.ok) {
      expect(r.message).toMatch(/^PyMOL could not read it/);
      expect(r.hint).toMatch(/PyMOL selection syntax/);
    }
    const empty = parseTarget('  ');
    expect(empty).toMatchObject({ ok: false, message: 'The target is empty.' });
  });

  it('fails an address-like text with the address message and the query one', () => {
    const r = compileTarget('B:12x3');
    expect(r.ok).toBe(false);
    if (!r.ok) {
      expect(r.message).toContain('12x3');
      expect(r.hint).toMatch(/Not read as PyMOL either/);
    }
  });
});

describe('selectors', () => {
  it('compile as MVS does: schema items, lists, { molql }', () => {
    expect(compileTarget({ auth_asym_id: 'A', auth_seq_id: 104 }).ok).toBe(true);
    expect(compileTarget([{ auth_asym_id: 'A' }, { auth_asym_id: 'B' }]).ok).toBe(true);
    expect(compileTarget({ molql: KEYWORD_EXPRESSIONS.water }).ok).toBe(true);
  });

  it('fail with the problem and the fields', () => {
    const unknown = compileTarget({ chain: 'A' } as never);
    expect(unknown).toMatchObject({ ok: false, message: '"chain" is not a field of a selector.' });
    expect(!unknown.ok && unknown.hint).toMatch(/auth_asym_id/);
    expect(compileTarget({ auth_seq_id: '104' } as never)).toMatchObject({
      ok: false,
      message: expect.stringMatching(/wrong type/),
    });
    expect(compileTarget([])).toMatchObject({ ok: false, message: 'The list of selectors is empty.' });
    expect(compileTarget({} as never).ok).toBe(false);
    expect(compileTarget({ molql: { head: { name: 'not.a.symbol' } } }).ok).toBe(false);
    expect(compileTarget({ molql: 'x' }).ok).toBe(false);
  });
});

describe('MolQL as it is kept', () => {
  const stored = (target: string) => {
    const c = compileTarget(target);
    if (!c.ok) throw new Error(c.message);
    return toMolqlJson(c.molql);
  };

  it('is plain JSON that survives a round trip and selects the same', async () => {
    const structure = await tinyStructure();
    for (const target of ['ligand', 'A:100-104, B:200', 'chain A and resn HIS', 'C:1@P', 'element:ZN']) {
      const json = stored(target);
      const text = JSON.stringify(json);
      const back = fromMolqlJson(text);
      expect(back.ok, target).toBe(true);
      expect(back.ok && JSON.parse(JSON.stringify(back.molql)), target).toEqual(json.molql);
      const a = evaluateTarget(structure, target);
      const b = evaluateTarget(structure, json);
      expect(a.ok && b.ok && b.selected, target).toEqual(a.ok && a.selected);
    }
  });

  it('copies the expression and takes { molql }, the expression, or either as text', () => {
    const expression = MS.struct.generator.all();
    const json = toMolqlJson(expression);
    expect(json.molql).toEqual(expression);
    expect(json.molql).not.toBe(expression);
    for (const form of [json, expression, JSON.stringify(json), JSON.stringify(expression)]) {
      expect(fromMolqlJson(form)).toEqual({ ok: true, molql: expression });
    }
  });

  it('fails with a hint on what is not MolQL', () => {
    expect(fromMolqlJson('{not json')).toMatchObject({ ok: false, message: expect.stringMatching(/not JSON/) });
    expect(fromMolqlJson({ molql: [1, [2]] })).toMatchObject({
      ok: false,
      message: 'The stored MolQL is not an expression.',
    });
    const bad = fromMolqlJson({ molql: { head: { name: 'not.a.symbol' } } });
    expect(bad.ok).toBe(false);
    expect(!bad.ok && bad.hint).toMatch(/molql/);
  });

  it('is shown as MolScript, of a target or of what was stored', () => {
    expect(formatTarget('ligand')).toMatch(/^\(structure-query\.modifier\.union/);
    const json = stored('A:104');
    expect(formatTarget(json)).toBe(formatMolQL(json.molql));
    expect(formatTarget('A:104')).toBe(formatMolQL(json.molql));
    // What does not compile is shown as it was written.
    expect(formatTarget('B:12x3')).toBe('B:12x3');
    expect(formatTarget('selection')).toBe('selection');
  });
});

describe('what it selected', () => {
  let structure: Structure;
  beforeAll(async () => {
    structure = await tinyStructure();
  });
  const select = (target: Parameters<typeof evaluateTarget>[1], options?: Parameters<typeof evaluateTarget>[2]) => {
    const r = evaluateTarget(structure, target, options);
    if (!r.ok) throw new Error(`${r.message} ${r.hint}`);
    return r.selected;
  };

  it('counts atoms and residues, names chains, folds residues by author number into ranges', () => {
    const s = select('A:100-104, A:110A, A:112, B');
    expect(s.chains).toEqual(['A', 'B']);
    expect(s.byChain.map((c) => [c.chain, c.residues])).toEqual([
      ['A', 7],
      // chain B: its 10 residues and the zinc ion of the author chain
      ['B', 11],
    ]);
    expect(s.byChain[0]?.ranges).toEqual(['100-104', '110A', '112']);
    expect(s.byChain[1]?.ranges).toEqual(['200-209', '301']);
    expect(s.addresses).toEqual(['A:100-104', 'A:110A', 'A:112', 'B:200-209', 'B:301']);
    expect(s.residues).toBe(18);
    expect(s.line).toBe(`18 residues (${s.atoms} atoms) in chains A and B`);
    expect(s.instances).toEqual([]);
  });

  it('is the same however the target was written', () => {
    const a = select('A:102-104');
    expect(select({ auth_asym_id: 'A', beg_auth_seq_id: 102, end_auth_seq_id: 104 })).toEqual(a);
    expect(select('chain A and resi 102-104')).toEqual(a);
    expect(a.byChain).toEqual([{ chain: 'A', residues: 3, atoms: 19, ranges: ['102-104'], more: 0 }]);
  });

  it('knows the keywords of this structure', () => {
    const count = (keyword: string) => select(keyword).residues;
    expect(count('all')).toBe(39);
    expect(count('protein')).toBe(31);
    expect(count('nucleic')).toBe(4);
    expect(count('polymer')).toBe(35);
    expect(count('ligand')).toBe(1);
    expect(count('ion')).toBe(1);
    expect(count('water')).toBe(2);
    expect(count('branched')).toBe(0);
    expect(count('coarse')).toBe(0);
    expect(select('branched').line).toBe('selects nothing');
  });

  it('selects atoms by the atom forms', () => {
    expect(select('A:103@OG1')).toMatchObject({ atoms: 1, residues: 1, addresses: ['A:103'] });
    expect(select('A:103 element:O')).toMatchObject({ atoms: 2, residues: 1 });
    expect(select('A:104 element:S').atoms).toBe(1);
    expect(select('element:FE')).toMatchObject({ atoms: 1, addresses: ['A:500'] });
    expect(select('element:ZN')).toMatchObject({ atoms: 1, addresses: ['B:301'] });
    expect(select('id:1').atoms).toBe(1);
    expect(select('@#1').atoms).toBe(1);
    expect(select('id:1, id:2, id:3').atoms).toBe(3);
    expect(select('id:9999')).toMatchObject({ atoms: 0, line: 'selects nothing', chains: [], addresses: [] });
    // Atom 1 is the first (id 1, index 0) atom of the file.
    expect(select('index:0').addresses).toEqual(select('id:1').addresses);
    expect(select('index:0').addresses).toEqual(['A:100']);
    const instance = structure.units[0]?.conformation.operator.instanceId as string;
    expect(select(`A[${instance}]:100`)).toMatchObject({ residues: 1, addresses: ['A:100'] });
    expect(select('A[nowhere]:100').atoms).toBe(0);
  });

  it('limits the ranges it lists', () => {
    const s = select('A:100, A:102, A:104, A:106, A:108', { maxRangesPerChain: 2 });
    expect(s.byChain[0]).toMatchObject({ residues: 5, ranges: ['100', '102'], more: 3 });
    expect(s.addresses).toEqual(['A:100', 'A:102']);
  });

  it('comes back as a failure, not a throw, when the target does not compile', () => {
    expect(evaluateTarget(structure, 'B:12x3')).toMatchObject({ ok: false });
    expect(evaluateTarget(structure, 'selection')).toMatchObject({ ok: false });
    const person = evaluateTarget(structure, 'selection', { resolvePerson: () => KEYWORD_EXPRESSIONS.ion });
    expect(person.ok && person.selected.addresses).toEqual(['B:301']);
  });

  it('gives the loci and the target it read', () => {
    const r = evaluateTarget(structure, 'ligand');
    expect(r.ok && r.target.parsed).toEqual({ kind: 'keyword', keyword: 'ligand' });
    expect(r.ok && r.loci.elements.length).toBe(1);
  });
});
