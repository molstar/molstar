/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure language on a small fixture recorded from 3SGJ (the Fc's Asn297 and its N-glycan, `structures.ts`)
// and on 1EHZ read by the mmCIF runtime (the quick representative; 4HHB, 1FIN and 3SGJ whole are the medium lane's,
// `elements-pdb.test.ts`): `elements` at each depth, over a selection as its prefix, `groups` by level, `filter`,
// `expand` and `linked`; a lacking property a warning and nothing, never a failure (monadic).

import { parseAddress, printAddress } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { openEntry } from '../../runtimes/_test/mmcif-entries.js';
import { codes, count, ex, glycanSpec, keysOf, loc, prop, sourceRun, specRun } from './structures.js';

const g = specRun(glycanSpec());
const D = '/data/c/D';
const DA = `${D}/(atoms)`;

describe('elements on the 3SGJ glycan fixture', () => {
  it('an atom field: positions on the atoms system, one C1 a residue, a heavy atom each (EC-64)', () => {
    expect(count(g.run('/data/glycan/(elements :atom C1)'))).toBe(8);
    expect(loc(g.run(`${D}/(elements :atom C1)`), DA)).toEqual([
      0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19, 21, 22,
    ]);
    // the long name reads as the short one
    expect(g.run(`${D}/(elements :label_atom_id C1)`).result).toEqual(g.run(`${D}/(elements :atom C1)`).result);
  });

  it('over a selection: the C1s of the residues it holds (monadic over a Location prefix)', () => {
    expect(loc(g.run(`${D}/(range :from 2 :to 4)/(elements :atom C1)`), DA)).toEqual([3, 4, 6, 7, 9, 10]);
    expect(loc(g.run(`${D}/(elements :component MAN)/(elements :atom C1)`), DA)).toEqual([9, 10, 15, 16]);
  });

  it('no residue or atom field: the matching units’ whole residues systems (EC-35); residue fields: those positions; author ends in the numbering’s order (EC-2, EC-6)', () => {
    expect(loc(g.run('/data/c/(elements :chain D)'), D)).toEqual([0, 8]);
    expect(loc(g.run('/data/c/(elements :chain A :auth-seq 297)'), '/data/c/A')).toEqual([4, 5]);
    expect(loc(g.run('/data/c/(elements :chain A :beg 3 :end 5)'), '/data/c/A')).toEqual([2, 5]);
    expect(loc(g.run('/data/c/(elements :auth-chain A :auth-beg 295 :auth-end 298)'), '/data/c/A')).toEqual([2, 6]);
  });

  it('columns and records: one table read alike (SoA ≡ AoS)', () => {
    const cols = g.run('/data/c/(elements :chain [A D] :seq [5 1])');
    const recs = g.run('/data/c/(elements :items [(fields :chain A :seq 5) (fields :chain D :seq 1)])');
    expect(cols.result).toEqual(recs.result);
    expect([loc(cols, '/data/c/A'), loc(cols, D)]).toEqual([
      [4, 5],
      [0, 1],
    ]);
  });

  it('linked: the Asn and what the links join to it, found from either end (EC-62), then expanded to the glycan’s unit', () => {
    const asn = '/data/c/(elements :chain A :auth-seq 297)';
    const l = g.run(`${asn}/(linked :via /data/c/links)`);
    expect([loc(l, '/data/c/A'), loc(l, D)]).toEqual([
      [4, 5],
      [0, 1],
    ]);
    const whole = g.run(`${asn}/(linked :via /data/c/links :type glycosidic)/(expand :to unit)`);
    expect([loc(whole, '/data/c/A'), loc(whole, D)]).toEqual([
      [0, 9],
      [0, 8],
    ]);
    expect(loc(g.run(`${asn}/(linked :via /data/c/links :type disulfide)`), D)).toEqual([]);
  });

  it('expand: atoms to their residues', () => {
    expect(loc(g.run(`${D}/(elements :atom C6)/(expand :to residue)`), D)).toEqual([7, 8]);
  });

  it('a lacking property: a warning naming it and the count, the rest kept', () => {
    const r = g.run('/data/c/(elements :auth-seq 297)');
    expect(loc(r, '/data/c/A')).toEqual([4, 5]);
    expect(codes(r)).toEqual(['seq.query.property.absent']);
    expect(r.issues[0].message).toMatch(/1 element lack auth-seq/);
  });

  it('filter: per element, a predicate with positional operands (built in code until the grammar reads them)', () => {
    // a positional operand on a form that is no expression
    expect(codes(g.form(`${D}/(elements)`, ex('filter', ex('=', prop('element'), 'C'))))).toEqual([
      'seq.query.form.param',
    ]);
    const c6 = g.form(`${D}/(elements :atom [C1 C6])`, {
      kind: 'filter',
      params: { test: ex('=', prop('atom'), 'C6') },
    });
    expect(loc(c6, DA)).toEqual([23, 24]);
    // over nodes: a unit passes by its unit-level test
    const units = g.form('/(merge :of [/data/c/A /data/c/D])', {
      kind: 'filter',
      params: { test: ex('=', prop('chain'), 'D') },
    });
    expect(units.result?.type === 'nodes' && units.result.nodes.map((n) => n.address.join('/'))).toEqual(['data/c/D']);
    const over = g.form(`${D}/(elements)`, {
      kind: 'filter',
      params: { test: ex('in', prop('component'), ['MAN', 'FUC']) },
    });
    expect(loc(over, D)).toEqual([3, 4, 5, 6, 7, 8]);
  });

  it('as text: filter and groups with predicates parse, resolve and print back the same', () => {
    const c6 = `${D}/(elements :atom [C1 C6])/(filter :test (= (prop atom) C6))`;
    expect(loc(g.run(c6), DA)).toEqual([23, 24]);
    expect(printAddress(parseAddress(c6).pattern!.segments)).toBe(c6);
    const man = '/data/c/(groups :by unit :residue-test (in (prop component) [MAN FUC]))';
    expect(keysOf(g.run(man))).toEqual(['/data/c#D']);
    // a positional operand on a form that is no expression
    expect(codes(g.run(`${D}/(elements C1)`))).toEqual(['seq.query.form.param']);
  });

  it('groups: a key a residue, an atom, a unit; a lower level’s property in a test refused', () => {
    const man = { kind: 'groups', params: { 'residue-test': ex('=', prop('component'), 'MAN') } };
    expect(keysOf(g.form('/data/c', man))).toEqual([`${D}#3`, `${D}#5`]);
    const units = { kind: 'groups', params: { 'residue-test': ex('=', prop('component'), 'MAN'), by: 'unit' } };
    expect(keysOf(g.form('/data/c', units))).toEqual(['/data/c#D']);
    const atoms = { kind: 'groups', params: { 'atom-test': ex('=', prop('atom'), 'C6'), by: 'atom' } };
    expect(keysOf(g.form('/data/c', atoms))).toEqual([`${DA}#23`]);
    const wrong = { kind: 'groups', params: { 'unit-test': ex('=', prop('component'), 'MAN') } };
    expect(codes(g.form('/data/c', wrong))).toEqual(['seq.query.property.level']);
  });
});

describe('elements on 1EHZ (read by the mmCIF runtime)', async () => {
  const source = await openEntry('1EHZ');
  const t = sourceRun('1ehz', source.coverage.nodes ?? []);
  const unit = '/data/1ehz/(node :key complex)/A';

  it('modified residues: a key each, the 14 rows of `pdbx_struct_mod_residue`', () => {
    const r = t.run('/data/1ehz/(groups :residue-test (= (prop modified) true))');
    expect(keysOf(r)).toEqual([9, 15, 16, 25, 31, 33, 36, 38, 39, 45, 48, 53, 54, 57].map((p) => `${unit}#${p}`));
    expect(count(t.run('/data/1ehz/(elements :modified true)'))).toBe(14);
    // Y37 is YYG: by its component
    expect(loc(t.run('/data/1ehz/(elements :component YYG)'), unit)).toEqual([36, 37]);
  });

  it('a nucleotide by its author number on the unit', () => {
    expect(loc(t.run('/data/1ehz/(elements :auth-chain A :auth-seq 37)'), unit)).toEqual([36, 37]);
  });
});
