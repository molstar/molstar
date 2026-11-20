/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure language's examples (locations-and-queries.md, The language) on PDB entries read whole by the mmCIF
// runtime: 4HHB's irons and its α chains' proximal His, 1FIN's glycine loop and the properties at each level, 3SGJ's
// glycans and the Asn297 they hang on. Medium lane: each file is a few hundred kilobytes parsed; the quick
// representatives are the 3SGJ glycan fixture and 1EHZ (`elements.test.ts`).
// @module-tag medium

import { parseAddress, printAddress } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { openEntry } from '../../runtimes/_test/mmcif-entries.js';
import { codes, count, ex, keysOf, loc, prop, sourceRun } from './structures.js';

describe('the structure language on PDB entries', { tags: ['medium'] }, () => {
  it('4HHB: four irons, a heme each; α’s proximal His by author chain and number in both copies (EC-64, EC-65)', async () => {
    const t = sourceRun('4hhb', (await openEntry('4HHB')).coverage.nodes ?? []);
    const c = '/data/4hhb/(node :key complex)';
    const fe = t.run('/data/4hhb/(elements :component HEM :atom FE)');
    expect(fe.result?.type === 'location' && fe.result.location.parts.map((p) => p.system)).toEqual(
      ['E', 'G', 'H', 'J'].map((u) => `${c}/${u}/(atoms)`),
    );
    expect(count(fe)).toBe(4);
    // the chains' residues are no HEM: their atoms, not held, are never asked for (no warning)
    expect(codes(fe)).toEqual([]);
    const his = t.run('/data/4hhb/(elements :auth-chain [A C] :auth-seq 87)');
    expect([loc(his, `${c}/A`), loc(his, `${c}/C`)]).toEqual([
      [86, 87],
      [86, 87],
    ]);
    expect(count(t.run('/data/4hhb/(elements :component HIS :auth-chain [A C] :auth-seq 87)'))).toBe(2);
    // the same by the long names, and by the entity
    expect(count(t.run('/data/4hhb/(elements :auth_asym_id [A C] :auth_seq_id 87 :label_entity_id 1)'))).toBe(2);
  });

  it('1FIN: the glycine loop by label numbers on chain A; a unit’s properties; groups by unit', async () => {
    const t = sourceRun('1fin', (await openEntry('1FIN')).coverage.nodes ?? []);
    const a = '/data/1fin/(node :key complex)/A';
    expect(loc(t.run('/data/1fin/(elements :chain A :beg 10 :end 20)'), a)).toEqual([9, 20]);
    // a polymer's atoms read by its residues' standard components: the eleven residues with a CA (no warning)
    const ca = t.run('/data/1fin/(elements :chain A :beg 10 :end 20 :atom CA)');
    expect([loc(ca, a), codes(ca)]).toEqual([[9, 20], []]);
    const units = t.form('/data/1fin', {
      kind: 'groups',
      params: { 'unit-test': ex('=', prop('auth-chain'), 'A'), by: 'unit' },
    });
    expect(keysOf(units)).toEqual(['/data/1fin/(node :key complex)#A', '/data/1fin/(node :key complex)#E']);
    // the ATP's atoms by element, through the molecule's class
    const p = t.form('/data/1fin/(elements :molecule-class ligand :element P)', {
      kind: 'filter',
      params: { test: ex('=', prop('component'), 'ATP') },
    });
    expect(count(p)).toBe(6);
  });

  it('1FIN: b-factor read where its channel lies, a residue’s from the mmCIF runtime: residues filtered as text', async () => {
    const source = await openEntry('1FIN');
    const unit = '/data/1fin/(node :key complex)/A';
    const read = await source.read(
      { system: 'A', start: 0, end: 298, level: 1, shapes: ['values'], field: 'b_factor' },
      new AbortController().signal,
    );
    if (!read.ok || read.tables.values?.form !== 'positions') throw new Error('no values');
    const values = read.tables.values;
    const channel = { kind: 'values', id: 'bf', params: { target: unit, source: '/data/1fin', field: 'b_factor' } };
    const data = {
      coverage: () => source.coverage,
      keysOf: () => ['A'],
      get: () => ({ status: 'loaded', tables: { values } }),
    };
    const t = sourceRun('1fin', [...(source.coverage.nodes ?? []), { key: 'bf', node: channel }], undefined, {
      data: data as never,
    });
    const text = '/data/1fin/(elements :chain A)/(filter :test (> (prop b-factor) 60))';
    expect(printAddress(parseAddress(text).pattern!.segments)).toBe(text);
    const r = t.run(text);
    expect(codes(r)).toEqual([]);
    const expected: number[] = [];
    for (let i = 0; i < values.value.length; i++) if (values.value[i] > 60) expected.push(values.start + i);
    expect(expected.length).toBeGreaterThan(0);
    const got = loc(r, unit) as number[];
    const positions: number[] = [];
    for (let k = 0; k < got.length; k += 2) for (let p = got[k]; p < got[k + 1]; p++) positions.push(p);
    expect(positions).toEqual(expected);
  });

  it('3SGJ: a C1 a glycan residue; the Asn297s and what their links join (EC-62)', async () => {
    const t = sourceRun('3sgj', (await openEntry('3SGJ')).coverage.nodes ?? []);
    const c = '/data/3sgj/(node :key complex)';
    expect(count(t.run(`${c}/D/(elements :atom C1)`))).toBe(8);
    expect(count(t.run('/data/3sgj/(elements :molecule-class branched :atom C1)'))).toBe(8 + 8 + 5 + 2);
    const asn = '/data/3sgj/(elements :auth-chain A :auth-seq 297)';
    expect(count(t.run(asn))).toBe(1);
    const glycan = t.run(`${asn}/(linked :via ${c}/links :type glycosidic)/(expand :to unit)`);
    expect(loc(glycan, `${c}/D`)).toEqual([0, 8]);
    const modified = t.form('/data/3sgj', {
      kind: 'groups',
      params: { 'molecule-test': ex('=', prop('molecule-class'), 'branched'), by: 'molecule' },
    });
    expect(keysOf(modified)).toEqual([`${c}#3`, `${c}#4`, `${c}#5`]);
  });
});
