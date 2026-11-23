/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The mmCIF runtime (_spec/README.md "Data", mmCIF): 1EHZ read against the counts recorded
// from the file (the quick representative of the PDB entries, `mmcif-pdb.test.ts`), and a small hand-written file
// for what no entry has: a polymer without canonical letters, an altloc, a hydrogen, an assembly of two operators,
// waters, pLDDT, the source's answers and its refusals.

import { gzipSync } from 'node:zlib';
import { PCG } from '@molstar/core/data/util/hash-functions';
import type { CategoryColumn, Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { DefaultOpenLimits, type OpenContext, type OpenedSource } from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';
import { assemblyOperators, operatorGroups } from '../mmcif-description.js';
import { type MmcifOptions, MmcifRuntime, mmcifRuntime, openMmcif } from '../mmcif.js';
import {
  countedFromFile,
  expectLinksAsCounted,
  expectValid,
  frameOf,
  openEntry,
  Recorded,
  summary,
} from './mmcif-entries.js';

const signal = () => new AbortController().signal;
const ctx = { signal: signal(), limits: DefaultOpenLimits };
const Address = '/data/t';
const values = (c: unknown) => {
  const col = c as CategoryColumn;
  return Array.from({ length: col.length }, (_, i) => col.get(i));
};

/**
 * A small entry: a peptide `MET MSE XYZ GLY` (no canonical letters; XYZ unknown; GLY 4A; XYZ unmodeled), a ligand
 * `LIG` (C1 at two alternate locations, a hydrogen), two waters, a covalent link from GLY's C to the ligand, the
 * peptide's own backbone bond at MSE, a metal on a water, assembly 1 of operators 1 and 2, and pLDDT.
 */
const Small = `data_TEST
_entry.id TEST
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'A peptide'
2 non-polymer 'A ligand'
3 water water
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 MET n
1 2 MSE n
1 3 XYZ n
1 4 GLY n
loop_
_pdbx_entity_nonpoly.entity_id
_pdbx_entity_nonpoly.comp_id
2 LIG
3 HOH
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
C 3
loop_
_chem_comp.id
_chem_comp.type
_chem_comp.name
GLY 'peptide linking' GLYCINE
HOH non-polymer WATER
LIG non-polymer 'A LIGAND'
MET 'L-peptide linking' METHIONINE
MSE 'L-peptide linking' SELENOMETHIONINE
XYZ 'L-peptide linking' 'AN UNKNOWN ONE'
loop_
_pdbx_struct_mod_residue.label_comp_id
_pdbx_struct_mod_residue.parent_comp_id
MSE MET
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM N N . MET A 1 1 ? 1.0 10.0 1 P 1
ATOM C CA . MET A 1 1 ? 1.0 20.0 1 P 1
ATOM H H . MET A 1 1 ? 1.0 90.0 1 P 1
ATOM N N . MSE A 1 2 ? 1.0 30.0 2 P 1
ATOM N N . GLY A 1 4 A 0.5 40.0 3 P 1
ATOM C C . GLY A 1 4 A 1.0 50.0 3 P 1
HETATM C C1 A LIG B 2 . ? 0.6 5.0 101 P 1
HETATM C C1 B LIG B 2 . ? 0.4 7.0 101 P 1
HETATM O O1 . LIG B 2 . ? 1.0 5.0 101 P 1
HETATM H H1 . LIG B 2 . ? 1.0 5.0 101 P 1
HETATM O O . HOH C 3 . ? 1.0 30.0 201 P 1
HETATM O O . HOH C 3 . ? 1.0 31.0 202 P 1
ATOM N N . MET A 1 1 ? 1.0 99.0 1 P 2
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_symmetry
_struct_conn.pdbx_role
covale A 1 C 1 1_555 A 2 N 2 1_555 ?
covale A 4 C 3 1_555 B . C1 101 1_555 ?
metalc B . O1 101 1_555 C . O 201 1_555 ?
disulf A 1 SG 1 1_555 A 4 SG 3 2_655 ?
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1-2)' A,B
loop_
_ma_qa_metric.id
_ma_qa_metric.name
_ma_qa_metric.type
1 pLDDT pLDDT
loop_
_ma_qa_metric_local.label_asym_id
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.metric_id
_ma_qa_metric_local.metric_value
A 1 1 91.5
A 2 1 72.0
`;

const bytes = (text = Small) => new TextEncoder().encode(text);
const open = (options: MmcifOptions = {}, text = Small) =>
  openMmcif(bytes(text), ctx, { address: Address, format: 'mmcif' }, options);
const node = (s: OpenedSource, key: string) => s.coverage.nodes!.find((n) => n.key === key)?.node as Node;
const units = (s: OpenedSource) => node(s, 'complex').children!.filter((c) => c.kind === 'unit');
const linkRows = (s: OpenedSource) => {
  const rows = node(s, 'complex').children!.find((c) => c.id === 'links')?.params?.rows as Record<string, any>;
  if (!rows) return [];
  return Array.from(
    { length: rows.a_start.length },
    (_, i) => `${rows.a_unit.get(i)}:${rows.a_start[i]} ${rows.type.get(i)} ${rows.b_unit.get(i)}:${rows.b_start[i]}`,
  );
};
const opened = async (p: Promise<OpenedSource>) => {
  try {
    await p;
  } catch (e) {
    if (e instanceof SourceOpenError) return e.issues.map((i) => i.code);
    throw e;
  }
  throw new Error('opened');
};

describe('the mmCIF runtime on 1EHZ (the quick representative)', () => {
  it('every entity, chain, metal link, modified residue and author number as recorded, counted again from the file', async () => {
    const source = await openEntry('1EHZ');
    expect(summary(source)).toEqual(Recorded['1EHZ']);
    const counted = countedFromFile(await frameOf('1EHZ'));
    const s = Recorded['1EHZ'];
    expect(Object.keys(s.sequences).length + Object.keys(s.molecules).length).toBe(counted.entities);
    expect(Object.keys(s.units).length).toBe(counted.chains);
    expect(counted.lengths).toEqual({ '1': 76 });
    // 47 metalc rows, 36 on a water; 25 covale rows, each within the chain (its backbone at a modified residue)
    expect([counted.conn.metalc, counted.conn.covale]).toEqual([
      { rows: 47, onWater: 36, oneChain: 0 },
      { rows: 25, onWater: 0, oneChain: 25 },
    ]);
    expectLinksAsCounted(s, counted);
    expectValid(source.coverage.nodes!);
  });

  it('the modified nucleotides: each its component, its parent from pdbx_struct_mod_residue', async () => {
    const source = await openEntry('1EHZ');
    const modified = node(source, '1').params!.modified as { position: Int32Array; component: CategoryColumn };
    expect(Array.from(modified.position).slice(0, 4)).toEqual([9, 15, 16, 25]);
    expect(values(modified.component).slice(0, 4)).toEqual(['2MG', 'H2U', 'H2U', 'M2G']);
    const compounds = node(source, 'compounds').children!;
    expect(compounds.find((c) => c.id === 'YYG')?.params?.parents).toEqual(['G']);
    expect(compounds.find((c) => c.id === 'MG')?.params?.parents).toBeUndefined();
  });
});

describe('the mmCIF runtime', () => {
  it('reads a polymer without canonical letters from its components: a parent’s letter, X for none, said (an entity without canonical letters)', async () => {
    const s = await open();
    expect(node(s, '1').params).toMatchObject({ alphabet: 'protein', letters: 'MMXG' });
    const modified = node(s, '1').params!.modified as { position: Int32Array; component: CategoryColumn };
    expect(Array.from(modified.position)).toEqual([1, 2]);
    expect(values(modified.component)).toEqual(['MSE', 'XYZ']);
    expect(s.issues?.map((i) => i.code)).toContain('seq.io.structure.letters');
  });

  it('a unit’s observed runs and auth numbering from atom_site, an insertion code, the first model only (EC-118, EC-2)', async () => {
    const [a, b] = units(await open());
    expect(a.params).toMatchObject({ chain: 'A', auth_chain: 'P' });
    const observed = a.params!.observed as { start: Int32Array; end: Int32Array };
    expect([Array.from(observed.start), Array.from(observed.end)]).toEqual([
      [0, 3],
      [2, 4],
    ]);
    expect(values(a.children![0].params!.labels)).toEqual(['1', '2', null, '3A']);
    expect(b.id).toBe('B');
    expect(b.params?.observed).toBeUndefined();
  });

  it('a molecule’s atoms: the first label_alt_id of each name, the other dropped and said; a hydrogen refused (EC-64)', async () => {
    const s = await open();
    const atoms = node(s, '2').params!.atoms as { name: CategoryColumn; element: CategoryColumn };
    expect(values(atoms.name)).toEqual(['C1', 'O1']);
    expect(s.issues?.map((i) => i.code)).toEqual(
      expect.arrayContaining(['seq.io.mmcif.altloc', 'seq.io.structure.hydrogen']),
    );
    expect(s.issues?.find((i) => i.code === 'seq.io.mmcif.altloc')?.message).toMatch(/^1 atoms/);
    expect(s.issues?.every((i) => i.address === Address)).toBe(true);
    expect(node(s, '2').params).toMatchObject({ class: 'ligand', name: 'A ligand' });
  });

  it('links: a covalent one kept, a backbone bond, a symmetry mate’s and a water’s not, said once', async () => {
    const s = await open();
    expect(linkRows(s)).toEqual(['A:3 covalent B:0']);
    expect(s.issues?.find((i) => i.code === 'seq.io.mmcif.links')?.message).toMatch(/^2 struct_conn/);
  });

  it('waters only when asked: a molecule of a water a position, its unit, the metal link on it', async () => {
    const dry = await open();
    expect(node(dry, '3')).toBeUndefined();
    const wet = await open({ water: true });
    expect(node(wet, '3').params).toMatchObject({ class: 'water' });
    expect(values(node(wet, '3').params!.residues)).toEqual(['HOH', 'HOH']);
    expect(units(wet).map((u) => u.id)).toEqual(['A', 'B', 'C']);
    expect(linkRows(wet)).toEqual(['A:3 covalent B:0', 'B:0 metal C:0']);
    expectValid(wet.coverage.nodes!);
  });

  it('an assembly: its chains under each operator, a unit each keyed by chain and operator, links within each copy (EC-65)', async () => {
    const s = await open({ assembly: '1' });
    expect(units(s).map((u) => u.id)).toEqual(['A-ASM-1', 'A-ASM-2', 'B-ASM-1', 'B-ASM-2']);
    expect(units(s)[1].params).toMatchObject({ chain: 'A', operator: 'ASM-2' });
    expect(linkRows(s)).toEqual(['A-ASM-1:3 covalent B-ASM-1:0', 'A-ASM-2:3 covalent B-ASM-2:0']);
    expectValid(s.coverage.nodes!);
    // a source's own `assembly` wins over the runtime's
    const own = await openMmcif(bytes(), ctx, { address: Address, format: 'mmcif', assembly: '1' } as never);
    expect(own.coverage.systems.length).toBe(4);
  });

  it('refuses an assembly the file lacks, with its own code', async () => {
    expect(await opened(open({ assembly: '9' }))).toEqual(['seq.io.mmcif.assembly']);
  });

  it('operator expressions as Mol*’s instance ids: lists, ranges, products, bounded', () => {
    expect(operatorGroups('1,2', 10)).toEqual({ groups: [['1', '2']], bounded: false });
    expect(operatorGroups('(1-3)', 10)).toEqual({ groups: [['1', '2', '3']], bounded: false });
    expect(assemblyOperators('(X0)(1-2)', 10)).toEqual({ operators: ['ASM-X0-1', 'ASM-X0-2'], bounded: false });
    expect(assemblyOperators('(1-60)(61-88)', 5)).toMatchObject({ bounded: true });
    expect(assemblyOperators('(1-60)(61-88)', 5).operators.length).toBe(5);
  });

  const crafted = Array.from({ length: 20 }, () => '1-999999').join(',');

  it('a crafted operator expression is bounded as it is expanded, never made whole', () => {
    expect(crafted.length).toBe(179);
    const { groups, bounded } = operatorGroups(crafted, 10_001);
    expect([groups.flat().length, bounded]).toEqual([10_001, true]);
  });

  // the open makes its 10,000 units (a quarter of a second): medium
  it('a crafted operator expression is bounded as the file opens, and said', { tags: ['medium'] }, async () => {
    const s = await open({ assembly: '1' }, Small.replace("1 '(1-2)' A,B", `1 '${crafted}' A,B`));
    expect(s.issues?.map((i) => i.code)).toContain('seq.io.mmcif.operators');
    expect(s.issues?.map((i) => i.code)).toContain('seq.io.structure.bounds');
    expect(s.coverage.systems.length).toBe(10_000);
  });

  it('serves values by field over a unit’s positions: pLDDT the default, B-factor the heavy atoms’ mean, NaN where unmodeled', async () => {
    const s = await open();
    expect(s.coverage.levels).toEqual([1]);
    expect(s.coverage.shapes[0].fields?.map((f) => f.name)).toEqual(['plddt', 'b_factor', 'occupancy']);
    const read = async (field?: string, start = 0, end = 4) => {
      const a = await s.read(
        { system: 'A', start, end, level: 1, shapes: ['values'], ...(field ? { field } : {}) },
        signal(),
      );
      if (!a.ok || a.tables.values?.form !== 'positions') throw new Error(JSON.stringify(a));
      return Array.from(a.tables.values.value);
    };
    expect(await read()).toEqual([91.5, 72, Number.NaN, Number.NaN]);
    expect(await read('b_factor')).toEqual([15, 30, Number.NaN, 45]);
    expect(await read('occupancy', 3, 9)).toEqual([0.75]);
  });

  it('answers not-served for a unit, a level or a field it lacks; empty tables for shapes it does not serve', async () => {
    const s = await open();
    const ask = (system: string, level = 1, field?: string) =>
      s.read({ system, start: 0, end: 2, level, shapes: ['values'], ...(field ? { field } : {}) }, signal());
    for (const a of [await ask('Z'), await ask('A', 2), await ask('A', 1, 'charge')])
      expect(a).toMatchObject({ ok: false, kind: 'not-served' });
    expect(await s.read({ system: 'A', start: 0, end: 2, level: 1, shapes: ['letters'] }, signal())).toEqual({
      ok: true,
      tables: {},
    });
    const aborted = new AbortController();
    aborted.abort();
    await expect(
      s.read({ system: 'A', start: 0, end: 2, level: 1, shapes: ['values'] }, aborted.signal),
    ).rejects.toThrow();
  });

  it('a branched molecule’s name: its description, or its composition when asked', async () => {
    const text = Small.replace("2 non-polymer 'A ligand'", "2 branched 'A ligand'").replace(
      'loop_\n_pdbx_entity_nonpoly.entity_id\n_pdbx_entity_nonpoly.comp_id\n2 LIG\n3 HOH',
      'loop_\n_pdbx_entity_nonpoly.entity_id\n_pdbx_entity_nonpoly.comp_id\n3 HOH\nloop_\n_pdbx_entity_branch_list.entity_id\n_pdbx_entity_branch_list.num\n_pdbx_entity_branch_list.comp_id\n2 1 NAG\n2 2 NAG\n2 3 BMA',
    );
    expect(node(await open({}, text), '2').params).toMatchObject({ class: 'branched', name: 'A ligand' });
    expect(node(await open({ name: 'composition' }, text), '2').params?.name).toBe('NAG2 BMA');
  });

  it('reads gzipped bytes as the text; refuses bytes past the bound, text that does not parse, a source without a uri', async () => {
    const plain = await open();
    const zipped = await openMmcif(gzipSync(bytes()), ctx, { address: Address, format: 'mmcif' });
    expect(zipped.coverage).toEqual(plain.coverage);
    expect(
      await opened(openMmcif(bytes(), { signal: signal(), limits: { ...DefaultOpenLimits, maxBytes: 100 } })),
    ).toEqual(['seq.source.unreadable']);
    expect(
      await opened(
        openMmcif(gzipSync(bytes()), { signal: signal(), limits: { ...DefaultOpenLimits, maxBytes: 1000 } }),
      ),
    ).toEqual(['seq.source.unreadable']);
    for (const broken of ['data_x\nloop_\n_a.b\n_a.c\n1 2 3\n', 'not a CIF file'])
      expect(await opened(open({}, broken))).toEqual(['seq.source.malformed']);
    const none: OpenContext = {
      assets: { readText: async () => '', readBytes: async () => new Uint8Array() },
      signal: signal(),
      limits: DefaultOpenLimits,
    };
    expect(await opened(MmcifRuntime.open({ address: Address, format: 'mmcif' }, none))).toEqual([
      'seq.source.no_data',
    ]);
  });

  it('opens a source by uri through the asset port, with the runtime’s own params', async () => {
    const asked: string[] = [];
    const assets: OpenContext = {
      assets: {
        readText: async () => '',
        readBytes: async (uri) => {
          asked.push(uri);
          if (uri !== 'small.cif') throw new Error('no such file');
          return bytes();
        },
      },
      signal: signal(),
      limits: DefaultOpenLimits,
    };
    const runtime = mmcifRuntime({ water: true });
    expect([runtime.name, runtime.formats]).toEqual(['mmcif', ['mmcif']]);
    const s = await runtime.open({ address: Address, format: 'mmcif', uri: 'small.cif' }, assets);
    expect(units(s).map((u) => u.id)).toEqual(['A', 'B', 'C']);
    expect(await opened(runtime.open({ address: Address, format: 'mmcif', uri: 'gone.cif' }, assets))).toEqual([
      'seq.source.unreadable',
    ]);
    expect(asked).toEqual(['small.cif', 'gone.cif']);
  });

  // a short fuzz quick, the full one medium
  for (const [runs, tags] of [
    [12, []],
    [60, ['medium']],
  ] as [number, string[]][])
    it(`a fuzz of ${runs}: seeded malformed files (categories dropped, columns short, rows repeated, huge numbers) never throw past a coded refusal, and give valid nodes`, {
      tags,
    }, async () => {
      const rng = new PCG(5);
      const pick = (n: number) => rng.int() % n;
      const lines = Small.split('\n');
      for (let run = 0; run < runs; run++) {
        const out = [...lines];
        for (let k = 0, n = 1 + pick(6); k < n; k++) {
          const i = 1 + pick(out.length - 1);
          const what = pick(5);
          if (what === 0)
            out.splice(i, 1); // a line gone: a column short, a category's header lost
          else if (what === 1)
            out.splice(i, 0, out[i]); // a row or a field name repeated
          else if (what === 2)
            out[i] = out[i].replace(/\b\d+\b/g, '2000000000'); // huge numbers
          else if (what === 3)
            out[i] = out[i].replace(/\S+$/, '?'); // a value unknown
          else
            out.splice(
              i,
              0,
              'loop_',
              '_pdbx_struct_assembly_gen.assembly_id',
              '_pdbx_struct_assembly_gen.oper_expression',
              '_pdbx_struct_assembly_gen.asym_id_list',
              "1 '(1-2000000000)(X)' A,B,Q",
            );
        }
        const text = out.join('\n');
        for (const options of [{}, { assembly: '1', water: true }] as MmcifOptions[]) {
          try {
            const s = await open(options, text);
            expectValid(s.coverage.nodes!);
            for (const i of s.issues ?? []) expect(i.code).toMatch(/^seq\.(io|source)\./);
          } catch (e) {
            if (!(e instanceof SourceOpenError)) throw new Error(`run ${run}: ${String(e)}`);
            for (const i of e.issues) expect(i.code).toMatch(/^seq\.(io|source)\./);
            // a refusal is the reading's (the parser's, the bound's), never a crash of the converter caught
            expect(e.issues[0].message, `run ${run}`).not.toMatch(/Cannot read|is not a function|is not iterable/);
          }
        }
      }
    });

  it('an aborted open rejects with an AbortError, never an issue', async () => {
    const c = new AbortController();
    c.abort();
    await expect(openMmcif(bytes(), { signal: c.signal, limits: DefaultOpenLimits })).rejects.toMatchObject({
      name: 'AbortError',
    });
  });
});
