/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `molstar:secondary_structure` over a recorded structure (2OCJ, p53's core domain): helices, strands and turns as
// ranges on the polymer units, one row per element of the assignment: the model's own, or the drawn structure's computed
// one (DSSP for a file without one) where the host gives it; no network.

import { readFileSync } from 'node:fs';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { compile, type Loaders } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { categoryTable, readStructure } from '../index.js';
import { structureOf } from '../source.js';
import { bytesOf, loaders, registries } from './fixtures.js';

/**
 * The annotation's rows by unit (or molecule, `copies: "one"`): `start-stop category`, 0-based and half-open on the
 * unit's label_seq_id - 1.
 */
async function secondary(source: object, load: Loaders, copies?: 'one') {
  const b = createSeqViewSpec();
  const m = b.data.complex({ ref: 'c', provider: 'molstar:structure', source } as never);
  const options = { of: m, ...(copies ? { copies } : {}) };
  b.data.computed({ ref: 'ss', provider: 'molstar:secondary_structure', options } as never);
  b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler();
  const c = await compile(b.getDocument(), registries(), load);
  const set = c.set;
  const a = set.annotations.find((x) => x && set.strings.get(x.key).endsWith('ss'));
  const byUnit = new Map<string, string[]>();
  for (let r = 0; a && r < a.count; r++) {
    const u = a.target.unit[r];
    const unit = set.strings.get(u >= 0 ? set.units.key[u] : (set.molecules[a.target.molecule[r]]?.key ?? -1));
    const rows = byUnit.get(unit) ?? [];
    rows.push(`${a.start[r]}-${a.stop[r]} ${set.strings.get(a.category?.[r] ?? -1)}`);
    byUnit.set(unit, rows);
  }
  return { issues: c.issues, annotation: a, byUnit };
}

const PDBE_2OCJ = { service: 'pdbe', resource: 'structure', id: '2ocj' };
const kindOf = (confType: string) => (confType.startsWith('HELX') ? 'helix' : confType === 'STRN' ? 'strand' : 'turn');

describe('molstar:secondary_structure', () => {
  it("gives one range per element of the model's assignment, as the file's struct_conf has them", async () => {
    const { issues, byUnit } = await secondary(PDBE_2OCJ, loaders());
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    // The file's own struct_conf (2OCJ lists every strand there too, beside struct_sheet_range): each row is one range
    // on its chain's unit, label_seq_id beg..end as [beg - 1, end).
    const src = await readStructure(bytesOf('2ocj.bcif'), 'bcif');
    const conf = categoryTable(src.frame as NonNullable<typeof src.frame>, 'struct_conf');
    const expected = new Map<string, string[]>();
    for (let r = 0; conf && r < conf.length; r++) {
      const unit = `c/${conf.columns.beg_label_asym_id[r]}`;
      const rows = expected.get(unit) ?? [];
      const beg = Number(conf.columns.beg_label_seq_id[r]);
      const end = Number(conf.columns.end_label_seq_id[r]);
      rows.push(`${beg - 1}-${end} ${kindOf(conf.columns.conf_type_id[r] ?? '')}`);
      expected.set(unit, rows);
    }
    expect([...expected.keys()]).toEqual(['c/A', 'c/B', 'c/C', 'c/D']);
    for (const [unit, rows] of expected) expect(byUnit.get(unit)).toEqual(rows);
  });

  it('keeps elements that touch apart: a bend beside a turn, a helix beside a turn, in each chain', async () => {
    const { byUnit } = await secondary(PDBE_2OCJ, loaders());
    const window = (unit: string, lo: number, hi: number) =>
      (byUnit.get(unit) ?? []).filter((r) => {
        const start = Number(r.split('-')[0]);
        return start >= lo && start < hi;
      });
    // Chain A, label_seq_id 71–90: BEND 71–72, TURN 73–77, BEND 78, the PPII helix 79–80, the α-helix 84–88, BEND
    // 89–90.
    expect(window('c/A', 68, 90)).toEqual([
      '70-72 turn',
      '72-77 turn',
      '77-78 turn',
      '78-80 helix',
      '83-88 helix',
      '88-90 turn',
    ]);
    // Chain B has a 3₁₀ helix 73–75 where A has a turn, and its α-helix ends at 87: a turn 88, a bend 89–90 after
    // it.
    expect(window('c/B', 68, 90)).toEqual([
      '70-72 turn',
      '72-75 helix',
      '75-77 turn',
      '77-78 turn',
      '78-80 helix',
      '83-87 helix',
      '87-88 turn',
      '88-90 turn',
    ]);
    // The C-terminal α-helix, PRO 185–GLU 194, in every chain.
    for (const unit of ['c/A', 'c/B', 'c/C', 'c/D']) expect(byUnit.get(unit)).toContain('184-194 helix');
    // Chain D: a one-residue strand (VAL 79) right before a one-residue PPII helix (VAL 80).
    expect(window('c/D', 78, 80)).toEqual(['78-79 strand', '79-80 helix']);
  });

  it('gives a chain with no assignment no rows: the zinc ions here, every chain of a model without one', async () => {
    const { byUnit } = await secondary(PDBE_2OCJ, loaders());
    for (const unit of ['c/E', 'c/F', 'c/G', 'c/H']) expect(byUnit.has(unit)).toBe(false);
    // ins1.cif has no struct_conf or struct_sheet_range: the annotation is there, empty, and nothing is reported.
    const ins1 = bytesOf('ins1.cif');
    const none = await secondary({ uri: 'https://example.org/ins1.cif', format: 'mmcif' }, {
      load: async () => ins1,
    } as Loaders);
    expect(none.issues).toEqual([]);
    expect(none.annotation?.count).toBe(0);
  });

  it("reads the 3D view's assignment where the host gives its drawn structure: DSSP where the file has none", async () => {
    // 1CRN without its helices and strands (`struct_conf`, `struct_sheet_range`), as a predicted model's file is.
    const text = readFileSync(new URL('../../../../../../../data/examples/1crn.cif', import.meta.url), 'utf8');
    const stripped = text
      .split(/^#\s*$/m)
      .filter((block) => !/^\s*(loop_\s*)?_struct_(conf|sheet_range)\./.test(block))
      .join('#\n');
    expect(stripped).not.toContain('_struct_conf.');
    const bytes = new TextEncoder().encode(stripped);
    const src = await readStructure(bytes, 'mmcif');
    const host = {
      load: async () => bytes,
      parsed: async () => ({ kind: 'object', value: src }),
    } as unknown as Loaders;
    const source = { host_file: 'structure:1crn', format: 'mmcif' };
    // No drawn structure: the file's assignment, none.
    expect((await secondary(source, host)).annotation?.count).toBe(0);
    // The host's drawn structure: Mol*'s computed assignment (DSSP here, as the cartoon draws it).
    src.drawn = await structureOf(src, 0, null);
    const { byUnit } = await secondary(source, host);
    expect([...byUnit.values()].flat().some((r) => r.endsWith('helix'))).toBe(true);
  });

  it('copies "one": one row set per molecule, on it, residue by residue what most of its chains have', async () => {
    const each = (await secondary(PDBE_2OCJ, loaders())).byUnit;
    const { issues, byUnit } = await secondary(PDBE_2OCJ, loaders(), 'one');
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    expect([...byUnit.keys()]).toEqual(['c/entity:1']);
    const rows = byUnit.get('c/entity:1') ?? [];
    const window = (lo: number, hi: number) =>
      rows.filter((r) => Number(r.split('-')[0]) >= lo && Number(r.split('-')[0]) < hi);
    // 71–90 (the chains above): B's 3₁₀ helix 73–75 is outvoted by the turns of A, C and D, a turn B starts at 76
    // runs on (one chain of four starts it); A's α-helix to 88 is the others' turn at 88; all four start a turn at 89.
    // D's strand 79 is outvoted by the PPII helix.
    expect(window(68, 90)).toEqual([
      '70-72 turn',
      '72-77 turn',
      '77-78 turn',
      '78-80 helix',
      '83-87 helix',
      '87-88 turn',
      '88-90 turn',
    ]);
    // Two of four: the first chain's (A has the turn 46, C and D have none).
    expect(each.get('c/A')).toContain('45-46 turn');
    expect(each.get('c/C')).not.toContain('45-46 turn');
    expect(rows).toContain('45-46 turn');
    // One of four: D's turn 169 is out; three of four: A's turn from 136 is the others' from 135.
    expect(rows.some((r) => r.startsWith('168-'))).toBe(false);
    expect(rows).toContain('134-136 turn');
    expect(rows).toContain('184-194 helix');
    // From 94 on, chain B's rows: it has what most chains have there.
    expect(rows.filter((r) => Number(r.split('-')[0]) >= 93)).toEqual(
      (each.get('c/B') ?? []).filter((r) => Number(r.split('-')[0]) >= 93),
    );
  });

  it("keeps a file's bends and drops DSSP's, and ignores a drawn model whose residues are not the source's", async () => {
    // 2OCJ drawn: Mol*'s computed assignment is the file's own, bends and all.
    const bytes = bytesOf('2ocj.bcif');
    const src = await readStructure(bytes, 'bcif');
    const host = {
      load: async () => bytes,
      parsed: async () => ({ kind: 'object', value: src }),
    } as unknown as Loaders;
    const source = { host_file: 'structure:2ocj', format: 'bcif' };
    const own = (await secondary(source, host)).byUnit.get('c/A');
    src.drawn = await structureOf(src, 0, null);
    expect((await secondary(source, host)).byUnit.get('c/A')).toEqual(own);
    // Another file's structure as the drawn one (another hierarchy): its assignment is not read onto these residues.
    const crn = readFileSync(new URL('../../../../../../../data/examples/1crn.cif', import.meta.url));
    const other = await readStructure(new Uint8Array(crn), 'mmcif');
    src.drawn = await structureOf(other, 0, null);
    expect((await secondary(source, host)).byUnit.get('c/A')).toEqual(own);
  });
});
