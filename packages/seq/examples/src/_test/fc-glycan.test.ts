/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { decodeSeq } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import {
  Asn297,
  Disulfides,
  FcChains,
  FcFile,
  FcLength,
  FcSystem,
  fcAuth,
  fcPosition,
  GlycanLinks,
  GlycanResidues,
  GlycanWurcs,
  ZincHistidines,
} from '../antibody/fc-glycan.js';
import { PackageFolder } from '../node-assets.js';
import { structureSeq } from '../p53/structure-record.js';

// The Antibody's Fc glycan, PDB 3AVE: what the view holds of the entry against the Fc's
// committed `seq` file (quick), and against the whole mmCIF (medium): the glycan's residues, links and WURCS, the units'
// observed runs and author numbers, the complex's links, the file written again.

const text = (f: string) => readFileSync(`${PackageFolder}${f}`, 'utf8');
const fc = decodeSeq(text(FcFile)).data!;
const letters = (() => {
  const a = fc.pairs[0].answer;
  return a.ok ? (a.tables.letters as { letters: string }).letters : '';
})();

describe('the Fc glycan (3AVE)', () => {
  it('the Fc: 223 residues; N297 in its sequon N-S-T, C at each disulfide’s ends, H at the zinc’s', () => {
    expect(fc.coverage.systems).toEqual([{ name: FcSystem, length: FcLength, alphabet: 'protein' }]);
    const at = (auth: number) => letters[fcPosition(auth)];
    expect(letters.slice(fcPosition(Asn297), fcPosition(Asn297) + 3)).toBe('NST');
    expect(Disulfides.flat().map(at).join('')).toBe('CCCC');
    expect(ZincHistidines.map(at).join('')).toBe('HH');
    expect(fcAuth(0)).toBe(225);
  });

  it('sources.json: the file’s bytes and checksum, from the recorded mmCIF', () => {
    const row = (
      JSON.parse(text('fixtures/antibody/sources.json')) as {
        file: string;
        bytesFixture: number;
        sha256Fixture: string;
        sha256Download: string;
      }[]
    ).find((r) => `fixtures/${r.file}` === FcFile)!;
    const entry = (
      JSON.parse(text('fixtures/structures/sources.json')) as { file: string; sha256Inflated: string }[]
    ).find((e) => e.file === 'structures/3AVE.cif.gz')!;
    const b = readFileSync(`${PackageFolder}${FcFile}`);
    expect([row.bytesFixture, row.sha256Fixture, row.sha256Download]).toEqual([
      b.byteLength,
      createHash('sha256').update(b).digest('hex'),
      entry.sha256Inflated,
    ]);
  });

  // the whole entry parsed: slow by nature, the medium lane (the committed file above is the quick representative)
  it('the entry’s glycan, links, observed runs and author numbers are the view’s; the file written again', {
    tags: ['medium'],
  }, async () => {
    // the parser imported here: only this medium test needs it (the quick lane does not load @molstar/io)
    const { readCif } = await import('../../scripts/cif.js');
    const cif = await readCif(gunzipSync(readFileSync(`${PackageFolder}fixtures/structures/3AVE.cif.gz`)));
    expect(structureSeq(cif, FcSystem)).toBe(text(FcFile));
    // an IgG1 Fc (IGHG1_HUMAN) with the glycan on Asn297
    expect(cif.rows('struct_ref')[0].db_code).toBe('IGHG1_HUMAN');
    const list = cif.rows('pdbx_entity_branch_list').filter((r) => r.entity_id === '2');
    expect(list.map((r) => r.comp_id)).toEqual([...GlycanResidues]);
    const links = cif.rows('pdbx_entity_branch_link').map((r) => ({
      a: Number(r.entity_branch_list_num_1) - 1,
      b: Number(r.entity_branch_list_num_2) - 1,
      atom: r.atom_id_2,
    }));
    expect(links).toEqual(GlycanLinks.map((l) => ({ ...l })));
    expect(cif.rows('pdbx_entity_branch_descriptor').find((r) => r.type === 'WURCS')!.descriptor).toBe(GlycanWurcs);
    const conn = (type: string) =>
      cif
        .rows('struct_conn')
        .filter((r) => r.conn_type_id === type)
        .map((r) => [
          r.ptnr1_label_asym_id,
          Number(r.ptnr1_auth_seq_id),
          r.ptnr2_label_asym_id,
          Number(r.ptnr2_auth_seq_id),
        ]);
    expect(conn('covale').filter((r) => r[0] !== r[2])).toEqual(FcChains.map((c) => [c.chain, Asn297, c.glycan, 1]));
    expect(conn('disulf')).toEqual(FcChains.flatMap((c) => Disulfides.map(([x, y]) => [c.chain, x, c.chain, y])));
    expect(conn('metalc').filter((r) => r[0] === 'A')).toEqual(ZincHistidines.map((h) => ['A', h, 'E', 448]));
    // each chain's modeled residues and author numbers (atom_site, model 1)
    for (const ch of FcChains) {
      const seen = new Map<number, string>();
      for (const r of cif.rows('atom_site'))
        if (r.label_asym_id === ch.chain && r.pdbx_PDB_model_num === '1')
          seen.set(Number(r.label_seq_id) - 1, r.auth_seq_id);
      const at = [...seen.keys()].sort((x, y) => x - y);
      expect([at[0], at[at.length - 1] + 1, at.length]).toEqual([
        ch.observed[0],
        ch.observed[1],
        ch.observed[1] - ch.observed[0],
      ]);
      for (const [p, auth] of seen) expect(Number(auth)).toBe(fcAuth(p));
    }
  });
});
