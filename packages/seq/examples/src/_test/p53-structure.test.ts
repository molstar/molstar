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
import { PackageFolder } from '../node-assets.js';
import { P53 } from '../p53/data.js';
import {
  AlignedSegments,
  ChainTableFile,
  Chains,
  EntityLength,
  EntityStart,
  Observed,
  positionOf,
  StructureFile,
  ZincSites,
} from '../p53/structure.js';
import {
  alignedSegments,
  chainTable,
  chainTableJson,
  readChainTable,
  StructureSystem,
  structureSeq,
} from '../p53/structure-record.js';

// p53's structure, PDB 2OCJ: what the view holds of the entry (its units' observed runs and author
// numbers, the zincs and their ligands, the alignment's segments) against the two committed files the record script
// writes from the one mmCIF (the entity's `seq` and the fake structure viewer's chain table); the files written again
// from the whole mmCIF in the medium lane.

const text = (f: string) => readFileSync(`${PackageFolder}${f}`, 'utf8');
const table = readChainTable(text(ChainTableFile));
const entity = decodeSeq(text(StructureFile)).data!;
const letters = (() => {
  const a = entity.pairs[0].answer;
  return a.ok ? (a.tables.letters as { letters: string }).letters : '';
})();

describe('p53’s structure (2OCJ)', () => {
  it('the entity: 219 residues, P04637 94–312 letter for letter; its global alignment one run, the view’s segments', () => {
    expect(entity.coverage.systems).toEqual([{ name: StructureSystem, length: EntityLength, alphabet: 'protein' }]);
    expect(letters).toBe(P53.slice(EntityStart, EntityStart + EntityLength));
    expect(alignedSegments(P53, letters)).toEqual({
      from: [...AlignedSegments.from],
      to: [...AlignedSegments.to],
      length: [...AlignedSegments.length],
    });
  });

  it('the chain table: four chains modeled 96–289 (EC-118), author numbers UniProt’s (EC-3); a zinc each, 501–504', () => {
    expect(table.chains.map((c) => [c.label_asym_id, c.auth_asym_id, c.entity])).toEqual([
      ...Chains.map((c) => [c.chain, c.chain, '1']),
      ...Chains.map((c) => [c.zinc, c.chain, '2']),
    ]);
    for (const { chain, zinc, zincAuth } of Chains) {
      const c = table.chains.find((x) => x.label_asym_id === chain)!;
      expect(c.residues).toHaveLength(EntityLength);
      const observed = c.residues.flatMap((r, i) => (r.observed ? [i] : []));
      expect([observed[0], observed[observed.length - 1] + 1, observed.length]).toEqual([
        Observed[0],
        Observed[1],
        Observed[1] - Observed[0],
      ]);
      for (const r of c.residues) expect(Number(r.auth)).toBe(EntityStart + r.label!);
      const z = table.chains.find((x) => x.label_asym_id === zinc)!;
      expect(z.residues.map((r) => [r.comp, r.auth, r.observed])).toEqual([['ZN', zincAuth, true]]);
    }
  });

  it('each zinc’s ligands are C176, H179, C238 and C242, inside every chain’s modeled range', () => {
    expect(ZincSites.map((s) => letters[positionOf(s.auth)]).join('')).toBe('CHCC');
    for (const s of ZincSites) {
      const p = positionOf(s.auth);
      expect(p >= Observed[0] && p < Observed[1]).toBe(true);
    }
  });

  it('sources.json: each file’s bytes and checksum are the file’s, from the recorded mmCIF', () => {
    const rows = JSON.parse(text('fixtures/p53/sources.json')) as {
      file: string;
      bytesFixture: number;
      sha256Fixture: string;
      sha256Download: string;
    }[];
    const entry = (
      JSON.parse(text('fixtures/structures/sources.json')) as { file: string; sha256Inflated: string }[]
    ).find((e) => e.file === 'structures/2OCJ.cif.gz')!;
    for (const f of [StructureFile, ChainTableFile]) {
      const r = rows.find((x) => `fixtures/${x.file}` === f)!;
      const b = readFileSync(`${PackageFolder}${f}`);
      expect([r.bytesFixture, r.sha256Fixture, r.sha256Download]).toEqual([
        b.byteLength,
        createHash('sha256').update(b).digest('hex'),
        entry.sha256Inflated,
      ]);
    }
  });

  // the whole entry parsed: slow by nature, the medium lane (the committed files above are the quick representative)
  it('the record rules write the committed files again from the mmCIF; its metal links are the view’s', {
    tags: ['medium'],
  }, async () => {
    // the parser imported here: only this medium test needs it (the quick lane does not load @molstar/io)
    const { readCif } = await import('../../scripts/cif.js');
    const cif = await readCif(gunzipSync(readFileSync(`${PackageFolder}fixtures/structures/2OCJ.cif.gz`)));
    expect(structureSeq(cif)).toBe(text(StructureFile));
    expect(chainTableJson(chainTable(cif, '2OCJ'))).toBe(text(ChainTableFile));
    const metal = cif
      .rows('struct_conn')
      .filter((r) => r.conn_type_id === 'metalc')
      .map((r) => [r.ptnr1_auth_asym_id, Number(r.ptnr1_auth_seq_id), r.ptnr1_label_atom_id, r.ptnr2_label_asym_id]);
    expect(metal).toEqual(Chains.flatMap(({ chain, zinc }) => ZincSites.map((s) => [chain, s.auth, s.atom, zinc])));
  });
});
