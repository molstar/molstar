/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { decodeSeq, encodeSeq } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { HemoglobinFile, SeedFile } from '../globins/data.js';
import {
  cifChains,
  cifHelices,
  columnsFromUniprot,
  cutSeed,
  hemeContacts,
  hemoglobinSeq,
  SeedRows,
  seedVersion,
} from '../globins/recording.js';
import { PackageFolder } from '../node-assets.js';

// The record script's rules for the Globins (`scripts/record.ts globins`), held without the network: the seed's cut
// over the committed fixture, the columns file from a UniProt TSV, and 4HHB's chains, helices and heme contacts over a
// hand-written mmCIF.

const text = (f: string) => readFileSync(`${PackageFolder}${f}`, 'utf8');

const letters = (system: string) => {
  const p = decodeSeq(text(HemoglobinFile)).data!.pairs.find((x) => x.request.system === system)!.answer;
  return p.ok ? (p.tables.letters as { letters: string }).letters : '';
};

/** A mini mmCIF: chain A (Val-His-Leu) and chain B (Gly-His), a helix each, a heme on A beside A's His. */
const MiniCif = `data_MINI
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.auth_seq_num
A VAL 1
A HIS 2
A LEU 3
B GLY 1
B HIS 2
#
loop_
_struct_conf.id
_struct_conf.pdbx_PDB_helix_id
_struct_conf.beg_label_asym_id
_struct_conf.beg_label_seq_id
_struct_conf.end_label_seq_id
HELX_P1 AA A 1 3
HELX_P2 BA B 1 2
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
ATOM 1 N N VAL A 1 10.0 0.0 0.0 A
ATOM 2 N NE2 HIS A 2 2.0 0.0 0.0 A
ATOM 3 H H HIS A 2 0.5 0.0 0.0 A
ATOM 4 C CA LEU A 3 0.0 9.0 0.0 A
ATOM 5 N NE2 HIS B 2 1.0 0.0 0.0 B
HETATM 6 FE FE HEM E . 0.0 0.0 0.0 A
HETATM 7 C CHA HEM E . 0.0 0.0 3.5 A
#
`;

describe('the Globins recorder’s rules', () => {
  it('the seed’s cut over the committed fixture: already 60 rows, it comes back as it is', () => {
    const seed = text(SeedFile);
    const cut = cutSeed(seed, [letters('hba'), letters('hbb')]);
    expect(cut.removed).toEqual([]);
    expect(cut.kept.length).toBe(SeedRows);
    expect(cut.text).toBe(seed);
    expect(seed).toMatch(/^#=GF SQ\s+60$/m);
    expect(seedVersion(seed)).toBe('PF00042.29');
  });

  it('a cut says the rows kept in #=GF SQ and drops the rows removed with their #=GS and #=GR lines', () => {
    const rows = Array.from({ length: SeedRows + 1 }, (_, k) => `R${k}/1-5`);
    const seed = [
      '# STOCKHOLM 1.0',
      '#=GF AC   PF00000.1',
      `#=GF SQ   ${rows.length}`,
      ...rows.map((r) => `#=GS ${r} AC P${k(r)}.1`),
      `#=GR ${rows[SeedRows]} SS ....`,
      ...rows.map((r, i) => `${r} ${i === SeedRows ? 'WWWWW' : 'VLSPA'}`),
      '//',
    ].join('\n');
    const cut = cutSeed(seed, ['VLSPA']);
    expect(cut.removed).toEqual([rows[SeedRows]]);
    expect(cut.text).toMatch(/^#=GF SQ\s+60$/m);
    expect(cut.text).not.toContain(rows[SeedRows]);
    expect(cut.text.split('\n').filter((l) => l.startsWith('#=GS')).length).toBe(SeedRows);
  });

  it('the columns file from UniProt’s TSV: keyed by the row whose AC names the entry, an organism without a class empty', () => {
    const seed = ['#=GS HBA1_BOSGR/26-136 AC P01966.1', '#=GS MYG_ALLMI/27-143 AC P02200.1'].join('\n');
    const tsv = [
      'Entry\tEntry Name\tOrganism\tLength\tTaxonomic lineage',
      'P01966\tHBA1_BOSGR\tBos grunniens (Wild yak)\t141\tEukaryota (domain), Mammalia (class), Bos (genus)',
      'P02200\tMYG_ALLMI\tAlligator mississippiensis\t154\tEukaryota (domain), Archosauria (clade)',
    ].join('\n');
    expect(columnsFromUniprot(tsv, seed)).toBe(
      'id\tOS\tclass\tuniprot_length\nHBA1_BOSGR/26-136\tBos grunniens (Wild yak)\tMammalia\t141\nMYG_ALLMI/27-143\tAlligator mississippiensis\t\t154\n',
    );
    expect(() => columnsFromUniprot(tsv.split('\n').slice(0, 2).join('\n'), seed)).toThrow(/P02200/);
  });

  it('a mini mmCIF: its chains, helices and heme contacts (heavy atoms within 4 Å, the heme’s chain), and its seq', () => {
    expect(cifChains(MiniCif)).toEqual([
      { asym: 'A', letters: 'VHL', authFirst: 1 },
      { asym: 'B', letters: 'GH', authFirst: 1 },
    ]);
    expect(cifHelices(MiniCif)).toEqual([
      { asym: 'A', name: 'AA', start: 1, end: 3 },
      { asym: 'B', name: 'BA', start: 1, end: 2 },
    ]);
    // A's His at 2 Å; its H (at 0.5 Å) not heavy, so Val and Leu out of reach stay out; B's His is not the heme's chain
    expect([...hemeContacts(MiniCif)]).toEqual([['A', [2]]]);
    const seq = hemoglobinSeq(MiniCif);
    expect(seq.coverage.systems).toEqual([
      { name: 'hba', length: 3, alphabet: 'protein' },
      { name: 'hbb', length: 2, alphabet: 'protein' },
      { name: 'hem', length: 1, alphabet: 'protein' },
    ]);
    const decoded = decodeSeq(encodeSeq(seq)).data!;
    const hba = decoded.pairs[0].answer;
    expect(hba.ok && (hba.tables.letters as { letters: string }).letters).toBe('VHL');
    const ranges = hba.ok ? (hba.tables.ranges as { start: Int32Array; end: Int32Array }) : undefined;
    expect([...ranges!.start, ...ranges!.end]).toEqual([0, 3]);
  });
});

function k(r: string) {
  return r.slice(1, r.indexOf('/'));
}
