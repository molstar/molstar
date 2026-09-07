/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ryan DiRisio <rjdiris@gmail.com>
 */

import { TokenBuilder, Tokenizer } from '../../../mol-io/reader/common/text/tokenizer';
import { guessElementSymbolTokens } from '../util';
import { pdbToMmCif } from '../pdb/to-cif';
import { PdbFile } from '../../../mol-io/reader/pdb/schema';

/** Helper: build a PdbFile from a raw PDB string. */
function makePdb(pdbText: string): PdbFile {
    const lines = Tokenizer.readAllLines(pdbText);
    return { lines, variant: 'pdb' };
}

const records = [
    ['ATOM     19 HD23 LEU A   1     151.940 143.340 155.670  0.00  0.00', 'H'],
    ['ATOM     38  CA  SER A   3     146.430 138.150 162.270  0.00  0.00', 'C'],
    ['ATOM     38 NA   SER A   3     146.430 138.150 162.270  0.00  0.00', 'NA'],
    ['ATOM     38  NAA SER A   3     146.430 138.150 162.270  0.00  0.00', 'N'],
];

describe('PDB to-cif', () => {
    it('guess-element-symbol', () => {
        for (let i = 0, il = records.length; i < il; ++i) {
            const [data, element] = records[i];
            const tokens = TokenBuilder.create(data, 2);
            guessElementSymbolTokens(tokens, data, 12, 16);
            expect(data.substring(tokens.indices[0], tokens.indices[1])).toBe(element);
        }
    });
});

describe('PDB SEQRES-to-label_seq_id alignment', () => {
    it('assigns label_seq_id using SEQRES alignment for a complete chain', async () => {
        // SEQRES declares ALA GLY VAL; ATOM records observe ALA GLY VAL at auth_seq_id 1-3.
        // All residues match, so label_seq_id should be 1, 2, 3.
        const pdb = makePdb([
            'SEQRES   1 A    3  ALA GLY VAL                                               ',
            'ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  ',
            'ATOM      2  CA  GLY A   2       4.000   5.000   6.000  1.00  0.00           C  ',
            'ATOM      3  CA  VAL A   3       7.000   8.000   9.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelSeqId = atomSite.getField('label_seq_id')!;

        expect(labelSeqId.int(0)).toBe(1);
        expect(labelSeqId.int(1)).toBe(2);
        expect(labelSeqId.int(2)).toBe(3);
    });

    it('offsets label_seq_id when leading SEQRES residues are missing from ATOM', async () => {
        // SEQRES: MET ALA GLY VAL (4 residues)
        // ATOM: only ALA GLY VAL at auth_seq_id 2-4.
        // MET is unobserved, so ALA aligns to SEQRES position 2, GLY→3, VAL→4.
        const pdb = makePdb([
            'SEQRES   1 A    4  MET ALA GLY VAL                                           ',
            'ATOM      1  CA  ALA A   2       1.000   2.000   3.000  1.00  0.00           C  ',
            'ATOM      2  CA  GLY A   3       4.000   5.000   6.000  1.00  0.00           C  ',
            'ATOM      3  CA  VAL A   4       7.000   8.000   9.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelSeqId = atomSite.getField('label_seq_id')!;

        expect(labelSeqId.int(0)).toBe(2); // ALA -> SEQRES pos 2
        expect(labelSeqId.int(1)).toBe(3); // GLY -> SEQRES pos 3
        expect(labelSeqId.int(2)).toBe(4); // VAL -> SEQRES pos 4
    });

    it('handles an internal gap in observed residues', async () => {
        // SEQRES: ALA GLY VAL LEU (4 residues)
        // ATOM: ALA VAL LEU (GLY is missing in the middle)
        // Alignment: ALA→1, VAL→3, LEU→4
        const pdb = makePdb([
            'SEQRES   1 A    4  ALA GLY VAL LEU                                           ',
            'ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  ',
            'ATOM      2  CA  VAL A   3       4.000   5.000   6.000  1.00  0.00           C  ',
            'ATOM      3  CA  LEU A   4       7.000   8.000   9.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelSeqId = atomSite.getField('label_seq_id')!;

        expect(labelSeqId.int(0)).toBe(1); // ALA -> SEQRES pos 1
        expect(labelSeqId.int(1)).toBe(3); // VAL -> SEQRES pos 3
        expect(labelSeqId.int(2)).toBe(4); // LEU -> SEQRES pos 4
    });

    it('assigns label_seq_id for multiple chains independently', async () => {
        // Chain A SEQRES: MET ALA, ATOM: ALA (at pos 2)
        // Chain B SEQRES: GLY VAL LEU, ATOM: GLY VAL LEU (at pos 1,2,3)
        const pdb = makePdb([
            'SEQRES   1 A    2  MET ALA                                                   ',
            'SEQRES   1 B    3  GLY VAL LEU                                                ',
            'ATOM      1  CA  ALA A   2       1.000   2.000   3.000  1.00  0.00           C  ',
            'TER       2      ALA A   2                                                     ',
            'ATOM      3  CA  GLY B   1       4.000   5.000   6.000  1.00  0.00           C  ',
            'ATOM      4  CA  VAL B   2       7.000   8.000   9.000  1.00  0.00           C  ',
            'ATOM      5  CA  LEU B   3      10.000  11.000  12.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelSeqId = atomSite.getField('label_seq_id')!;

        // Chain A: ALA aligns to SEQRES position 2
        expect(labelSeqId.int(0)).toBe(2);
        // Chain B: GLY→1, VAL→2, LEU→3
        expect(labelSeqId.int(1)).toBe(1);
        expect(labelSeqId.int(2)).toBe(2);
        expect(labelSeqId.int(3)).toBe(3);
    });

    it('populates pdbx_unobs_or_zero_occ_residues for missing SEQRES residues', async () => {
        // SEQRES: MET ALA GLY VAL (4 residues), ATOM: ALA GLY (pos 2,3)
        // MET (pos 1) and VAL (pos 4) should appear in unobs
        const pdb = makePdb([
            'SEQRES   1 A    4  MET ALA GLY VAL                                           ',
            'ATOM      1  N   ALA A   2       1.000   2.000   3.000  1.00  0.00           N  ',
            'ATOM      2  CA  ALA A   2       1.500   2.500   3.500  1.00  0.00           C  ',
            'ATOM      3  C   ALA A   2       2.000   3.000   4.000  1.00  0.00           C  ',
            'ATOM      4  N   GLY A   3       4.000   5.000   6.000  1.00  0.00           N  ',
            'ATOM      5  CA  GLY A   3       4.500   5.500   6.500  1.00  0.00           C  ',
            'ATOM      6  C   GLY A   3       5.000   6.000   7.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const unobs = cif.categories['pdbx_unobs_or_zero_occ_residues'];
        expect(unobs).toBeDefined();

        const unobsSeqId = unobs.getField('label_seq_id')!;
        const unobsCompId = unobs.getField('label_comp_id')!;
        expect(unobsSeqId.rowCount).toBe(2);

        // Positions 1 (MET) and 4 (VAL) are unobserved
        expect(unobsSeqId.int(0)).toBe(1);
        expect(unobsCompId.str(0)).toBe('MET');
        expect(unobsSeqId.int(1)).toBe(4);
        expect(unobsCompId.str(1)).toBe('VAL');

        // entity_poly_seq should list all 4 SEQRES residues with 1-based numbering
        const entityPolySeq = cif.categories['entity_poly_seq'];
        expect(entityPolySeq).toBeDefined();
        const epsNum = entityPolySeq.getField('num')!;
        const epsMonId = entityPolySeq.getField('mon_id')!;
        expect(epsNum.rowCount).toBe(4);
        expect(epsNum.int(0)).toBe(1);
        expect(epsMonId.str(0)).toBe('MET');
        expect(epsNum.int(1)).toBe(2);
        expect(epsMonId.str(1)).toBe('ALA');
        expect(epsNum.int(2)).toBe(3);
        expect(epsMonId.str(2)).toBe('GLY');
        expect(epsNum.int(3)).toBe(4);
        expect(epsMonId.str(3)).toBe('VAL');
    });

    it('handles multiple atoms per residue correctly', async () => {
        // SEQRES: MET ALA, ATOM: ALA with N, CA, C atoms
        // All 3 atoms of ALA should get label_seq_id = 2
        const pdb = makePdb([
            'SEQRES   1 A    2  MET ALA                                                   ',
            'ATOM      1  N   ALA A   2       1.000   2.000   3.000  1.00  0.00           N  ',
            'ATOM      2  CA  ALA A   2       1.500   2.500   3.500  1.00  0.00           C  ',
            'ATOM      3  C   ALA A   2       2.000   3.000   4.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelSeqId = atomSite.getField('label_seq_id')!;

        expect(labelSeqId.int(0)).toBe(2);
        expect(labelSeqId.int(1)).toBe(2);
        expect(labelSeqId.int(2)).toBe(2);
    });

    it('produces linear label_seq_id with overlapping auth_seq_id and insertion codes (1NSA-like, regression #1730)', async () => {
        // Regression test for https://github.com/molstar/molstar/issues/1730
        // 1NSA has residues numbered 7A, 8A, ... (with insertion code) followed by
        // 4, 5, 6, 7, 8, ... (without insertion code). The overlapping auth_seq_id
        // values (e.g. 7A vs 7) previously caused interleaved label_seq_id, breaking
        // residue ordering. The fix ensures label_seq_id is strictly sequential when
        // insertion codes are present, even without SEQRES.
        const pdb = makePdb([
            'ATOM      1  CA  HIS A   7A      1.000   2.000   3.000  1.00  0.00           C  ',
            'ATOM      2  CA  PHE A   8A      4.000   5.000   6.000  1.00  0.00           C  ',
            'ATOM      3  CA  GLY A   9A      7.000   8.000   9.000  1.00  0.00           C  ',
            'ATOM      4  CA  ALA A   4      10.000  11.000  12.000  1.00  0.00           C  ',
            'ATOM      5  CA  VAL A   5      13.000  14.000  15.000  1.00  0.00           C  ',
            'ATOM      6  CA  LEU A   6      16.000  17.000  18.000  1.00  0.00           C  ',
            'ATOM      7  CA  ILE A   7      19.000  20.000  21.000  1.00  0.00           C  ',
            'ATOM      8  CA  PRO A   8      22.000  23.000  24.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelSeqId = atomSite.getField('label_seq_id')!;

        // label_seq_id must be strictly increasing: 1, 2, 3, 4, 5, 6, 7, 8
        for (let i = 0; i < labelSeqId.rowCount; ++i) {
            expect(labelSeqId.int(i)).toBe(i + 1);
        }
    });

    it('produces linear label_seq_id with insertion codes when SEQRES is present', async () => {
        // Same 1NSA-like scenario but now with a SEQRES record.
        // The alignment should still produce strictly sequential label_seq_id.
        const pdb = makePdb([
            'SEQRES   1 A    8  HIS PHE GLY ALA VAL LEU ILE PRO                           ',
            'ATOM      1  CA  HIS A   7A      1.000   2.000   3.000  1.00  0.00           C  ',
            'ATOM      2  CA  PHE A   8A      4.000   5.000   6.000  1.00  0.00           C  ',
            'ATOM      3  CA  GLY A   9A      7.000   8.000   9.000  1.00  0.00           C  ',
            'ATOM      4  CA  ALA A   4      10.000  11.000  12.000  1.00  0.00           C  ',
            'ATOM      5  CA  VAL A   5      13.000  14.000  15.000  1.00  0.00           C  ',
            'ATOM      6  CA  LEU A   6      16.000  17.000  18.000  1.00  0.00           C  ',
            'ATOM      7  CA  ILE A   7      19.000  20.000  21.000  1.00  0.00           C  ',
            'ATOM      8  CA  PRO A   8      22.000  23.000  24.000  1.00  0.00           C  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelSeqId = atomSite.getField('label_seq_id')!;

        // All 8 residues observed, aligned to SEQRES positions 1–8
        for (let i = 0; i < labelSeqId.rowCount; ++i) {
            expect(labelSeqId.int(i)).toBe(i + 1);
        }
    });
});

describe('PDB label_atom_id unicity', () => {
    it('preserves unique atom names within a residue', async () => {
        const pdb = makePdb([
            'ATOM   2161  CG2 MOL A 303      31.377  65.493  27.278  1.00 32.12           C  ',
            'ATOM   2162  OG2 MOL A 303      31.727  64.432  26.692  1.00 30.76           O  ',
            'ATOM   2163  OG3 MOL A 303      30.170  65.084  27.156  1.00 32.03           O  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelAtomId = atomSite.getField('label_atom_id')!;

        expect(labelAtomId.str(0)).toBe('CG2');
        expect(labelAtomId.str(1)).toBe('OG2');
        expect(labelAtomId.str(2)).toBe('OG3');
    });

    it('deduplicates repeated atom names within a residue', async () => {
        const pdb = makePdb([
            'ATOM   2161  C   MOL A 303      31.377  65.493  27.278  1.00 32.12           C  ',
            'ATOM   2162  O   MOL A 303      31.727  64.432  26.692  1.00 30.76           O  ',
            'ATOM   2163  O   MOL A 303      30.170  65.084  27.156  1.00 32.03           O  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelAtomId = atomSite.getField('label_atom_id')!;

        expect(labelAtomId.str(0)).toBe('C');
        expect(labelAtomId.str(1)).toBe('O');
        expect(labelAtomId.str(2)).toBe('O_1'); // deduplicated
    });

    it('No suffix added to same atom names for altlocs)', async () => {
        const pdb = makePdb([
            'HETATM 2161  CG2 SRY A 303      31.377  65.493  27.278  1.00 32.12           C  ',
            'HETATM 2162  OG2ASRY A 303      31.727  64.432  26.692  0.50 30.76           O  ',
            'HETATM 2163  OG2BSRY A 303      30.170  65.084  27.156  0.50 32.03           O  ',
            'END                                                                             ',
        ].join('\n'));

        const cif = await pdbToMmCif(pdb);
        const atomSite = cif.categories['atom_site'];
        const labelAtomId = atomSite.getField('label_atom_id')!;
        const labelAltId = atomSite.getField('label_alt_id')!;

        expect(labelAtomId.str(0)).toBe('CG2');
        expect(labelAtomId.str(1)).toBe('OG2');
        expect(labelAtomId.str(2)).toBe('OG2');
        expect(labelAltId.str(1)).toBe('A');
        expect(labelAltId.str(2)).toBe('B');
    });
});