/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Kim Juho <juho_kim@outlook.com>
 */

import { CifCategory, CifField } from '../../../mol-io/reader/cif';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { Tokens } from '../../../mol-io/reader/common/text/tokenizer';
import { Column } from '../../../mol-data/db';

const HelixTypes: {[k: string]: mmCIF_Schema['struct_conf']['conf_type_id']['T']} = {
    // CLASS NUMBER
    // TYPE OF  HELIX                     (COLUMNS 39 - 40)
    // --------------------------------------------------------------
    // Right-handed alpha (default)                1
    // Right-handed omega                          2
    // Right-handed pi                             3
    // Right-handed gamma                          4
    // Right-handed 3 - 10                         5
    // Left-handed alpha                           6
    // Left-handed omega                           7
    // Left-handed gamma                           8
    // 2 - 7 ribbon/helix                          9
    // Polyproline                                10
    1: 'helx_rh_al_p',
    2: 'helx_rh_om_p',
    3: 'helx_rh_pi_p',
    4: 'helx_rh_ga_p',
    5: 'helx_rh_3t_p',
    6: 'helx_lh_al_p',
    7: 'helx_lh_om_p',
    8: 'helx_lh_ga_p',
    9: 'helx_rh_27_p', // TODO or left-handed???
    10: 'helx_rh_pp_p', // TODO or left-handed???
};
function getStructConfTypeId(type: string): mmCIF_Schema['struct_conf']['conf_type_id']['T'] {
    return HelixTypes[type] || 'helx_p';
}

interface PdbHelix {
    serNum: string,
    helixID: string,
    initResName: string,
    initChainID: string,
    initSeqNum: string,
    initICode: string,
    endResName: string,
    endChainID: string,
    endSeqNum: string,
    endICode: string,
    helixClass: string,
    comment: string,
    length: string
}

export function parseHelix(lines: Tokens, lineStart: number, lineEnd: number): CifCategory {
    const helices: PdbHelix[] = [];
    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);

    for (let i = lineStart; i < lineEnd; i++) {
        const line = getLine(i);
        // COLUMNS        DATA  TYPE     FIELD         DEFINITION
        // -----------------------------------------------------------------------------------
        //  1 -  6        Record name    "HELIX "
        //  8 - 10        Integer        serNum        Serial number of the helix. This starts
        //                                             at 1  and increases incrementally.
        // 12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
        //                                             number, each helix is given an
        //                                             alphanumeric character helix identifier.
        // 16 - 19        Residue name   initResName   Name of the initial residue.
        //                                             PDB spec defines 3-letter for residue name,
        //                                             but 4-letter are commonly used
        // 20             Character      initChainID   Chain identifier for the chain containing
        //                                             this  helix.
        // 22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
        // 26             AChar          initICode     Insertion code of the initial residue.
        // 28 - 31        Residue name   endResName    Name of the terminal residue of the helix.
        //                                             PDB spec defines 3-letter for residue name,
        //                                             but 4-letter are commonly used
        // 32             Character      endChainID    Chain identifier for the chain containing
        //                                             this  helix.
        // 34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
        // 38             AChar          endICode      Insertion code of the terminal residue.
        // 39 - 40        Integer        helixClass    Helix class (see below).
        // 41 - 70        String         comment       Comment about this helix.
        // 72 - 76        Integer        length        Length of this helix.
        helices.push({
            serNum: line.substring(7, 10).trim(),
            helixID: line.substring(11, 14).trim(),
            initResName: line.substring(15, 19).trim(),
            initChainID: line.substring(19, 20).trim(),
            initSeqNum: line.substring(21, 25).trim(),
            initICode: line.substring(25, 26).trim(),
            endResName: line.substring(27, 31).trim(),
            endChainID: line.substring(31, 34).trim(),
            endSeqNum: line.substring(33, 37).trim(),
            endICode: line.substring(37, 38).trim(),
            helixClass: line.substring(38, 40).trim(),
            comment: line.substring(40, 70).trim(),
            length: line.substring(71, 76).trim()
        });
    }

    const beg_auth_asym_id = CifField.ofStrings(helices.map(h => h.initChainID));
    const beg_auth_comp_id = CifField.ofStrings(helices.map(h => h.initResName));

    const end_auth_asym_id = CifField.ofStrings(helices.map(h => h.endChainID));
    const end_auth_comp_id = CifField.ofStrings(helices.map(h => h.endResName));

    const struct_conf: CifCategory.Fields<mmCIF_Schema['struct_conf']> = {
        beg_label_asym_id: beg_auth_asym_id,
        beg_label_comp_id: beg_auth_comp_id,
        beg_label_seq_id: CifField.ofUndefined(helices.length, Column.Schema.int),
        beg_auth_asym_id,
        beg_auth_comp_id,
        beg_auth_seq_id: CifField.ofStrings(helices.map(h => h.initSeqNum)),

        conf_type_id: CifField.ofStrings(helices.map(h => getStructConfTypeId(h.helixClass))),
        details: CifField.ofStrings(helices.map(h => h.comment)),

        end_label_asym_id: end_auth_asym_id,
        end_label_comp_id: end_auth_comp_id,
        end_label_seq_id: CifField.ofUndefined(helices.length, Column.Schema.int),
        end_auth_asym_id,
        end_auth_comp_id,
        end_auth_seq_id: CifField.ofStrings(helices.map(h => h.endSeqNum)),

        id: CifField.ofStrings(helices.map(h => h.serNum)),
        pdbx_beg_PDB_ins_code: CifField.ofStrings(helices.map(h => h.initICode)),
        pdbx_end_PDB_ins_code: CifField.ofStrings(helices.map(h => h.endICode)),
        pdbx_PDB_helix_class: CifField.ofStrings(helices.map(h => h.helixClass)),
        pdbx_PDB_helix_length: CifField.ofStrings(helices.map(h => h.length)),
        pdbx_PDB_helix_id: CifField.ofStrings(helices.map(h => h.helixID)),
    };
    return CifCategory.ofFields('struct_conf', struct_conf);
}

//

interface PdbSheet {
    strand: string,
    sheetID: string,
    numStrands: string,
    initResName: string,
    initChainID: string,
    initSeqNum: string,
    initICode: string,
    endResName: string,
    endChainID: string,
    endSeqNum: string,
    endICode: string,
    sense: string,
    curAtom: string,
    curResName: string,
    curChainId: string,
    curResSeq: string,
    curICode: string,
    prevAtom: string,
    prevResName: string,
    prevChainId: string,
    prevResSeq: string,
    prevICode: string,
}

export function parseSheet(lines: Tokens, lineStart: number, lineEnd: number): CifCategory {
    const sheets: PdbSheet[] = [];
    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);

    for (let i = lineStart; i < lineEnd; i++) {
        const line = getLine(i);
        // COLUMNS       DATA  TYPE     FIELD          DEFINITION
        // -------------------------------------------------------------------------------------
        //  1 -  6        Record name   "SHEET "
        //  8 - 10        Integer       strand         Strand  number which starts at 1 for each
        //                                             strand within a sheet and increases by one.
        // 12 - 14        LString(3)    sheetID        Sheet  identifier.
        // 15 - 16        Integer       numStrands     Number  of strands in sheet.
        // 18 - 21        Residue name  initResName    Residue  name of initial residue.
        //                                             PDB spec defines 3-letter for residue name,
        //                                             but 4-letter are commonly used
        // 22             Character     initChainID    Chain identifier of initial residue
        //                                             in strand.
        // 23 - 26        Integer       initSeqNum     Sequence number of initial residue
        //                                             in strand.
        // 27             AChar         initICode      Insertion code of initial residue
        //                                             in  strand.
        // 29 - 32        Residue name  endResName     Residue name of terminal residue.
        //                                             PDB spec defines 3-letter for residue name,
        //                                             but 4-letter are commonly used
        // 33             Character     endChainID     Chain identifier of terminal residue.
        // 34 - 37        Integer       endSeqNum      Sequence number of terminal residue.
        // 38             AChar         endICode       Insertion code of terminal residue.
        // 39 - 40        Integer       sense          Sense of strand with respect to previous
        //                                             strand in the sheet. 0 if first strand,
        //                                             1 if  parallel,and -1 if anti-parallel.
        // 42 - 45        Atom          curAtom        Registration.  Atom name in current strand.
        // 46 - 49        Residue name  curResName     Registration.  Residue name in current strand
        //                                             PDB spec defines 3-letter for residue name,
        //                                             but 4-letter are commonly used
        // 50             Character     curChainId     Registration. Chain identifier in
        //                                             current strand.
        // 51 - 54        Integer       curResSeq      Registration.  Residue sequence number
        //                                             in current strand.
        // 55             AChar         curICode       Registration. Insertion code in
        //                                             current strand.
        // 57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
        // 61 - 64        Residue name  prevResName    Registration.  Residue name in
        //                                             previous strand.
        //                                             PDB spec defines 3-letter for residue name,
        //                                             but 4-letter are commonly used
        // 65             Character     prevChainId    Registration.  Chain identifier in
        //                                             previous  strand.
        // 66 - 69        Integer       prevResSeq     Registration. Residue sequence number
        //                                             in previous strand.
        // 70             AChar         prevICode      Registration.  Insertion code in
        //                                             previous strand.
        sheets.push({
            strand: line.substring(7, 10).trim(),
            sheetID: line.substring(11, 14).trim(),
            numStrands: line.substring(14, 16).trim(),
            initResName: line.substring(17, 21).trim(),
            initChainID: line.substring(21, 22).trim(),
            initSeqNum: line.substring(22, 26).trim(),
            initICode: line.substring(26, 27).trim(),
            endResName: line.substring(28, 32).trim(),
            endChainID: line.substring(32, 33).trim(),
            endSeqNum: line.substring(33, 37).trim(),
            endICode: line.substring(37, 38).trim(),
            sense: line.substring(38, 40).trim(),
            curAtom: line.substring(41, 45).trim(),
            curResName: line.substring(45, 49).trim(),
            curChainId: line.substring(49, 50).trim(),
            curResSeq: line.substring(50, 54).trim(),
            curICode: line.substring(54, 55).trim(),
            prevAtom: line.substring(56, 60).trim(),
            prevResName: line.substring(60, 64).trim(),
            prevChainId: line.substring(64, 65).trim(),
            prevResSeq: line.substring(65, 69).trim(),
            prevICode: line.substring(69, 70).trim(),
        });
    }

    const beg_auth_asym_id = CifField.ofStrings(sheets.map(s => s.initChainID));
    const beg_auth_comp_id = CifField.ofStrings(sheets.map(s => s.initResName));
    const beg_auth_seq_id = CifField.ofStrings(sheets.map(s => s.initSeqNum));

    const end_auth_asym_id = CifField.ofStrings(sheets.map(s => s.endChainID));
    const end_auth_comp_id = CifField.ofStrings(sheets.map(s => s.endResName));
    const end_auth_seq_id = CifField.ofStrings(sheets.map(s => s.endSeqNum));

    const struct_sheet_range: CifCategory.Fields<mmCIF_Schema['struct_sheet_range']> = {
        beg_label_asym_id: beg_auth_asym_id,
        beg_label_comp_id: beg_auth_comp_id,
        beg_label_seq_id: beg_auth_seq_id,
        beg_auth_asym_id,
        beg_auth_comp_id,
        beg_auth_seq_id,

        end_label_asym_id: end_auth_asym_id,
        end_label_comp_id: end_auth_asym_id,
        end_label_seq_id: end_auth_seq_id,
        end_auth_asym_id,
        end_auth_comp_id,
        end_auth_seq_id,

        id: CifField.ofStrings(sheets.map(s => s.strand)),
        sheet_id: CifField.ofStrings(sheets.map(s => s.sheetID)), // TODO wrong, needs to point to _struct_sheet.id
        pdbx_beg_PDB_ins_code: CifField.ofStrings(sheets.map(s => s.initICode)),
        pdbx_end_PDB_ins_code: CifField.ofStrings(sheets.map(s => s.endICode)),
    };
    return CifCategory.ofFields('struct_sheet_range', struct_sheet_range);
}