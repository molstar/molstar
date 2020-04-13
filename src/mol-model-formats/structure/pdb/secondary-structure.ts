/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifCategory, CifField } from '../../../mol-io/reader/cif';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { Tokens } from '../../../mol-io/reader/common/text/tokenizer';

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
    1: 'HELX_RH_AL_P',
    2: 'HELX_RH_OM_P',
    3: 'HELX_RH_PI_P',
    4: 'HELX_RH_GA_P',
    5: 'HELX_RH_3T_P',
    6: 'HELX_LH_AL_P',
    7: 'HELX_LH_OM_P',
    8: 'HELX_LH_GA_P',
    9: 'HELX_RH_27_P', // TODO or left-handed???
    10: 'HELX_RH_PP_P', // TODO or left-handed???
};
function getStructConfTypeId(type: string): mmCIF_Schema['struct_conf']['conf_type_id']['T'] {
    return HelixTypes[type] || 'HELX_P';
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
        let line = getLine(i);
        // COLUMNS        DATA  TYPE     FIELD         DEFINITION
        // -----------------------------------------------------------------------------------
        // 1 -  6        Record name    "HELIX "
        // 8 - 10        Integer        serNum        Serial number of the helix. This starts
        //                                             at 1  and increases incrementally.
        // 12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
        //                                             number, each helix is given an
        //                                             alphanumeric character helix identifier.
        // 16 - 18        Residue name   initResName   Name of the initial residue.
        // 20             Character      initChainID   Chain identifier for the chain containing
        //                                             this  helix.
        // 22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
        // 26             AChar          initICode     Insertion code of the initial residue.
        // 28 - 30        Residue  name  endResName    Name of the terminal residue of the helix.
        // 32             Character      endChainID    Chain identifier for the chain containing
        //                                             this  helix.
        // 34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
        // 38             AChar          endICode      Insertion code of the terminal residue.
        // 39 - 40        Integer        helixClass    Helix class (see below).
        // 41 - 70        String         comment       Comment about this helix.
        // 72 - 76        Integer        length        Length of this helix.
        helices.push({
            serNum: line.substr(7, 3).trim(),
            helixID: line.substr(11, 3).trim(),
            initResName: line.substr(15, 3).trim(),
            initChainID: line.substr(19, 1).trim(),
            initSeqNum: line.substr(21, 4).trim(),
            initICode: line.substr(25, 1).trim(),
            endResName: line.substr(27, 3).trim(),
            endChainID: line.substr(31, 3).trim(),
            endSeqNum: line.substr(33, 4).trim(),
            endICode: line.substr(37, 1).trim(),
            helixClass: line.substr(38, 2).trim(),
            comment: line.substr(40, 30).trim(),
            length: line.substr(71, 5).trim()
        });
    }

    const beg_auth_asym_id = CifField.ofStrings(helices.map(h => h.initChainID));
    const beg_auth_comp_id = CifField.ofStrings(helices.map(h => h.initResName));
    const beg_auth_seq_id = CifField.ofStrings(helices.map(h => h.initSeqNum));

    const end_auth_asym_id = CifField.ofStrings(helices.map(h => h.endChainID));
    const end_auth_comp_id = CifField.ofStrings(helices.map(h => h.endResName));
    const end_auth_seq_id = CifField.ofStrings(helices.map(h => h.endSeqNum));

    const struct_conf: CifCategory.Fields<mmCIF_Schema['struct_conf']> = {
        beg_label_asym_id: beg_auth_asym_id,
        beg_label_comp_id: beg_auth_comp_id,
        beg_label_seq_id: beg_auth_seq_id,
        beg_auth_asym_id,
        beg_auth_comp_id,
        beg_auth_seq_id,

        conf_type_id: CifField.ofStrings(helices.map(h => getStructConfTypeId(h.helixClass))),
        details: CifField.ofStrings(helices.map(h => h.comment)),

        end_label_asym_id: end_auth_asym_id,
        end_label_comp_id: end_auth_asym_id,
        end_label_seq_id: end_auth_seq_id,
        end_auth_asym_id,
        end_auth_comp_id,
        end_auth_seq_id,

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
        let line = getLine(i);
        // COLUMNS       DATA  TYPE     FIELD          DEFINITION
        // -------------------------------------------------------------------------------------
        // 1 -  6        Record name   "SHEET "
        // 8 - 10        Integer       strand         Strand  number which starts at 1 for each
        //                                             strand within a sheet and increases by one.
        // 12 - 14        LString(3)    sheetID        Sheet  identifier.
        // 15 - 16        Integer       numStrands     Number  of strands in sheet.
        // 18 - 20        Residue name  initResName    Residue  name of initial residue.
        // 22             Character     initChainID    Chain identifier of initial residue
        //                                             in strand.
        // 23 - 26        Integer       initSeqNum     Sequence number of initial residue
        //                                             in strand.
        // 27             AChar         initICode      Insertion code of initial residue
        //                                             in  strand.
        // 29 - 31        Residue name  endResName     Residue name of terminal residue.
        // 33             Character     endChainID     Chain identifier of terminal residue.
        // 34 - 37        Integer       endSeqNum      Sequence number of terminal residue.
        // 38             AChar         endICode       Insertion code of terminal residue.
        // 39 - 40        Integer       sense          Sense of strand with respect to previous
        //                                             strand in the sheet. 0 if first strand,
        //                                             1 if  parallel,and -1 if anti-parallel.
        // 42 - 45        Atom          curAtom        Registration.  Atom name in current strand.
        // 46 - 48        Residue name  curResName     Registration.  Residue name in current strand
        // 50             Character     curChainId     Registration. Chain identifier in
        //                                             current strand.
        // 51 - 54        Integer       curResSeq      Registration.  Residue sequence number
        //                                             in current strand.
        // 55             AChar         curICode       Registration. Insertion code in
        //                                             current strand.
        // 57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
        // 61 - 63        Residue name  prevResName    Registration.  Residue name in
        //                                             previous strand.
        // 65             Character     prevChainId    Registration.  Chain identifier in
        //                                             previous  strand.
        // 66 - 69        Integer       prevResSeq     Registration. Residue sequence number
        //                                             in previous strand.
        // 70             AChar         prevICode      Registration.  Insertion code in
        //                                             previous strand.
        sheets.push({
            strand: line.substr(7, 3).trim(),
            sheetID: line.substr(11, 3).trim(),
            numStrands: line.substr(14, 2).trim(),
            initResName: line.substr(17, 3).trim(),
            initChainID: line.substr(21, 1).trim(),
            initSeqNum: line.substr(22, 4).trim(),
            initICode: line.substr(26, 1).trim(),
            endResName: line.substr(28, 3).trim(),
            endChainID: line.substr(32, 1).trim(),
            endSeqNum: line.substr(33, 4).trim(),
            endICode: line.substr(37, 1).trim(),
            sense: line.substr(38, 2).trim(),
            curAtom: line.substr(41, 4).trim(),
            curResName: line.substr(45, 3).trim(),
            curChainId: line.substr(49, 1).trim(),
            curResSeq: line.substr(50, 4).trim(),
            curICode: line.substr(54, 1).trim(),
            prevAtom: line.substr(56, 4).trim(),
            prevResName: line.substr(60, 3).trim(),
            prevChainId: line.substr(64, 1).trim(),
            prevResSeq: line.substr(65, 4).trim(),
            prevICode: line.substr(69, 1).trim(),
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