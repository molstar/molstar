
/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { SecondaryStructureType } from '../../../mol-model/structure/model/types';
import { AtomicHierarchy } from '../../../mol-model/structure/model/properties/atomic';
import { SecondaryStructure } from '../../../mol-model/structure/model/properties/secondary-structure';
import { Column, Table } from '../../../mol-data/db';
import { ChainIndex, ResidueIndex } from '../../../mol-model/structure/model/indexing';
import { FormatPropertyProvider } from '../common/property';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';

export { ModelSecondaryStructure };

type StructConf = Table<mmCIF_Schema['struct_conf']>
type StructSheetRange = Table<mmCIF_Schema['struct_sheet_range']>
type CoordinateType = 'label' | 'auth';

namespace ModelSecondaryStructure {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'model_secondary_structure',
    };

    export const Provider = FormatPropertyProvider.create<SecondaryStructure>(Descriptor);

    export function fromStruct(conf: StructConf, sheetRange: StructSheetRange, hierarchy: AtomicHierarchy): SecondaryStructure {
        const map: SecondaryStructureMap = new Map();
        const elements: SecondaryStructure.Element[] = [{ kind: 'none' }];

        const coordinates = getCoordinateType(conf, sheetRange);

        addHelices(conf, coordinates, map, elements);
        // must add Helices 1st because of 'key' value assignment.
        addSheets(sheetRange, coordinates, map, conf._rowCount, elements);

        const n = hierarchy.residues._rowCount;
        const getIndex = (rI: ResidueIndex) => rI;

        const secStruct: SecondaryStructureData = {
            type: new Int32Array(n) as any,
            key: new Int32Array(n) as any,
            elements
        };

        if (map.size > 0) assignSecondaryStructureRanges(hierarchy, coordinates, map, secStruct);
        return SecondaryStructure(secStruct.type, secStruct.key, secStruct.elements, getIndex);
    }
}

function getCoordinateType(conf: StructConf, sheetRange: StructSheetRange): CoordinateType {
    if (conf._rowCount > 0) {
        if (conf.beg_label_seq_id.valueKind(0) !== Column.ValueKinds.Present || conf.end_label_seq_id.valueKind(0) !== Column.ValueKinds.Present) return 'auth';
    } else if (sheetRange) {
        if (sheetRange.beg_label_seq_id.valueKind(0) !== Column.ValueKinds.Present || sheetRange.end_label_seq_id.valueKind(0) !== Column.ValueKinds.Present) return 'auth';
    }
    return 'label';
}

type SecondaryStructureEntry = {
    startSeqId: number,
    startInsCode: string | null,
    endSeqId: number,
    endInsCode: string | null,
    type: SecondaryStructureType,
    key: number
}
type SecondaryStructureMap = Map<string, Map<number, SecondaryStructureEntry[]>>
type SecondaryStructureData = { type: SecondaryStructureType[], key: number[], elements: SecondaryStructure.Element[] }

function addHelices(cat: StructConf, coordinates: CoordinateType, map: SecondaryStructureMap, elements: SecondaryStructure.Element[]) {
    if (!cat._rowCount) return;

    const { beg_label_asym_id, beg_label_seq_id, beg_auth_seq_id, pdbx_beg_PDB_ins_code } = cat;
    const { end_label_seq_id, end_auth_seq_id, pdbx_end_PDB_ins_code } = cat;
    const { pdbx_PDB_helix_class, conf_type_id, details } = cat;

    const beg_seq_id = coordinates === 'label' ? beg_label_seq_id : beg_auth_seq_id;
    const end_seq_id = coordinates === 'label' ? end_label_seq_id : end_auth_seq_id;

    for (let i = 0, _i = cat._rowCount; i < _i; i++) {
        const type = SecondaryStructureType.create(pdbx_PDB_helix_class.valueKind(i) === Column.ValueKinds.Present
            ? SecondaryStructureType.SecondaryStructurePdb[pdbx_PDB_helix_class.value(i)]
            : conf_type_id.valueKind(i) === Column.ValueKinds.Present
                ? SecondaryStructureType.SecondaryStructureMmcif[conf_type_id.value(i)]
                : SecondaryStructureType.Flag.NA);

        const element: SecondaryStructure.Helix = {
            kind: 'helix',
            flags: type,
            type_id: conf_type_id.valueKind(i) === Column.ValueKinds.Present ? conf_type_id.value(i) : 'helx_p',
            helix_class: pdbx_PDB_helix_class.value(i),
            details: details.valueKind(i) === Column.ValueKinds.Present ? details.value(i) : void 0
        };
        const entry: SecondaryStructureEntry = {
            startSeqId: beg_seq_id.value(i),
            startInsCode: pdbx_beg_PDB_ins_code.value(i),
            endSeqId: end_seq_id.value(i),
            endInsCode: pdbx_end_PDB_ins_code.value(i),
            type,
            key: elements.length
        };

        elements[elements.length] = element;

        const asymId = beg_label_asym_id.value(i)!;
        if (map.has(asymId)) {
            const entries = map.get(asymId)!;
            if (entries.has(entry.startSeqId)) {
                entries.get(entry.startSeqId)!.push(entry);
            } else {
                entries.set(entry.startSeqId, [entry]);
            }
        } else {
            map.set(asymId, new Map([[entry.startSeqId, [entry]]]));
        }
    }
}

function addSheets(cat: StructSheetRange, coordinates: CoordinateType, map: SecondaryStructureMap, sheetCount: number, elements: SecondaryStructure.Element[]) {
    if (!cat._rowCount) return;

    const { beg_label_asym_id, beg_label_seq_id, beg_auth_seq_id, pdbx_beg_PDB_ins_code } = cat;
    const { end_label_seq_id, end_auth_seq_id, pdbx_end_PDB_ins_code } = cat;
    const { sheet_id } = cat;

    const beg_seq_id = coordinates === 'label' ? beg_label_seq_id : beg_auth_seq_id;
    const end_seq_id = coordinates === 'label' ? end_label_seq_id : end_auth_seq_id;

    const sheet_id_key = new Map<string, number>();
    let currentKey = sheetCount + 1;

    for (let i = 0, _i = cat._rowCount; i < _i; i++) {
        const id = sheet_id.value(i);
        let key: number;
        if (sheet_id_key.has(id)) key = sheet_id_key.get(id)!;
        else {
            key = currentKey++;
            sheet_id_key.set(id, key);
        }

        const type = SecondaryStructureType.create(SecondaryStructureType.Flag.Beta | SecondaryStructureType.Flag.BetaSheet);
        const element: SecondaryStructure.Sheet = {
            kind: 'sheet',
            flags: type,
            sheet_id: id,
            symmetry: void 0
        };
        const entry: SecondaryStructureEntry = {
            startSeqId: beg_seq_id.value(i),
            startInsCode: pdbx_beg_PDB_ins_code.value(i),
            endSeqId: end_seq_id.value(i),
            endInsCode: pdbx_end_PDB_ins_code.value(i),
            type,
            key: elements.length
        };

        elements[elements.length] = element;

        const asymId = beg_label_asym_id.value(i)!;
        if (map.has(asymId)) {
            const entries = map.get(asymId)!;
            if (entries.has(entry.startSeqId)) {
                entries.get(entry.startSeqId)!.push(entry);
            } else {
                entries.set(entry.startSeqId, [entry]);
            }
        } else {
            map.set(asymId, new Map([[entry.startSeqId, [entry]]]));
        }
    }

    return;
}

function assignSecondaryStructureEntry(hierarchy: AtomicHierarchy, coordinates: CoordinateType, entry: SecondaryStructureEntry, resStart: ResidueIndex, resEnd: ResidueIndex, data: SecondaryStructureData) {
    const { auth_seq_id, label_seq_id, pdbx_PDB_ins_code } = hierarchy.residues;
    const { endSeqId, endInsCode, key, type } = entry;

    const seq_id = coordinates === 'label' ? label_seq_id : auth_seq_id;

    let rI = resStart;
    while (rI < resEnd) {
        const seqNumber = seq_id.value(rI);
        data.type[rI] = type;
        data.key[rI] = key;

        if ((seqNumber > endSeqId) ||
            (seqNumber === endSeqId && pdbx_PDB_ins_code.value(rI) === endInsCode)) {
            break;
        }

        rI++;
    }
}

function assignSecondaryStructureRanges(hierarchy: AtomicHierarchy, coordinates: CoordinateType, map: SecondaryStructureMap, data: SecondaryStructureData) {
    const { count: chainCount } = hierarchy.chainAtomSegments;
    const { label_asym_id } = hierarchy.chains;
    const { auth_seq_id, label_seq_id, pdbx_PDB_ins_code } = hierarchy.residues;
    const seq_id = coordinates === 'label' ? label_seq_id : auth_seq_id;

    for (let cI = 0 as ChainIndex; cI < chainCount; cI++) {
        const resStart = AtomicHierarchy.chainStartResidueIndex(hierarchy, cI), resEnd = AtomicHierarchy.chainEndResidueIndexExcl(hierarchy, cI);
        const asymId = label_asym_id.value(cI);
        if (map.has(asymId)) {
            const entries = map.get(asymId)!;

            for (let rI = resStart; rI < resEnd; rI++) {
                const seqId = seq_id.value(rI);
                if (entries.has(seqId)) {
                    const entryList = entries.get(seqId)!;
                    for (const entry of entryList) {
                        const insCode = pdbx_PDB_ins_code.value(rI);
                        if (entry.startInsCode !== insCode) continue;
                        assignSecondaryStructureEntry(hierarchy, coordinates, entry, rI, resEnd, data);
                    }
                }
            }
        }
    }
}