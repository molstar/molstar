/**
 * Copyright (c) 2021-23 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../mol-data/db';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { Model } from '../../mol-model/structure';
import { StructureElement } from '../../mol-model/structure/structure';
import { CustomModelProperty } from '../common/custom-model-property';

export { SIFTSMapping as SIFTSMapping };

interface SIFTSMappingMapping {
    readonly dbName: string[],
    readonly accession: string[],
    readonly num: string[],
    readonly residue: string[]
}

namespace SIFTSMapping {
    export const Provider: CustomModelProperty.Provider<{}, SIFTSMappingMapping> = CustomModelProperty.createProvider({
        label: 'SIFTS Mapping',
        descriptor: CustomPropertyDescriptor({
            name: 'sifts_sequence_mapping'
        }),
        type: 'static',
        defaultParams: {},
        getParams: () => ({}),
        isApplicable: (data: Model) => isAvailable(data),
        obtain: async (ctx, data) => {
            return { value: fromCif(data) };
        }
    });

    export function isAvailable(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) return false;

        const {
            pdbx_sifts_xref_db_name: db_name,
            pdbx_sifts_xref_db_acc: db_acc,
            pdbx_sifts_xref_db_num: db_num,
            pdbx_sifts_xref_db_res: db_res
        } = model.sourceData.data.db.atom_site;

        return db_name.isDefined && db_acc.isDefined && db_num.isDefined && db_res.isDefined;
    }

    export function getKey(loc: StructureElement.Location) {
        const model = loc.unit.model;
        const data = Provider.get(model).value;
        if (!data) return '';
        const rI = model.atomicHierarchy.residueAtomSegments.index[loc.element];
        return data.accession[rI];
    }

    export function getLabel(loc: StructureElement.Location) {
        const model = loc.unit.model;
        const data = Provider.get(model).value;
        if (!data) return;
        const rI = model.atomicHierarchy.residueAtomSegments.index[loc.element];
        const dbName = data.dbName[rI];

        if (!dbName) return;
        return `${dbName} ${data.accession[rI]} ${data.num[rI]} ${data.residue[rI]}`;
    }

    function fromCif(model: Model): SIFTSMappingMapping | undefined {
        if (!MmcifFormat.is(model.sourceData)) return;

        const {
            pdbx_sifts_xref_db_name: db_name,
            pdbx_sifts_xref_db_acc: db_acc,
            pdbx_sifts_xref_db_num: db_num,
            pdbx_sifts_xref_db_res: db_res
        } = model.sourceData.data.db.atom_site;

        if (!db_name.isDefined || !db_acc.isDefined || !db_num.isDefined || !db_res.isDefined) return;

        const { atomSourceIndex } = model.atomicHierarchy;
        const { count, offsets: residueOffsets } = model.atomicHierarchy.residueAtomSegments;
        const dbName = new Array<string>(count);
        const accession = new Array<string>(count);
        const num = new Array<string>(count);
        const residue = new Array<string>(count);

        for (let i = 0; i < count; i++) {
            const row = atomSourceIndex.value(residueOffsets[i]);

            if (db_name.valueKind(row) !== Column.ValueKinds.Present) {
                dbName[i] = '';
                accession[i] = '';
                num[i] = '';
                residue[i] = '';
                continue;
            }

            dbName[i] = db_name.value(row);
            accession[i] = db_acc.value(row);
            num[i] = db_num.value(row);
            residue[i] = db_res.value(row);
        }

        return { dbName, accession, num, residue };
    }
}
