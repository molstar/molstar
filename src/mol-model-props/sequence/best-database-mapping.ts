/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../mol-data/db';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { Model } from '../../mol-model/structure';
import { StructureElement } from '../../mol-model/structure/structure';
import { CustomModelProperty } from '../common/custom-model-property';

export { BestDatabaseSequenceMapping };

interface BestDatabaseSequenceMapping {
    readonly dbName: string[],
    readonly accession: string[],
    readonly num: number[],
    readonly residue: string[]
}

namespace BestDatabaseSequenceMapping {
    export const Provider: CustomModelProperty.Provider<{}, BestDatabaseSequenceMapping> = CustomModelProperty.createProvider({
        label: 'Best Database Sequence Mapping',
        descriptor: CustomPropertyDescriptor({
            name: 'molstar_best_database_sequence_mapping'
        }),
        type: 'static',
        defaultParams: {},
        getParams: () => ({}),
        isApplicable: (data: Model) => MmcifFormat.is(data.sourceData) && data.sourceData.data.frame.categories?.atom_site?.fieldNames.indexOf('db_name') >= 0,
        obtain: async (ctx, data) => {
            return { value: fromCif(data) };
        }
    });

    export function getKey(loc: StructureElement.Location) {
        const model = loc.unit.model;
        const data = Provider.get(model).value;
        if (!data) return '';
        const eI = loc.unit.elements[loc.element];
        const rI = model.atomicHierarchy.residueAtomSegments.offsets[eI];
        return data.accession[rI];
    }

    export function getLabel(loc: StructureElement.Location) {
        const model = loc.unit.model;
        const data = Provider.get(model).value;
        if (!data) return;
        const eI = loc.unit.elements[loc.element];
        const rI = model.atomicHierarchy.residueAtomSegments.offsets[eI];
        const dbName = data.dbName[rI];
        if (!dbName) return;
        return `${dbName} ${data.accession[rI]} ${data.num[rI]} ${data.residue[rI]}`;
    }

    function fromCif(model: Model): BestDatabaseSequenceMapping | undefined {
        if (!MmcifFormat.is(model.sourceData)) return;

        const { atom_site } = model.sourceData.data.frame.categories;
        const db_name = atom_site.getField('db_name');
        const db_acc = atom_site.getField('db_acc');
        const db_num = atom_site.getField('db_num');
        const db_res = atom_site.getField('db_res');

        if (!db_name || !db_acc || !db_num || !db_res) return;

        const { atomSourceIndex } = model.atomicHierarchy;
        const { count, offsets: residueOffsets } = model.atomicHierarchy.residueAtomSegments;
        const dbName = new Array<string>(count);
        const accession = new Array<string>(count);
        const num = new Array<number>(count);
        const residue = new Array<string>(count);

        for (let i = 0; i < count; i++) {
            const row = atomSourceIndex.value(residueOffsets[i]);

            if (db_name.valueKind(row) !== Column.ValueKind.Present) {
                dbName[row] = '';
                accession[row] = '';
                num[row] = 0;
                residue[row] = '';
                continue;
            }

            dbName[row] = db_name.str(row);
            accession[row] = db_acc.str(row);
            num[row] = db_num.int(row);
            residue[row] = db_res.str(row);
        }

        return { dbName, accession, num, residue };
    }
}
