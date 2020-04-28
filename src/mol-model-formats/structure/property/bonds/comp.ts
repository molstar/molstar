/**
 * Copyright (c) 2017-2020 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../../mol-model/structure/model/model';
import { BondType } from '../../../../mol-model/structure/model/types';
import { CustomPropertyDescriptor } from '../../../../mol-model/custom-property';
import { mmCIF_Schema } from '../../../../mol-io/reader/cif/schema/mmcif';
import { CifWriter } from '../../../../mol-io/writer/cif';
import { Table } from '../../../../mol-data/db';
import { FormatPropertyProvider } from '../../common/property';

export interface ComponentBond {
    readonly data: Table<mmCIF_Schema['chem_comp_bond']>
    readonly entries: ReadonlyMap<string, ComponentBond.Entry>
}

export namespace ComponentBond {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'chem_comp_bond',
        cifExport: {
            prefix: '',
            categories: [{
                name: 'chem_comp_bond',
                instance(ctx) {
                    const p = Provider.get(ctx.firstModel);
                    if (!p) return CifWriter.Category.Empty;
                    const chem_comp_bond = p.data;
                    if (!chem_comp_bond) return CifWriter.Category.Empty;

                    const comp_names = ctx.structures[0].uniqueResidueNames;
                    const { comp_id, _rowCount } = chem_comp_bond;
                    const indices: number[] = [];
                    for (let i = 0; i < _rowCount; i++) {
                        if (comp_names.has(comp_id.value(i))) indices[indices.length] = i;
                    }

                    return CifWriter.Category.ofTable(chem_comp_bond, indices);
                }
            }]
        }
    };

    export const Provider = FormatPropertyProvider.create<ComponentBond>(Descriptor);

    export function chemCompBondFromTable(model: Model, table: Table<mmCIF_Schema['chem_comp_bond']>): Table<mmCIF_Schema['chem_comp_bond']> {
        return Table.pick(table, mmCIF_Schema.chem_comp_bond, (i: number) => {
            return model.properties.chemicalComponentMap.has(table.comp_id.value(i));
        });
    }

    export function getEntriesFromChemCompBond(data: Table<mmCIF_Schema['chem_comp_bond']>) {
        const entries: Map<string, Entry> = new Map();

        function addEntry(id: string) {
            let e = new Entry(id);
            entries.set(id, e);
            return e;
        }

        const { comp_id, atom_id_1, atom_id_2, value_order, pdbx_aromatic_flag, _rowCount } = data;

        let entry = addEntry(comp_id.value(0)!);
        for (let i = 0; i < _rowCount; i++) {
            const id = comp_id.value(i)!;
            const nameA = atom_id_1.value(i)!;
            const nameB = atom_id_2.value(i)!;
            const order = value_order.value(i)!;
            const aromatic = pdbx_aromatic_flag.value(i) === 'Y';

            if (entry.id !== id) {
                entry = addEntry(id);
            }

            let flags: number = BondType.Flag.Covalent;
            let ord = 1;
            if (aromatic) flags |= BondType.Flag.Aromatic;
            switch (order.toLowerCase()) {
                case 'doub':
                case 'delo':
                    ord = 2;
                    break;
                case 'trip': ord = 3; break;
                case 'quad': ord = 4; break;
            }

            entry.add(nameA, nameB, ord, flags);
        }

        return entries;
    }

    export class Entry {
        readonly map: Map<string, Map<string, { order: number, flags: number }>> = new Map();

        add(a: string, b: string, order: number, flags: number, swap = true) {
            let e = this.map.get(a);
            if (e !== void 0) {
                let f = e.get(b);
                if (f === void 0) {
                    e.set(b, { order, flags });
                }
            } else {
                let map = new Map<string, { order: number, flags: number }>();
                map.set(b, { order, flags });
                this.map.set(a, map);
            }

            if (swap) this.add(b, a, order, flags, false);
        }

        constructor(public readonly id: string) { }
    }
}