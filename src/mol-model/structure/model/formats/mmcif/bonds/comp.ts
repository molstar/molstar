/**
 * Copyright (c) 2017-2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Model from '../../../model'
import { LinkType } from '../../../types'

export interface ComponentBond {
    entries: Map<string, ComponentBond.Entry>
}

export namespace ComponentBond {
    export class ComponentBondImpl implements ComponentBond {
        entries: Map<string, ComponentBond.Entry> = new Map();

        addEntry(id: string) {
            let e = new Entry(id);
            this.entries.set(id, e);
            return e;
        }
    }

    export class Entry implements Entry {
        map: Map<string, Map<string, { order: number, flags: number }>> = new Map();

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

        constructor(public id: string) {
        }
    }

    export const PropName = '__ComponentBond__';
    export function fromModel(model: Model): ComponentBond | undefined {
        if (model._staticPropertyData[PropName]) return model._staticPropertyData[PropName];

        if (model.sourceData.kind !== 'mmCIF') return
        const { chem_comp_bond } = model.sourceData.data;
        if (!chem_comp_bond._rowCount) return void 0;

        let compBond = new ComponentBondImpl();

        const { comp_id, atom_id_1, atom_id_2, value_order, pdbx_aromatic_flag, _rowCount: rowCount } = chem_comp_bond;

        let entry = compBond.addEntry(comp_id.value(0)!);

        for (let i = 0; i < rowCount; i++) {

            const id = comp_id.value(i)!;
            const nameA = atom_id_1.value(i)!;
            const nameB = atom_id_2.value(i)!;
            const order = value_order.value(i)!;
            const aromatic = pdbx_aromatic_flag.value(i) === 'Y';

            if (entry.id !== id) {
                entry = compBond.addEntry(id);
            }

            let flags: number = LinkType.Flag.Covalent;
            let ord = 1;
            if (aromatic) flags |= LinkType.Flag.Aromatic;
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

        model._staticPropertyData[PropName] = compBond;
        return compBond;
    }
}