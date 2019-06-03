/**
 * Copyright (c) 2017-2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../../mol-model/structure/model/model'
import { LinkType } from '../../../../mol-model/structure/model/types'
import { CustomPropertyDescriptor } from '../../../../mol-model/structure';
import { mmCIF_Database } from '../../../../mol-io/reader/cif/schema/mmcif';
import { Structure, Unit, StructureProperties, StructureElement } from '../../../../mol-model/structure';
import { Segmentation } from '../../../../mol-data/int';
import { CifWriter } from '../../../../mol-io/writer/cif'

export interface ComponentBond {
    entries: Map<string, ComponentBond.Entry>
}

export namespace ComponentBond {
    export const Descriptor: CustomPropertyDescriptor = {
        isStatic: true,
        name: 'chem_comp_bond',
        cifExport: {
            prefix: '',
            categories: [{
                name: 'chem_comp_bond',
                instance(ctx) {
                    const chem_comp_bond = getChemCompBond(ctx.structures[0].model);
                    if (!chem_comp_bond) return CifWriter.Category.Empty;

                    const comp_names = getUniqueResidueNames(ctx.structures[0]);
                    const { comp_id, _rowCount } = chem_comp_bond;
                    const indices: number[] = [];
                    for (let i = 0; i < _rowCount; i++) {
                        if (comp_names.has(comp_id.value(i))) indices[indices.length] = i;
                    }

                    return CifWriter.Category.ofTable(chem_comp_bond, indices)
                }
            }]
        }
    }

    export function attachFromMmCif(model: Model): boolean {
        if (model.customProperties.has(Descriptor)) return true;
        if (model.sourceData.kind !== 'mmCIF') return false;
        const { chem_comp_bond } = model.sourceData.data;
        if (chem_comp_bond._rowCount === 0) return false;

        model.customProperties.add(Descriptor);
        model._staticPropertyData.__ComponentBondData__ = chem_comp_bond;
        return true;
    }

    export function attachFromExternalData(model: Model, bonds: ComponentBond, force = false) {
        if (!force && model.customProperties.has(Descriptor)) return true;
        if (model._staticPropertyData.__ComponentBondData__) delete model._staticPropertyData.__ComponentBondData__;
        model.customProperties.add(Descriptor);
        model._staticPropertyData[PropName] = bonds;
        return true;
    }

    export class ComponentBondImpl implements ComponentBond {
        entries: Map<string, ComponentBond.Entry> = new Map();

        addEntry(id: string) {
            let e = new Entry(id);
            this.entries.set(id, e);
            return e;
        }
    }

    export class Entry {
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

    export function parseChemCompBond(data: mmCIF_Database['chem_comp_bond']): ComponentBond {
        const { comp_id, atom_id_1, atom_id_2, value_order, pdbx_aromatic_flag, _rowCount: rowCount } = data;

        const compBond = new ComponentBondImpl();
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

        return compBond;
    }

    function getChemCompBond(model: Model) {
        return model._staticPropertyData.__ComponentBondData__ as mmCIF_Database['chem_comp_bond'];
    }

    export const PropName = '__ComponentBond__';
    export function get(model: Model): ComponentBond | undefined {
        if (model._staticPropertyData[PropName]) return model._staticPropertyData[PropName];
        if (!model.customProperties.has(Descriptor)) return void 0;

        const chem_comp_bond = getChemCompBond(model);
        if (!chem_comp_bond) return void 0;

        const chemComp = parseChemCompBond(chem_comp_bond);
        model._staticPropertyData[PropName] = chemComp;
        return chemComp;
    }

    function getUniqueResidueNames(s: Structure) {
        const prop = StructureProperties.residue.label_comp_id;
        const names = new Set<string>();
        const loc = StructureElement.create();
        for (const unit of s.units) {
            if (!Unit.isAtomic(unit)) continue;
            const residues = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
            loc.unit = unit;
            while (residues.hasNext) {
                const seg = residues.move();
                loc.element = unit.elements[seg.start];
                names.add(prop(loc));
            }
        }
        return names;
    }
}