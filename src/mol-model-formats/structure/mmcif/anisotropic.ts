/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Table } from '../../../mol-data/db';
import { Model, CustomPropertyDescriptor } from '../../../mol-model/structure';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { CifWriter } from '../../../mol-io/writer/cif';

export interface AtomSiteAnisotrop {
    data: Table<AtomSiteAnisotrop.Schema['atom_site_anisotrop']>
    /** maps atom_site-index to atom_site_anisotrop-index */
    elementToAnsiotrop: Int32Array
}

export namespace AtomSiteAnisotrop {
    export function getAtomSiteAnisotrop(model: Model) {
        if (model.sourceData.kind !== 'mmCIF') return void 0;
        const { atom_site_anisotrop } = model.sourceData.data
        return Table.ofColumns(Schema.atom_site_anisotrop, atom_site_anisotrop);
    }

    export const PropName = '__AtomSiteAnisotrop__';
    export function get(model: Model): AtomSiteAnisotrop | undefined {
        if (model._staticPropertyData[PropName]) return model._staticPropertyData[PropName]
        if (!model.customProperties.has(Descriptor)) return void 0;

        const data = getAtomSiteAnisotrop(model);
        if (!data) return void 0;

        const prop = { data, elementToAnsiotrop: getElementToAnsiotrop(model, data) }
        set(model, prop)

        return prop;
    }
    function set(model: Model, prop: AtomSiteAnisotrop) {
        (model._staticPropertyData[PropName] as AtomSiteAnisotrop) = prop;
    }

    export const Schema = { atom_site_anisotrop: mmCIF_Schema['atom_site_anisotrop'] };
    export type Schema = typeof Schema

    export const Descriptor: CustomPropertyDescriptor = {
        name: 'atom_site_anisotrop',
        cifExport: {
            prefix: '',
            categories: [{
                name: 'atom_site_anisotrop',
                instance(ctx) {
                    const atom_site_anisotrop = getAtomSiteAnisotrop(ctx.firstModel);
                    if (!atom_site_anisotrop) return CifWriter.Category.Empty;
                    return CifWriter.Category.ofTable(atom_site_anisotrop);
                }
            }]
        }
    };

    function getElementToAnsiotrop(model: Model, data: Table<Schema['atom_site_anisotrop']>) {
        const { atomId } = model.atomicConformation
        const atomIdToElement = new Int32Array(atomId.rowCount)
        atomIdToElement.fill(-1)
        for (let i = 0, il = atomId.rowCount; i < il; i++) {
            atomIdToElement[atomId.value(i)] = i
        }

        const { id } = data
        const elementToAnsiotrop = new Int32Array(atomId.rowCount)
        elementToAnsiotrop.fill(-1)
        for (let i = 0, il = id.rowCount; i < il; ++i) {
            const ei = atomIdToElement[id.value(i)]
            if (ei !== -1) elementToAnsiotrop[ei] = i
        }

        return elementToAnsiotrop
    }

    export async function attachFromMmCif(model: Model) {
        if (model.customProperties.has(Descriptor)) return true;
        if (model.sourceData.kind !== 'mmCIF') return false;
        const atomSiteAnisotrop = getAtomSiteAnisotrop(model);
        if (!atomSiteAnisotrop || atomSiteAnisotrop._rowCount === 0) return false;

        model.customProperties.add(Descriptor);
        return true;
    }
}