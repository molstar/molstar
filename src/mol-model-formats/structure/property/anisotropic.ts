/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Table } from '../../../mol-data/db';
import { Model, CustomPropertyDescriptor } from '../../../mol-model/structure';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { CifWriter } from '../../../mol-io/writer/cif';
import { FormatPropertyProvider } from '../common/property';

export { AtomSiteAnisotrop }

type Anisotrop = Table<mmCIF_Schema['atom_site_anisotrop']>

interface AtomSiteAnisotrop {
    data: Anisotrop
    /** maps atom_site-index to atom_site_anisotrop-index */
    elementToAnsiotrop: Int32Array
}

namespace AtomSiteAnisotrop {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'atom_site_anisotrop',
        cifExport: {
            prefix: '',
            categories: [{
                name: 'atom_site_anisotrop',
                instance(ctx) {
                    const p = Provider.get(ctx.firstModel);
                    if (!p) return CifWriter.Category.Empty;
                    // TODO filter to write only data for elements that exist in model
                    return CifWriter.Category.ofTable(p.data);
                }
            }]
        }
    };

    export const Provider = FormatPropertyProvider.create<AtomSiteAnisotrop>(Descriptor)

    export function getElementToAnsiotrop(model: Model, data: Anisotrop) {
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
}