/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Table, Column } from '../../../mol-data/db';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { CifWriter } from '../../../mol-io/writer/cif';
import { FormatPropertyProvider } from '../common/property';
import { MmcifFormat } from '../mmcif';

export { AtomSiteAnisotrop };

const Anisotrop = {
    U: mmCIF_Schema.atom_site_anisotrop.U,
    U_esd: mmCIF_Schema.atom_site_anisotrop.U_esd
};
type Anisotrop = Table<typeof Anisotrop>

interface AtomSiteAnisotrop {
    data: Anisotrop
    /** maps atom_site-index to atom_site_anisotrop-index */
    elementToAnsiotrop: Int32Array
}

namespace AtomSiteAnisotrop {
    export const Schema = Anisotrop;

    export const Descriptor: CustomPropertyDescriptor = {
        name: 'atom_site_anisotrop',
        cifExport: {
            prefix: '',
            categories: [{
                name: 'atom_site_anisotrop',
                instance(ctx) {
                    const p = Provider.get(ctx.firstModel);
                    if (!p) return CifWriter.Category.Empty;
                    if (!MmcifFormat.is(ctx.firstModel.sourceData)) return CifWriter.Category.Empty;
                    // TODO filter to write only data for elements that exist in model
                    return CifWriter.Category.ofTable(ctx.firstModel.sourceData.data.db.atom_site_anisotrop);
                }
            }]
        }
    };

    export const Provider = FormatPropertyProvider.create<AtomSiteAnisotrop>(Descriptor);

    export function getElementToAnsiotrop(atomId: Column<number>, ansioId: Column<number>) {
        const atomIdToElement = new Int32Array(atomId.rowCount);
        atomIdToElement.fill(-1);
        for (let i = 0, il = atomId.rowCount; i < il; i++) {
            atomIdToElement[atomId.value(i)] = i;
        }

        const elementToAnsiotrop = new Int32Array(atomId.rowCount);
        elementToAnsiotrop.fill(-1);
        for (let i = 0, il = ansioId.rowCount; i < il; ++i) {
            const ei = atomIdToElement[ansioId.value(i)];
            if (ei !== -1) elementToAnsiotrop[ei] = i;
        }

        return elementToAnsiotrop;
    }

    export function getElementToAnsiotropFromLabel(atomLabel: Column<string>, ansioLabel: Column<string>) {
        const atomLabelToElement: { [k: string]: number | undefined } = {};
        for (let i = 0, il = atomLabel.rowCount; i < il; i++) {
            atomLabelToElement[atomLabel.value(i)] = i;
        }

        const elementToAnsiotrop = new Int32Array(atomLabel.rowCount);
        elementToAnsiotrop.fill(-1);
        for (let i = 0, il = ansioLabel.rowCount; i < il; ++i) {
            const ei = atomLabelToElement[ansioLabel.value(i)];
            if (ei !== undefined) elementToAnsiotrop[ei] = i;
        }

        return elementToAnsiotrop;
    }
}