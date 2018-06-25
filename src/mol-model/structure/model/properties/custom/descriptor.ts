/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { CifExportContext } from '../../../export/mmcif';

interface PropertyDescriptor {
    readonly isStatic: boolean,
    readonly name: string,

    /** Given a structure, returns a list of category providers used for export. */
    getCifCategories: (ctx: CifExportContext) => CifWriter.Category.Provider[]
}

export { PropertyDescriptor }