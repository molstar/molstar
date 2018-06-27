/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { CifExportContext } from '../../../export/mmcif';

interface ModelPropertyDescriptor {
    readonly isStatic: boolean,
    readonly name: string,

    cifExport: {
        /** used category names that can be used for "filtering" by the writer */
        readonly categoryNames: ReadonlyArray<string>,
        categoryProvider: (ctx: CifExportContext) => CifWriter.Category.Provider[]
    }
}

export { ModelPropertyDescriptor }