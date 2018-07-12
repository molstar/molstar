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
        categories: CifWriter.Category<CifExportContext>[]
    }
}

export { ModelPropertyDescriptor }