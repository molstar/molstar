/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { CifExportContext } from '../../../export/mmcif';
import { QuerySymbolRuntime } from 'mol-script/runtime/query/compiler';

interface ModelPropertyDescriptor<Symbols extends { [name: string]: QuerySymbolRuntime } = { }> {
    readonly isStatic: boolean,
    readonly name: string,

    cifExport?: {
        // Prefix enforced during export.
        prefix: string,
        categories: CifWriter.Category<CifExportContext>[]
    },

    // TODO: add aliases when lisp-like mol-script is done
    symbols?: Symbols
}

function ModelPropertyDescriptor<Desc extends ModelPropertyDescriptor>(desc: Desc) {
    return desc;
}

export { ModelPropertyDescriptor }