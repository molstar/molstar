/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { CifExportContext } from '../../../export/mmcif';
import { QuerySymbolRuntime } from 'mol-script/runtime/query/compiler';
import { UUID } from 'mol-util';

interface ModelPropertyDescriptor<ExportCtx = CifExportContext, Symbols extends { [name: string]: QuerySymbolRuntime } = { }> {
    readonly isStatic: boolean,
    readonly name: string,

    cifExport?: {
        // Prefix enforced during export.
        prefix: string,
        context?: (ctx: CifExportContext) => ExportCtx | undefined,
        categories: CifWriter.Category<ExportCtx>[]
    },

    // TODO: add aliases when lisp-like mol-script is done
    symbols?: Symbols
}

function ModelPropertyDescriptor<Ctx, Desc extends ModelPropertyDescriptor<Ctx>>(desc: Desc) {
    return desc;
}

namespace ModelPropertyDescriptor {
    export function getUUID(prop: ModelPropertyDescriptor): UUID {
        if (!(prop as any).__key) {
            (prop as any).__key = UUID.create();
        }
        return (prop as any).__key;
    }
}

export { ModelPropertyDescriptor }