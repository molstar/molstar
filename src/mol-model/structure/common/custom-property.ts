/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { CifExportContext } from '../export/mmcif';
import { QuerySymbolRuntime } from 'mol-script/runtime/query/compiler';
import { UUID } from 'mol-util';

export { CustomPropertyDescriptor, CustomProperties }

interface CustomPropertyDescriptor<ExportCtx = CifExportContext, Symbols extends { [name: string]: QuerySymbolRuntime } = { }> {
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

function CustomPropertyDescriptor<Ctx, Desc extends CustomPropertyDescriptor<Ctx>>(desc: Desc) {
    return desc;
}

namespace CustomPropertyDescriptor {
    export function getUUID(prop: CustomPropertyDescriptor): UUID {
        if (!(prop as any).__key) {
            (prop as any).__key = UUID.create22();
        }
        return (prop as any).__key;
    }
}

class CustomProperties {
    private _list: CustomPropertyDescriptor[] = [];
    private _set = new Set<CustomPropertyDescriptor>();

    get all(): ReadonlyArray<CustomPropertyDescriptor> {
        return this._list;
    }

    add(desc: CustomPropertyDescriptor<any>) {
        if (this._set.has(desc)) return;

        this._list.push(desc);
        this._set.add(desc);
    }

    has(desc: CustomPropertyDescriptor<any>): boolean {
        return this._set.has(desc);
    }
}