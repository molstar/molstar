/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from '../mol-io/writer/cif';
import { CifExportContext } from './structure/export/mmcif';
import { QuerySymbolRuntime } from '../mol-script/runtime/query/compiler';
import { UUID } from '../mol-util';
import { Asset } from '../mol-util/assets';

export { CustomPropertyDescriptor, CustomProperties };

interface CustomPropertyDescriptor<ExportCtx = CifExportContext, Symbols extends { [name: string]: QuerySymbolRuntime } = { }> {
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
    private _refs = new Map<CustomPropertyDescriptor, number>();
    private _assets = new Map<CustomPropertyDescriptor, Asset.Wrapper[]>();

    get all(): ReadonlyArray<CustomPropertyDescriptor> {
        return this._list;
    }

    add(desc: CustomPropertyDescriptor<any>) {
        if (this._set.has(desc)) return;

        this._list.push(desc);
        this._set.add(desc);
    }

    reference(desc: CustomPropertyDescriptor<any>, add: boolean) {
        let refs = this._refs.get(desc) || 0;
        refs += add ? 1 : -1;
        this._refs.set(desc, Math.max(refs, 0));
    }

    hasReference(desc: CustomPropertyDescriptor<any>) {
        return (this._refs.get(desc) || 0) > 0;
    }

    has(desc: CustomPropertyDescriptor<any>): boolean {
        return this._set.has(desc);
    }

    /** Sets assets for a prop, disposes of existing assets for that prop */
    assets(desc: CustomPropertyDescriptor<any>, assets?: Asset.Wrapper[]) {
        const prevAssets = this._assets.get(desc);
        if (prevAssets) {
            for (const a of prevAssets) a.dispose();
        }
        if (assets) this._assets.set(desc, assets);
        else this._assets.delete(desc);
    }

    /** Disposes of all assets of all props */
    dispose() {
        this._assets.forEach(assets => {
            for (const a of assets) a.dispose();
        });
    }
}