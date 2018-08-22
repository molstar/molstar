/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ModelPropertyDescriptor } from './descriptor'

export class CustomProperties {
    private _list: ModelPropertyDescriptor[] = [];
    private _set = new Set<ModelPropertyDescriptor>();

    get all(): ReadonlyArray<ModelPropertyDescriptor> {
        return this._list;
    }

    add(desc: ModelPropertyDescriptor) {
        if (this._set.has(desc)) return;

        this._list.push(desc);
        this._set.add(desc);
    }

    has(desc: ModelPropertyDescriptor): boolean {
        return this._set.has(desc);
    }
}