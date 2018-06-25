/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PropertyDescriptor } from './descriptor'

export class CustomProperties {
    private _list: PropertyDescriptor[] = [];
    private _set = new Set<PropertyDescriptor>();

    get all(): ReadonlyArray<PropertyDescriptor> {
        return this._list;
    }

    add(desc: PropertyDescriptor) {
        this._list.push(desc);
        this._set.add(desc);
    }

    has(desc: PropertyDescriptor): boolean {
        return this._set.has(desc);
    }
}