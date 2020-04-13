/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../db';

export interface Grouping<V, K> {
    map: Map<K, V[]>,
    keys: ReadonlyArray<K>,
    groups: ReadonlyArray<ReadonlyArray<V>>
}

class GroupingImpl<K, V> {
    readonly map = new Map<K, V[]>();
    readonly keys: K[] = [];
    readonly groups: V[][] = [];

    add(a: V) {
        const key = this.getKey(a) as any;
        if (!!this.map.has(key)) {
            const group = this.map.get(key)!;
            group[group.length] = a;
        } else {
            const group = [a];
            this.map.set(key, group);
            this.keys[this.keys.length] = key;
            this.groups[this.groups.length] = group;
        }
    }

    getGrouping(): Grouping<V, K> {
        return { keys: this.keys, groups: this.groups, map: this.map };
    }

    constructor(private getKey: (v: V) => K) { }
}

export function Grouper<V, K>(getKey: (x: V) => K) {
    return new GroupingImpl<K, V>(getKey);
}

export function groupBy<V, K>(values: ArrayLike<V> | Column<V>, getKey: (x: V) => K) {
    const gs = Grouper(getKey);
    if (Column.is(values)) {
        const v = values.value;
        for (let i = 0, _i = values.rowCount; i < _i; i++) gs.add(v(i));
    } else {
        for (let i = 0, _i = values.length; i < _i; i++) gs.add(values[i]);
    }
    return gs.getGrouping();
}