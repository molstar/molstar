/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface Grouping<V, K> {
    keys: ReadonlyArray<K>,
    groups: ReadonlyArray<ReadonlyArray<V>>
}

class GroupingImpl<K, V> {
    private byKey = new Map<K, V[]>();
    readonly keys: K[] = [];
    readonly groups: V[][] = [];

    add(a: V) {
        const key = this.getKey(a) as any;
        if (!!this.byKey.has(key)) {
            const group = this.byKey.get(key)!;
            group[group.length] = a;
        } else {
            const group = [a];
            this.byKey.set(key, group);
            this.keys[this.keys.length] = key;
            this.groups[this.groups.length] = group;
        }
    }

    getGrouping(): Grouping<V, K> {
        return { keys: this.keys, groups: this.groups };
    }

    constructor(private getKey: (v: V) => K) { }
}

export function Grouper<V, K>(getKey: (x: V) => K) {
    return new GroupingImpl<K, V>(getKey);
}

function groupBy<V, K>(values: ArrayLike<V>, getKey: (x: V) => K) {
    const gs = Grouper(getKey);
    for (let i = 0, _i = values.length; i < _i; i++) gs.add(values[i]);
    return gs.getGrouping();
}

export default groupBy;