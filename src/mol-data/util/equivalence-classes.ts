/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export class EquivalenceClassesImpl<K, V> {
    private id = 0;
    private byHash = new Map<number, { id: number, keys: K[], value: V }[]>();

    readonly groups: K[][] = [];

    private createGroup(key: K, value: V) {
        const id = this.id++;
        const keys = [key];
        this.groups[id] = keys;
        return { id, keys, value };
    }

    // Return the group representative.
    add(key: K, a: V) {
        const hash = this.getHash(a);
        if (this.byHash.has(hash)) {
            const groups = this.byHash.get(hash)!;
            for (let i = 0, _i = groups.length; i < _i; i++) {
                const group = groups[i];
                if (this.areEqual(a, group.value)) {
                    group.keys[group.keys.length] = key;
                    return group.value;
                }
            }
            const group = this.createGroup(key, a);
            groups[groups.length] = group;
            return group.value;
        } else {
            const group = this.createGroup(key, a);
            this.byHash.set(hash, [group]);
            return group.value;
        }
    }

    constructor(private getHash: (v: V) => any, private areEqual: (a: V, b: V) => boolean) { }
}

export function EquivalenceClasses<K, V>(getHash: (x: V) => any, areEqual: (a: V, b: V) => boolean) {
    return new EquivalenceClassesImpl<K, V>(getHash, areEqual);
}
