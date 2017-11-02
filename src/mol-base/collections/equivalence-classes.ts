/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

class EquivalenceClassesImpl<K, V> {
    private id = 0;
    private byHash: { [hash: number]: { id: number, keys: K[], value: V }[] } = Object.create(null);

    readonly groups: K[][] = [];

    private createGroup(key: K, value: V) {
        const id = this.id++;
        const keys = [key];
        this.groups[id] = keys;
        return { id, keys, value };
    }

    add(key: K, a: V) {
        const hash = this.getHash(a);
        if (!!this.byHash[hash]) {
            const groups = this.byHash[hash];
            for (let i = 0, _i = groups.length; i < _i; i++) {
                const group = groups[i];
                if (this.areEqual(a, group.value)) {
                    group.keys[group.keys.length] = key;
                    return group.id;
                }
            }
            const group = this.createGroup(key, a);
            groups[groups.length] = group;
            return group.id;
        } else {
            const group = this.createGroup(key, a);
            this.byHash[hash] = [group];
            return group.id;
        }
    }

    constructor(private getHash: (v: V) => any, private areEqual: (a: V, b: V) => boolean) { }
}

function EquivalenceClasses<K, V>(getHash: (x: V) => any, areEqual: (a: V, b: V) => boolean) {
    return new EquivalenceClassesImpl<K, V>(getHash, areEqual);
}

export default EquivalenceClasses;