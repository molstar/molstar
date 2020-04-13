/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LinkedList } from '../linked-list';

describe('linked list', () => {

    function toArray<T>(list: LinkedList<T>) {
        const ret: T[] = [];
        for (let t = list.first; !!t; t = t.next) {
            ret[ret.length] = t.value;
        }
        return ret;
    }

    function create<T>(xs: T[]) {
        const list = LinkedList<T>();
        for (const x of xs) list.addLast(x);
        return list;
    }

    it('add', () => {
        const list = LinkedList<number>();
        list.addFirst(1);
        list.addLast(2);
        list.addFirst(3);
        list.addFirst(4);
        list.addLast(5);
        expect(toArray(list)).toEqual([4, 3, 1, 2, 5]);
        expect(list.count).toBe(5);
    });

    it ('remove', () => {
        const list = create([1, 2, 3, 4]);
        let fst = list.removeFirst();
        expect(fst).toBe(1);
        expect(list.last!.value).toBe(4);
        expect(list.count).toBe(3);
        expect(toArray(list)).toEqual([2, 3, 4]);

        let last = list.removeLast();
        expect(last).toBe(4);
        expect(list.last!.value).toBe(3);
        expect(list.count).toBe(2);
        expect(toArray(list)).toEqual([2, 3]);

        let n3 = list.find(3)!;
        list.remove(n3);
        expect(list.first!.value).toBe(2);
        expect(list.last!.value).toBe(2);
        expect(list.count).toBe(1);
        expect(toArray(list)).toEqual([2]);

        list.removeFirst();
        expect(list.first).toBe(null);
        expect(list.last).toBe(null);
        expect(list.count).toBe(0);
        expect(toArray(list)).toEqual([]);
    });
});