/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { FibonacciHeap } from '../fibonacci-heap';

describe('fibonacci-heap', () => {
    it('basic', () => {
        const heap = new FibonacciHeap();
        heap.insert(1, 2);
        heap.insert(4);
        heap.insert(2);
        heap.insert(3);
        expect(heap.size()).toBe(4);
        const node = heap.extractMinimum();
        expect(node!.key).toBe(1);
        expect(node!.value).toBe(2);
        expect(heap.size()).toBe(3);
    });
});
