/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arraySetAdd } from '../../mol-util/array';

export { UndirectedGraph };

class UndirectedGraph<T extends string | number, S = T> {
    vertices = new Map<T, S | undefined>();
    edges = new Map<T, T[]>();

    addVertex(a: T, s?: S) {
        this.vertices.set(a, s);
        if (!this.edges.has(a)) {
            this.edges.set(a, []);
        }
    }

    addEdge(a: T, b: T) {
        arraySetAdd(this.edges.get(a)!, b);
        arraySetAdd(this.edges.get(b)!, a);
    }

    /**
     * Returns a connected component induced by a and
     * its unique representative.
     */
    getComponent(a: T) {
        let pivot = a;
        const component: T[] = [];

        const visited = new Set<T>();
        const stack = [a];
        while (stack.length > 0) {
            const e = stack.pop()!;
            component.push(e);
            visited.add(e);
            if (e < pivot) pivot = e;

            for (const b of this.edges.get(e)!) {
                if (visited.has(b)) continue;
                stack.push(b);
            }
        }

        return { pivot, component };
    }
}