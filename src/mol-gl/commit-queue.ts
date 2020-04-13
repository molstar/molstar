/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LinkedList } from '../mol-data/generic';
import { GraphicsRenderObject } from './render-object';

type N = LinkedList.Node<GraphicsRenderObject>

export class CommitQueue {
    private removeList = LinkedList<GraphicsRenderObject>();
    private removeMap = new Map<GraphicsRenderObject, N>();
    private addList = LinkedList<GraphicsRenderObject>();
    private addMap = new Map<GraphicsRenderObject, N>();

    get isEmpty() {
        return this.removeList.count === 0 && this.addList.count === 0;
    }

    add(o: GraphicsRenderObject) {
        if (this.removeMap.has(o)) {
            const a = this.removeMap.get(o)!;
            this.removeMap.delete(o);
            this.removeList.remove(a);
        }
        if (this.addMap.has(o)) return;
        const b = this.addList.addLast(o);
        this.addMap.set(o, b);
    }

    remove(o: GraphicsRenderObject) {
        if (this.addMap.has(o)) {
            const a = this.addMap.get(o)!;
            this.addMap.delete(o);
            this.addList.remove(a);
        }
        if (this.removeMap.has(o)) return;
        const b = this.removeList.addLast(o);
        this.removeMap.set(o, b);
    }

    tryGetRemove() {
        const o = this.removeList.removeFirst();
        if (o) this.removeMap.delete(o);
        return o;
    }

    tryGetAdd() {
        const o = this.addList.removeFirst();
        if (o) this.addMap.delete(o);
        return o;
    }
}
