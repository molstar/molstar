/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LinkedList } from '../mol-data/generic';
import { GraphicsRenderObject } from './render-object';
import { Renderable } from './renderable';

type RenderObjectItem = LinkedList.Node<GraphicsRenderObject>
type ReadyItem = { object: GraphicsRenderObject, renderable: Renderable<any> }

export class CommitQueue {
    private removeList = LinkedList<GraphicsRenderObject>();
    private removeMap = new Map<GraphicsRenderObject, RenderObjectItem>();
    private addList = LinkedList<GraphicsRenderObject>();
    private addMap = new Map<GraphicsRenderObject, RenderObjectItem>();
    private readyList = LinkedList<ReadyItem>();
    private readyMap = new Map<Renderable<any>, LinkedList.Node<ReadyItem>>();

    get isEmpty() {
        return this.removeList.count === 0 && this.addList.count === 0 && this.readyList.count === 0;
    }

    get size() {
        return this.removeMap.size + this.addMap.size;
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

    ready(e: ReadyItem) {
        if (this.readyMap.has(e.renderable)) return;
        const b = this.readyList.addLast(e);
        this.readyMap.set(e.renderable, b);
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

    tryGetReady() {
        const o = this.readyList.removeFirst();
        if (o) this.readyMap.delete(o.renderable);
        return o;
    }
}
