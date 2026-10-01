/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface LinkedList<T> {
    readonly count: number,
    readonly first: LinkedList.Node<T> | null,
    readonly last: LinkedList.Node<T> | null,
    addFirst(value: T): LinkedList.Node<T>,
    addLast(value: T): LinkedList.Node<T>,
    remove(node: LinkedList.Node<T>): void,
    removeFirst(): T | undefined,
    removeLast(): T | undefined,
    find(value: T): LinkedList.Node<T> | undefined
}

function LinkedList<T>(): LinkedList<T> {
    return new LinkedListImpl();
}

namespace LinkedList {
    export interface Node<T> {
        previous: Node<T> | null,
        next: Node<T> | null,
        inList: boolean,
        value: T
    }
}

function createListNode<T>(value: T): LinkedList.Node<T> {
    return { previous: null, next: null, inList: true, value };
}

class LinkedListImpl<T> implements LinkedList<T> {
    count: number = 0;
    first: LinkedList.Node<T> | null = null;
    last: LinkedList.Node<T> | null = null;

    addFirst(value: T) {
        const node = createListNode(value);
        node.inList = true;
        if (this.first) this.first.previous = node;
        node.next = this.first;
        this.first = node;
        this.count++;
        if (!this.last) this.last = node;
        return node;
    }

    addLast(value: T) {
        const node = createListNode(value);
        if (this.last !== null) {
            this.last.next = node;
        }
        node.previous = this.last;
        this.last = node;
        if (this.first === null) {
            this.first = node;
        }
        node.inList = true;
        this.count++;
        return node;
    }

    removeFirst(): T | undefined {
        const fst = this.first;
        if (fst) {
            this.remove(fst);
            return fst.value;
        }
        return void 0;
    }

    removeLast(): T | undefined {
        const last = this.last;
        if (last) {
            this.remove(last);
            return last.value;
        }
        return void 0;
    }

    remove(node: LinkedList.Node<T>) {
        if (!node.inList) return;

        node.inList = false;

        if (node.previous !== null) {
            node.previous.next = node.next;
        } else if (/* first == item*/ node.previous === null) {
            this.first = node.next;
        }

        if (node.next !== null) {
            node.next.previous = node.previous;
        } else if (/* last == item*/ node.next === null) {
            this.last = node.previous;
        }

        node.next = null;
        node.previous = null;
        this.count--;
    }

    find(value: T): LinkedList.Node<T> | undefined {
        let current = this.first;
        while (current !== null) {
            if (current.value === value) return current;
            current = current.next;
        }
        return void 0;
    }
}

export { LinkedList };