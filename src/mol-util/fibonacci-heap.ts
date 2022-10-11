/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 *
 * Adapted from https://github.com/gwtw/ts-fibonacci-heap, Copyright (c) 2014 Daniel Imms, MIT
 */

interface INode<K, V> {
    key: K;
    value?: V;
}

type CompareFunction<K, V> = (a: INode<K, V>, b: INode<K, V>) => number;

class Node<K, V> implements INode<K, V> {
    public key: K;
    public value: V | undefined;
    public prev: Node<K, V>;
    public next: Node<K, V>;
    public parent: Node<K, V> | null = null;
    public child: Node<K, V> | null = null;

    public degree: number = 0;
    public isMarked: boolean = false;

    constructor(key: K, value?: V) {
        this.key = key;
        this.value = value;
        this.prev = this;
        this.next = this;
    }
}

class NodeListIterator<K, V> {
    private _index: number;
    private _items: Node<K, V>[];
    private _len: number;
    /**
   * Creates an Iterator used to simplify the consolidate() method. It works by
   * making a shallow copy of the nodes in the root list and iterating over the
   * shallow copy instead of the source as the source will be modified.
   * @param start A node from the root list.
   */
    constructor(start?: Node<K, V>) {
        this._index = -1;
        this._items = [];
        this._len = 0;
        if (start) {
            let current = start, l = 0;
            do {
                this._items[l++] = current;
                current = current.next;
            } while (start !== current);
            this._len = l;
        }
    }

    /**
   * @return Whether there is a next node in the iterator.
   */
    public hasNext(): boolean {
        return this._index < this._len - 1;
    }

    /**
   * @return The next node.
   */
    public next(): Node<K, V> {
        return this._items[++this._index];
    }

    /**
   * @return Resets iterator to reuse it.
   */
    public reset(start: Node<K, V>) {
        this._index = -1;
        this._len = 0;
        let current = start, l = 0;
        do {
            this._items[l++] = current;
            current = current.next;
        } while (start !== current);
        this._len = l;
    }
}

const tmpIt = new NodeListIterator<any, any>();
/**
 * A Fibonacci heap data structure with a key and optional value.
*/
export class FibonacciHeap<K, V> {
    private _minNode: Node<K, V> | null = null;
    private _nodeCount: number = 0;
    private _compare: CompareFunction<K, V>;

    constructor(
        compare?: CompareFunction<K, V>
    ) {
        this._compare = compare ? compare : this._defaultCompare;
    }

    /**
   * Clears the heap's data, making it an empty heap.
   */
    public clear(): void {
        this._minNode = null;
        this._nodeCount = 0;
    }

    /**
   * Decreases a key of a node.
   * @param node The node to decrease the key of.
   * @param newKey The new key to assign to the node.
   */
    public decreaseKey(node: Node<K, V>, newKey: K): void {
        if (!node) {
            throw new Error('Cannot decrease key of non-existent node');
        }
        if (this._compare({ key: newKey }, { key: node.key }) > 0) {
            throw new Error('New key is larger than old key');
        }

        node.key = newKey;
        const parent = node.parent;
        if (parent && this._compare(node, parent) < 0) {
            this._cut(node, parent, <Node<K, V>> this._minNode);
            this._cascadingCut(parent, <Node<K, V>> this._minNode);
        }
        if (this._compare(node, <Node<K, V>> this._minNode) < 0) {
            this._minNode = node;
        }
    }

    /**
   * Deletes a node.
   * @param node The node to delete.
   */
    public delete(node: Node<K, V>): void {
    // This is a special implementation of decreaseKey that sets the argument to
    // the minimum value. This is necessary to make generic keys work, since there
    // is no MIN_VALUE constant for generic types.
        const parent = node.parent;
        if (parent) {
            this._cut(node, parent, <Node<K, V>> this._minNode);
            this._cascadingCut(parent, <Node<K, V>> this._minNode);
        }
        this._minNode = node;

        this.extractMinimum();
    }

    /**
   * Extracts and returns the minimum node from the heap.
   * @return The heap's minimum node or null if the heap is empty.
   */
    public extractMinimum(): Node<K, V> | null {
        const extractedMin = this._minNode;
        if (extractedMin) {
            // Set parent to null for the minimum's children
            if (extractedMin.child) {
                let child = extractedMin.child;
                do {
                    child.parent = null;
                    child = child.next;
                } while (child !== extractedMin.child);
            }

            let nextInRootList = null;
            if (extractedMin.next !== extractedMin) {
                nextInRootList = extractedMin.next;
            }
            // Remove min from root list
            this._removeNodeFromList(extractedMin);
            this._nodeCount--;

            // Merge the children of the minimum node with the root list
            this._minNode = this._mergeLists(nextInRootList, extractedMin.child);
            if (this._minNode) {
                this._minNode = this._consolidate(this._minNode);
            }
        }
        return extractedMin;
    }

    /**
   * Returns the minimum node from the heap.
   * @return The heap's minimum node or null if the heap is empty.
   */
    public findMinimum(): Node<K, V> | null {
        return this._minNode;
    }

    /**
   * Inserts a new key-value pair into the heap.
   * @param key The key to insert.
   * @param value The value to insert.
   * @return node The inserted node.
   */
    public insert(key: K, value?: V): Node<K, V> {
        const node = new Node(key, value);
        this._minNode = this._mergeLists(this._minNode, node);
        this._nodeCount++;
        return node;
    }

    /**
   * @return Whether the heap is empty.
   */
    public isEmpty(): boolean {
        return this._minNode === null;
    }

    /**
   * @return The size of the heap.
   */
    public size(): number {
        if (this._minNode === null) {
            return 0;
        }
        return this._getNodeListSize(this._minNode);
    }

    /**
   * Joins another heap to this heap.
   * @param other The other heap.
   */
    public union(other: FibonacciHeap<K, V>): void {
        this._minNode = this._mergeLists(this._minNode, other._minNode);
        this._nodeCount += other._nodeCount;
    }

    /**
   * Compares two nodes with each other.
   * @param a The first key to compare.
   * @param b The second key to compare.
   * @return -1, 0 or 1 if a < b, a == b or a > b respectively.
   */
    private _defaultCompare(a: INode<K, V>, b: INode<K, V>): number {
        if (a.key > b.key) {
            return 1;
        }
        if (a.key < b.key) {
            return -1;
        }
        return 0;
    }

    /**
   * Cut the link between a node and its parent, moving the node to the root list.
   * @param node The node being cut.
   * @param parent The parent of the node being cut.
   * @param minNode The minimum node in the root list.
   * @return The heap's new minimum node.
   */
    private _cut(node: Node<K, V>, parent: Node<K, V>, minNode: Node<K, V>): Node<K, V> | null {
        node.parent = null;
        parent.degree--;
        if (node.next === node) {
            parent.child = null;
        } else {
            parent.child = node.next;
        }
        this._removeNodeFromList(node);
        const newMinNode = this._mergeLists(minNode, node);
        node.isMarked = false;
        return newMinNode;
    }

    /**
   * Perform a cascading cut on a node; mark the node if it is not marked,
   * otherwise cut the node and perform a cascading cut on its parent.
   * @param node The node being considered to be cut.
   * @param minNode The minimum node in the root list.
   * @return The heap's new minimum node.
   */
    private _cascadingCut(node: Node<K, V>, minNode: Node<K, V> | null): Node<K, V> | null {
        const parent = node.parent;
        if (parent) {
            if (node.isMarked) {
                minNode = this._cut(node, parent, <Node<K, V>>minNode);
                minNode = this._cascadingCut(parent, minNode);
            } else {
                node.isMarked = true;
            }
        }
        return minNode;
    }

    /**
   * Merge all trees of the same order together until there are no two trees of
   * the same order.
   * @param minNode The current minimum node.
   * @return The new minimum node.
   */
    private _consolidate(minNode: Node<K, V>): Node<K, V> | null {

        const aux = [];
        tmpIt.reset(minNode);
        while (tmpIt.hasNext()) {
            let current = tmpIt.next();

            // If there exists another node with the same degree, merge them
            let auxCurrent = aux[current.degree];
            while (auxCurrent) {
                if (this._compare(current, auxCurrent) > 0) {
                    const temp = current;
                    current = auxCurrent;
                    auxCurrent = temp;
                }
                this._linkHeaps(auxCurrent, current);
                aux[current.degree] = null;
                current.degree++;
                auxCurrent = aux[current.degree];
            }

            aux[current.degree] = current;
        }

        let newMinNode = null;
        for (let i = 0; i < aux.length; i++) {
            const node = aux[i];
            if (node) {
                // Remove siblings before merging
                node.next = node;
                node.prev = node;
                newMinNode = this._mergeLists(newMinNode, node);
            }
        }
        return newMinNode;
    }

    /**
   * Removes a node from a node list.
   * @param node The node to remove.
   */
    private _removeNodeFromList(node: Node<K, V>): void {
        const prev = node.prev;
        const next = node.next;
        prev.next = next;
        next.prev = prev;
        node.next = node;
        node.prev = node;
    }

    /**
   * Links two heaps of the same order together.
   *
   * @private
   * @param max The heap with the larger root.
   * @param min The heap with the smaller root.
   */
    private _linkHeaps(max: Node<K, V>, min: Node<K, V>): void {
        this._removeNodeFromList(max);
        min.child = this._mergeLists(max, min.child);
        max.parent = min;
        max.isMarked = false;
    }

    /**
   * Merge two lists of nodes together.
   *
   * @private
   * @param a The first list to merge.
   * @param b The second list to merge.
   * @return The new minimum node from the two lists.
   */
    private _mergeLists(a: Node<K, V> | null, b: Node<K, V> | null): Node<K, V> | null {
        if (!a) {
            if (!b) {
                return null;
            }
            return b;
        }
        if (!b) {
            return a;
        }

        const temp = a.next;
        a.next = b.next;
        a.next.prev = a;
        b.next = temp;
        b.next.prev = b;

        return this._compare(a, b) < 0 ? a : b;
    }

    /**
   * Gets the size of a node list.
   * @param node A node within the node list.
   * @return The size of the node list.
   */
    private _getNodeListSize(node: Node<K, V>): number {
        let count = 0;
        let current = node;

        do {
            count++;
            if (current.child) {
                count += this._getNodeListSize(current.child);
            }
            current = current.next;
        } while (current !== node);

        return count;
    }
}
