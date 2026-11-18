/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// One reusable cursor for every collection (README §5): a walk of a trie on a fixed stack, or of a vector's leaves.
// `next()` allocates nothing; a cursor passed to another collection's `cursor(c)` is reset to it.

import { bitCount, KindChamp, type Root, type TrieNode } from './champ.js';
import type { Cursor, MapCursor } from './config.js';

/** Leaf arrays a vector cursor reads: the vector's `arrayFor`. */
export interface Leaves {
  readonly size: number;
  arrayFor(index: number): any[];
}

/** A removed entry's slot in an ordered collection's vector. */
export const Hole: unique symbol = Symbol('Hole') as any;

const ModeDone = 0;
const ModeTrie = 1;
/** a vector of values */
const ModeVector = 2;
/** a vector of keys with holes: an ordered set */
const ModeKeys = 3;
/** a vector of interleaved keys and values with holes: an ordered map */
const ModePairs = 4;

/** The trie's deepest path: 7 levels of 5 bits, then a collision node. */
const MaxDepth = 9;
const NoLeaf: any[] = [];

export class CursorImpl implements MapCursor<any, any>, Cursor<any> {
  key: any = undefined;
  value: any = undefined;
  private mode = ModeDone;
  private depth = -1;
  private readonly nodes: any[] = new Array(MaxDepth).fill(undefined);
  private readonly pos = new Int32Array(MaxDepth);
  private leaves: Leaves | undefined = undefined;
  private index = 0;
  private leaf: any[] = NoLeaf;

  private clear() {
    this.key = undefined;
    this.value = undefined;
    for (let i = 0; i <= this.depth && i < MaxDepth; i++) this.nodes[i] = undefined;
    this.depth = -1;
    this.leaves = undefined;
    this.leaf = NoLeaf;
    this.index = 0;
    this.mode = ModeDone;
  }

  resetTrie(root: Root) {
    this.clear();
    this.mode = ModeTrie;
    this.depth = 0;
    this.nodes[0] = root;
    this.pos[0] = 0;
  }

  resetVector(leaves: Leaves) {
    this.clear();
    this.mode = ModeVector;
    this.leaves = leaves;
  }

  resetKeys(leaves: Leaves) {
    this.clear();
    this.mode = ModeKeys;
    this.leaves = leaves;
  }

  resetPairs(leaves: Leaves) {
    this.clear();
    this.mode = ModePairs;
    this.leaves = leaves;
  }

  next(): boolean {
    switch (this.mode) {
      case ModeTrie:
        return this.nextTrie();
      case ModeVector:
      case ModeKeys:
        return this.nextVector();
      case ModePairs:
        return this.nextPair();
      default:
        return false;
    }
  }

  private element(i: number) {
    if ((i & 31) === 0 || this.leaf === NoLeaf) this.leaf = this.leaves!.arrayFor(i);
    return this.leaf[i & 31];
  }

  private nextVector(): boolean {
    const v = this.leaves!;
    while (this.index < v.size) {
      const x = this.element(this.index++);
      if (this.mode === ModeKeys && x === Hole) continue;
      this.key = x;
      this.value = x;
      return true;
    }
    this.clear();
    return false;
  }

  private nextPair(): boolean {
    const v = this.leaves!;
    while (this.index < v.size) {
      const k = this.element(this.index);
      const x = this.element(this.index + 1);
      this.index += 2;
      if (k === Hole) continue;
      this.key = k;
      this.value = x;
      return true;
    }
    this.clear();
    return false;
  }

  private nextTrie(): boolean {
    const nodes = this.nodes;
    const pos = this.pos;
    let d = this.depth;
    while (d >= 0) {
      const n: Root | TrieNode = nodes[d];
      const p = pos[d];
      if (n.kind !== KindChamp) {
        const e = n.entries;
        if (p < e.length) {
          this.key = e[p];
          this.value = e[p + 1];
          pos[d] = p + 2;
          this.depth = d;
          return true;
        }
        nodes[d--] = undefined;
        continue;
      }
      const c = n.content;
      const data = 2 * bitCount(n.dataMap);
      if (p < data) {
        this.key = c[p];
        this.value = c[p + 1];
        pos[d] = p + 2;
        this.depth = d;
        return true;
      }
      if (p < c.length) {
        pos[d] = p + 1;
        nodes[++d] = c[c.length - 1 - (p - data)];
        pos[d] = 0;
        continue;
      }
      nodes[d--] = undefined;
    }
    this.clear();
    return false;
  }
}

export function asCursor(c: unknown): CursorImpl {
  return c instanceof CursorImpl ? c : new CursorImpl();
}
