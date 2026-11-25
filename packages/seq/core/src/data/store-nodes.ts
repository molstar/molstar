/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { KeyedState } from '@molstar/core/reactive';
import { nodesEqual } from '@molstar/core/repr/tree';
import { nodeKey } from './keys.js';
import type { DynamicNode } from './protocol.js';
import type { NodeEntry } from './store-types.js';

// The data store's nodes (_spec/README.md "Data", Store, Nodes): held, not cached, by each data entry that carried
// them, by their source's coverage and by `hold`; with no holder they go quietly.

/** What the nodes ask of their store: a new version, a conflict said, whether a key is an old reading's. */
export interface NodesHost {
  /** A node came or went: the store's and the nodes' versions move to `version`. */
  bump(version: number): void;
  conflict(source: string): void;
  takeStale(key: string): boolean;
}

export class StoreNodes {
  readonly entries = new KeyedState<string, NodeEntry>();
  private readonly _holders = new Map<string, number>();

  constructor(
    private readonly host: NodesHost,
    private readonly next: () => number,
  ) {}

  hold(k: string): void {
    this._holders.set(k, (this._holders.get(k) ?? 0) + 1);
  }

  /**
   * Holds the nodes an answer (or a coverage) brought; the same key with other content keeps the held node, or, for a
   * new reading (`replace`), takes the new content under its identity.
   */
  intake(source: string, nodes: readonly DynamicNode[] | undefined, replace = false): readonly string[] {
    if (!nodes?.length) return [];
    const keys: string[] = [];
    const seen = new Set<string>();
    for (const n of nodes) {
      const k = nodeKey(source, n.key);
      if (seen.has(k)) continue;
      seen.add(k);
      keys.push(k);
      this.hold(k);
      const held = this.entries.get(k);
      // a key an old reading's tile brought, served again: the new reading's content under its identity
      const fresh = this.host.takeStale(k);
      if (held && nodesEqual(held.node, n.node)) continue;
      if (!held || replace || fresh) {
        const version = this.next();
        this.entries.set(k, { source, key: n.key, node: n.node, version });
        this.host.bump(version);
      } else this.host.conflict(source);
    }
    return keys;
  }

  /** A node a new reading no longer brings goes now; its holders' releases count down to nothing after it. */
  remove(k: string): void {
    if (!this.entries.has(k)) return;
    this.entries.delete(k);
    this.host.bump(this.next());
  }

  /** A node with no holder left goes, quietly. */
  release(k: string): void {
    const n = (this._holders.get(k) ?? 0) - 1;
    if (n > 0) {
      this._holders.set(k, n);
      return;
    }
    this._holders.delete(k);
    if (this.entries.has(k)) {
      this.entries.delete(k);
      this.host.bump(this.next());
    }
  }
}
