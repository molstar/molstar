/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A jump that waits for a node data makes (_spec/README.md "AddressModel", A focus that waits): the node is asked of
// its source, never guessed. A source that finds names (`findNames`) says where its key is, and the region there is
// read through one wanter, so the node arrives even when nothing shows its region (it was evicted, or never read); a
// source that finds none of the key says it has no such node, once a reading. A source that finds no names leaves the
// jump waiting for its data, as before. Nothing is held while nothing waits: a waiting ref holds no data.

import { type DataStore, nodeKey, type OpenedSource, type Want, type Wanter } from '@molstar/seqx-core/data';
import type { WaitedNode } from '@molstar/seqx-core/query';

export class NodeSeek {
  private _wanter: Wanter | undefined;
  private _abort: AbortController | undefined;
  /** The nodes asked for the jump that waits now, by `<source>#<key>`. */
  private readonly _sought = new Set<string>();
  private _wants: Want[] = [];
  /** The keys a source said it has no node of, with the reading that said it (a new reading is asked again). */
  private readonly _unknown = new Map<string, OpenedSource>();

  constructor(
    private readonly store: DataStore,
    /** Told once a reading when a source says it has no node of a key. */
    private readonly refused: (node: WaitedNode) => void,
  ) {}

  private _source(address: string): OpenedSource | undefined {
    for (const s of this.store.named()) if (s.address === address) return s.source;
    return undefined;
  }

  /** Whether the node's source, as attached now, said it has no node of its key. */
  unknown(node: WaitedNode): boolean {
    const by = this._unknown.get(nodeKey(node.source, node.key));
    return by !== undefined && by === this._source(node.source);
  }

  /**
   * Asks for the first of `nodes` not held (one held inside it is asked once it arrives: its source is that node): its
   * source finds where its key is, and that region is read. Each node is asked once a wait.
   */
  seek(nodes: readonly WaitedNode[]): void {
    const n = nodes.find((x) => !this.store.getNode(x.source, x.key));
    if (!n) return;
    const k = nodeKey(n.source, n.key);
    const source = this._source(n.source);
    if (this._sought.has(k) || !source?.findNames) return;
    this._sought.add(k);
    this._abort ??= new AbortController();
    const signal = this._abort.signal;
    source.findNames(n.key, { max: 20 }, signal).then(
      (hits) => {
        if (signal.aborted) return;
        const hit = hits.find((h) => h.id === n.key);
        if (!hit) {
          this._unknown.set(k, source);
          this.refused(n);
          return;
        }
        // the finest level: a region's answer brings its nodes there
        const level = this.store.coverage(n.source)?.levels[0] ?? 1;
        this._wants = [...this._wants, { source: n.source, system: hit.system, start: hit.start, end: hit.end, level }];
        (this._wanter ??= this.store.wanter()).want(this._wants);
      },
      // a source that fails to find leaves the jump waiting for its data
      () => {},
    );
  }

  /** The wait ended (the jump ran, another replaced it, or it was dropped): what it asked for is let go. */
  stop(): void {
    this._abort?.abort();
    this._abort = undefined;
    this._sought.clear();
    if (this._wants.length) {
      this._wants = [];
      this._wanter?.want([]);
    }
  }

  dispose(): void {
    this.stop();
    this._wanter?.dispose();
    this._wanter = undefined;
  }
}
