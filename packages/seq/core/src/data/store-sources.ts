/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import { channelOf, channelSource, fieldChannel } from './channels.js';
import type { OpenedSource } from './protocol.js';
import { checkCoverage } from './seq-format.js';
import { type SummarySettings, summaryBytes, withSummaries } from './summaries.js';

// The data store's sources (_spec/README.md "Data", Store): each attached source checked and summarized, its systems'
// names, the nodes its coverage brought; a source's values channels made when a tile of one is read (Channels); and
// what a new reading (Reattach) leaves of the old one's tiles: the keys they brought, kept by identity while held and
// replaced when served again, removed when the tile that brought one is read without it.

export interface AttachedSource {
  readonly opened: OpenedSource;
  /** As attached, before its summaries: what a values channel of it reads. */
  readonly raw: OpenedSource;
  readonly names: ReadonlyMap<string, string>;
  /** The node keys its coverage brought, held while it is attached. */
  readonly nodes: readonly string[];
}

export type PreparedSource = Omit<AttachedSource, 'nodes'>;

export class StoreSources {
  private readonly _all = new Map<string, AttachedSource>();
  /** A node key an old reading's tiles brought → the tiles (data keys) that brought it, until read again. */
  private readonly _stale = new Map<string, { readonly source: string; readonly tiles: Set<string> }>();

  constructor(private readonly settings: { readonly summaries: SummarySettings | false; readonly maxRows: number }) {}

  /** A source checked and summarized, its systems' names; else the coverage's issue. */
  prepare(source: OpenedSource): PreparedSource | Issue {
    const bad = checkCoverage(source.coverage);
    if (bad) return bad;
    // a level-1 source's values get a pyramid, so a values layer zoomed out reads a coarse level
    const summaries = this.settings.summaries;
    const opened = summaries ? withSummaries(source, summaries, this.settings.maxRows) : source;
    const names = new Map<string, string>();
    for (const s of opened.coverage.systems) {
      names.set(s.name, s.name);
      for (const a of s.aliases ?? []) if (!names.has(a)) names.set(a, s.name);
    }
    return { opened, raw: source, names };
  }

  /** An attached source by its own address (no channel). */
  attached(address: string): AttachedSource | undefined {
    return this._all.get(address);
  }

  /** Every attached source by its own address (no channel). */
  all(): IterableIterator<[string, AttachedSource]> {
    return this._all.entries();
  }

  set(address: string, source: AttachedSource): void {
    this._all.set(address, source);
  }

  /**
   * What an address reads: an attached source, or a channel of one. Read only (`make` false), a channel answers with its
   * source's names and coverage, made of nothing; read for a tile (`make`), it is made once, summarized on its own, and
   * goes with its source.
   */
  get(address: string, make = false): AttachedSource | undefined {
    const at = this._all.get(address);
    if (at) return at;
    const channel = channelOf(address);
    const base = channel && this._all.get(channel.source);
    if (!channel || !base || channelOf(channel.source)) return undefined;
    if (!make) return base;
    const prepared = this.prepare(fieldChannel(base.raw, channel.field)) as PreparedSource;
    const made = { ...prepared, nodes: [] };
    this._all.set(address, made);
    return made;
  }

  /** Removes a source and its channels (theirs disposed; the source's is the caller's). */
  delete(address: string): void {
    this._all.delete(address);
    this.dropChannels(address);
    for (const [k, s] of [...this._stale]) if (s.source === address) this._stale.delete(k);
  }

  /** Disposes and removes a source's channels. */
  dropChannels(address: string): void {
    for (const [a, s] of [...this._all])
      if (a !== address && channelSource(a) === address) {
        this._all.delete(a);
        s.opened.dispose();
      }
  }

  addresses(): string[] {
    return [...this._all.keys()];
  }

  get summaryBytes(): number {
    let n = 0;
    for (const s of this._all.values()) n += summaryBytes(s.opened);
    return n;
  }

  /** Marks a key an old reading's tile brought (`tile`), not in the new reading's coverage: replaced when served. */
  markStale(source: string, key: string, tile: string): void {
    let at = this._stale.get(key);
    if (!at) this._stale.set(key, (at = { source, tiles: new Set() }));
    at.tiles.add(tile);
  }

  /** Whether a key comes from an old reading: served again, it takes the new content under its identity. */
  takeStale(key: string): boolean {
    return this._stale.delete(key);
  }

  /** The stale keys a tile read again without them no longer brings: when no other tile may, they go. */
  goneWith(tile: string, brought: ReadonlySet<string>): string[] {
    const out: string[] = [];
    for (const [k, { tiles }] of this._stale) {
      if (!tiles.has(tile) || brought.has(k)) continue;
      tiles.delete(tile);
      if (!tiles.size) {
        this._stale.delete(k);
        out.push(k);
      }
    }
    return out;
  }
}
