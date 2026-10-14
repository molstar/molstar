/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The prefetch story's document, provider and request log (`Seq*/Streaming/Prefetch`), kept out of the story file so
// the smoke test can open the document. The chromosome and its providers are the Streamed data story's (a committed
// synthetic one, offline); the log records every region the viewer asks of them, with the priority it asked it at.
// No React here.

import type { Loaders, StreamRequest } from "@molstar/seq-core/format";
import { STREAMED_SECTION, streamedDoc, streamedLoaders } from "./seqstar-docs-streamed";

export const PREFETCH_SECTION = STREAMED_SECTION;
/** The snapshots: a title, and the 100 kb of the chromosome it focuses on (far apart, so one's data is not the next's). */
export const PREFETCH_SNAPSHOTS: readonly (readonly [string, number, number])[] = [
  ["A", 2_000_000, 2_100_000],
  ["B", 9_000_000, 9_100_000],
  ["C", 16_000_000, 16_100_000],
];
export const prefetchDoc = () => streamedDoc({ snapshots: PREFETCH_SNAPSHOTS });

/** One request as the log shows it. */
export interface Asked {
  readonly n: number;
  readonly source: string;
  readonly level: number;
  readonly start: number;
  readonly stop: number;
  readonly priority: StreamRequest["priority"];
}
/** The requests in the order made, newest kept `KEEP` of, observable for a panel. */
export class AskedLog {
  static readonly KEEP = 200;
  entries: readonly Asked[] = [];
  private count = 0;
  private readonly listeners = new Set<() => void>();
  push(source: string, r: StreamRequest) {
    const next = [
      ...this.entries,
      { n: ++this.count, source, level: r.level, start: r.start, stop: r.stop, priority: r.priority },
    ];
    this.entries = next.slice(-AskedLog.KEEP);
    for (const f of this.listeners) f();
  }
  clear() {
    this.entries = [];
    for (const f of this.listeners) f();
  }
  readonly subscribe = (f: () => void) => {
    this.listeners.add(f);
    return () => this.listeners.delete(f);
  };
  readonly snapshot = () => this.entries;
}

/** The streamed chromosome's loaders, each request logged as it is made. */
export const prefetchLoaders = (log: AskedLog): Loaders => streamedLoaders((source, r) => log.push(source, r));
