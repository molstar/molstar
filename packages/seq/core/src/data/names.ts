/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { abortError, type FindNamesOptions, isAbortError, type NameHit, type OpenedSource } from './protocol.js';

// Finding items by name across the opened sources (_spec/README.md "Data", Names; EC-228): a function over the sources a
// viewer's source sync holds, not the store's.

/** An opened source and the address its spec gives it. */
export interface NamedSource {
  readonly address: string;
  readonly source: OpenedSource;
}

/** A hit and the address of the source that found it (where `(item :id …)` addresses it). */
export interface SourceNameHit extends NameHit {
  readonly source: string;
}

/** Hits a find gives at most, all sources together, unless asked otherwise. */
export const DefaultMaxNameHits = 20;

/**
 * Asks every source that serves `findNames` (in parallel) and merges their hits in source order, each source's as it
 * gave them (its exact hits first), at most `options.max` (20) together. A source that fails is left out of the
 * answer, never failing it; an aborted `signal` rejects with an `AbortError`. Blank text finds nothing.
 */
export async function findNames(
  sources: Iterable<NamedSource>,
  text: string,
  options: FindNamesOptions,
  signal: AbortSignal,
): Promise<SourceNameHit[]> {
  if (signal.aborted) throw abortError();
  const query = text.trim();
  const max = Math.max(0, Math.floor(options.max ?? DefaultMaxNameHits));
  if (!query || max === 0) return [];
  const asked = { max, prefix: options.prefix ?? false };
  const pending: Promise<readonly SourceNameHit[]>[] = [];
  for (const { address, source } of sources) {
    if (!source.findNames) continue;
    pending.push(
      source.findNames(query, asked, signal).then(
        (hits) => hits.slice(0, max).map((h) => ({ ...h, source: address })),
        (e) => {
          if (isAbortError(e) || signal.aborted) throw abortError();
          return [];
        },
      ),
    );
  }
  const all = await Promise.all(pending);
  if (signal.aborted) throw abortError();
  const out: SourceNameHit[] = [];
  for (const hits of all) {
    for (const h of hits) {
      if (out.length >= max) return out;
      out.push(h);
    }
  }
  return out;
}
