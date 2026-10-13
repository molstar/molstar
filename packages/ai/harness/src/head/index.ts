/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { digest } from '@molstar/ai-agent';
import type { Move } from '../contracts/ports.js';
import type { EntryId, HeadRecord, ScopeId, Segment } from '../contracts/records.js';

/**
 * The shared head scheme (harness.md §14). A head is a list of segments; a segment an immutable run of entry
 * hashes. Both are content-addressed. Segments are kept in size classes like a binary counter: an
 * append adds one segment and merges it into the last while the last is no larger, so a head holds
 * about log n segments and a hash is copied about log n times over its life. An edit rebuilds the
 * head from chunks of decreasing powers of two. The runner calls this with its cache (`read`).
 */

async function segmentOf(entries: readonly EntryId[]): Promise<Segment> {
  return { id: await digest({ segment: [...entries] }), entries };
}

async function headOf(segments: readonly Segment[]): Promise<HeadRecord> {
  const ids = segments.map((segment) => segment.id);
  const length = segments.reduce((sum, segment) => sum + segment.entries.length, 0);
  return { id: await digest({ head: ids }), segments: ids, length };
}

/** The empty head every harness starts with. */
export function emptyHead(): Promise<HeadRecord> {
  return headOf([]);
}

/** The head's entries, in order. */
export function entriesOf(head: HeadRecord, read: (ids: readonly string[]) => readonly Segment[]): EntryId[] {
  return read(head.segments).flatMap((segment) => segment.entries);
}

/** Splits `entries` into chunks of decreasing powers of two (the binary representation of the count). */
function chunks(entries: readonly EntryId[]): EntryId[][] {
  const out: EntryId[][] = [];
  let at = 0;
  for (let size = 2 ** Math.floor(Math.log2(Math.max(1, entries.length))); size >= 1; size /= 2)
    if (entries.length - at >= size) {
      out.push(entries.slice(at, at + size));
      at += size;
    }
  return out;
}

export async function buildMove(input: {
  readonly pointer: ScopeId;
  readonly from?: HeadRecord;
  readonly edit: { readonly append: readonly EntryId[] } | { readonly entries: readonly EntryId[] };
  read(segments: readonly string[]): readonly Segment[];
}): Promise<Move> {
  const created: Segment[] = [];
  let segments: Segment[];
  if ('append' in input.edit) {
    segments = input.from === undefined ? [] : [...input.read(input.from.segments)];
    if (input.edit.append.length > 0) {
      let tail: EntryId[] = [...input.edit.append];
      for (
        let last = segments.at(-1);
        last !== undefined && last.entries.length <= tail.length;
        last = segments.at(-1)
      ) {
        tail = [...last.entries, ...tail];
        segments.pop();
      }
      const merged = await segmentOf(tail);
      created.push(merged);
      segments.push(merged);
    }
  } else {
    segments = await Promise.all(chunks(input.edit.entries).map(segmentOf));
    created.push(...segments);
  }
  const head = await headOf(segments);
  return {
    pointer: input.pointer,
    ...(input.from !== undefined ? { from: input.from.id } : {}),
    head,
    segments: created,
  };
}
