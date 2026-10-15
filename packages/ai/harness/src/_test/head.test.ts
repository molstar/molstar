/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import fc from 'fast-check';
import type { HeadRecord, Segment } from '../contracts/records.js';
import { buildMove, emptyHead, entriesOf } from '../head/index.js';

/** A segment cache, as a runner keeps one. */
function cache() {
  const segments = new Map<string, Segment>();
  const read = (ids: readonly string[]) =>
    ids.map((id) => {
      const segment = segments.get(id);
      if (segment === undefined) throw new Error(`missing segment ${id}`);
      return segment;
    });
  return { segments, read };
}

describe('the head scheme', () => {
  it('reads back like a flat list after any sequence of appends and edits', async () => {
    await fc.assert(
      fc.asyncProperty(
        fc.array(
          fc.oneof(
            fc.record({ append: fc.array(fc.nat(50), { maxLength: 9 }) }),
            fc.record({ keep: fc.array(fc.boolean(), { maxLength: 60 }), add: fc.array(fc.nat(50), { maxLength: 4 }) }),
          ),
          { maxLength: 30 },
        ),
        async (operations) => {
          const { segments, read } = cache();
          let head: HeadRecord = await emptyHead();
          let flat: string[] = [];
          for (const op of operations) {
            const edit: { append: string[] } | { entries: string[] } =
              'append' in op
                ? { append: op.append.map(String) }
                : { entries: [...flat.filter((_, i) => op.keep[i] ?? true), ...op.add.map(String)] };
            const move = await buildMove({ pointer: 'main', from: head, edit, read });
            for (const segment of move.segments) segments.set(segment.id, segment);
            flat = 'append' in edit ? [...flat, ...edit.append] : edit.entries;
            head = move.head;
            expect(move.from).toBeDefined();
            expect(entriesOf(head, read)).toEqual(flat);
            expect(head.length).toBe(flat.length);
          }
        },
      ),
      { numRuns: 100 },
    );
  });

  it('keeps about log n segments and copies about n log n hashes over 2,000 appends', async () => {
    const { segments, read } = cache();
    let head = await emptyHead();
    let copied = 0;
    for (let i = 0; i < 2_000; i += 1) {
      const move = await buildMove({ pointer: 'main', from: head, edit: { append: [`e${i}`] }, read });
      for (const segment of move.segments) {
        segments.set(segment.id, segment);
        copied += segment.entries.length;
      }
      head = move.head;
      expect(head.segments.length).toBeLessThanOrEqual(Math.ceil(Math.log2(i + 2)) + 1);
    }
    expect(copied).toBeLessThanOrEqual(2_000 * (Math.ceil(Math.log2(2_000)) + 1));
  });

  it('is content-addressed: the same head built twice has the same id', async () => {
    const { read } = cache();
    const a = await buildMove({ pointer: 'main', edit: { entries: ['x', 'y', 'z'] }, read });
    const b = await buildMove({ pointer: 'main', edit: { entries: ['x', 'y', 'z'] }, read });
    expect(a.head.id).toBe(b.head.id);
    expect((await emptyHead()).id).toBe((await emptyHead()).id);
  });
});
