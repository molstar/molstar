/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  createThreadHost,
  createMemoryRows,
  THREAD_LIMIT_HARD_BOUNDS,
  type ThreadFailure,
  type ThreadTable,
  threadCommandId,
  threadEntryId,
  threadId,
  threadOwnerId,
} from '../index.js';

/** A read's page; a failure fails the case, naming why. */
function read<T extends object>(result: T | ThreadFailure): T {
  if ('status' in result && 'reason' in result) throw new Error(`read failed: ${result.status} ${result.reason}`);
  return result as T;
}

const id = threadId('story-1');

/** Review cursor-1. A thread host must see another host's writes on the same rows, as the harness host does. */
describe('review cursor-1: thread hosts share rows', () => {
  it("a host that already loaded the thread sees another host's turn, and its own write is not a conflict", async () => {
    const rows = createMemoryRows<ThreadTable>();
    const a = createThreadHost({ rows });
    const b = createThreadHost({ rows });
    await a.create({
      commandId: threadCommandId('create'),
      threadId: id,
      ownerId: threadOwnerId('me'),
      limits: THREAD_LIMIT_HARD_BOUNDS,
    });
    expect(read(await b.window({ threadId: id, limit: 8, direction: 'backward' }))).toBeTruthy();
    await a.begin({
      commandId: threadCommandId('begin:0'),
      threadId: id,
      entryId: threadEntryId('turn-0'),
      descriptor: { codec: 'json', version: 1, value: { prompt: 'p0' } },
    });
    const seen = read(await b.window({ threadId: id, limit: 8, direction: 'backward' }));
    let write = 'accepted';
    try {
      await b.begin({
        commandId: threadCommandId('begin:1'),
        threadId: id,
        entryId: threadEntryId('turn-1'),
        descriptor: { codec: 'json', version: 1, value: { prompt: 'p1' } },
      });
    } catch (e) {
      write = e instanceof Error ? e.message : String(e);
    }
    expect({ entries: seen?.entries.map((e) => e.entryId), write }).toEqual({
      entries: ['turn-0'],
      write: 'accepted',
    });
  });
});
