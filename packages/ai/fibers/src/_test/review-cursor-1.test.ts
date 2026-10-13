/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  createFiberHost,
  createMemoryRows,
  FIBER_LIMIT_HARD_BOUNDS,
  type FiberTable,
  fiberCommandId,
  fiberEntryId,
  fiberId,
  fiberOwnerId,
} from '../index.js';

const id = fiberId('story-1');

/** Review cursor-1. A fiber host must see another host's writes on the same rows, as the harness host does. */
describe('review cursor-1: fiber hosts share rows', () => {
  it("a host that already loaded the fiber sees another host's turn, and its own write is not a conflict", async () => {
    const rows = createMemoryRows<FiberTable>();
    const a = createFiberHost({ rows });
    const b = createFiberHost({ rows });
    await a.create({
      commandId: fiberCommandId('create'),
      fiberId: id,
      ownerId: fiberOwnerId('me'),
      limits: FIBER_LIMIT_HARD_BOUNDS,
    });
    expect(await b.window({ fiberId: id, limit: 8, direction: 'backward' })).toBeTruthy();
    await a.begin({
      commandId: fiberCommandId('begin:0'),
      fiberId: id,
      entryId: fiberEntryId('turn-0'),
      descriptor: { codec: 'json', version: 1, value: { prompt: 'p0' } },
    });
    const seen = await b.window({ fiberId: id, limit: 8, direction: 'backward' });
    let write = 'accepted';
    try {
      await b.begin({
        commandId: fiberCommandId('begin:1'),
        fiberId: id,
        entryId: fiberEntryId('turn-1'),
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
