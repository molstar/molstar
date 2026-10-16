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
  type ThreadHost,
  type ThreadPayloadEnvelope,
  type ThreadRows,
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
const json = (value: unknown): ThreadPayloadEnvelope => ({ codec: 'json', version: 1, value });
const create = (host: ThreadHost) =>
  host.create({
    commandId: threadCommandId('create'),
    threadId: id,
    ownerId: threadOwnerId('me'),
    limits: THREAD_LIMIT_HARD_BOUNDS,
  });

async function turn(host: ThreadHost, n: number, settle = true) {
  const entryId = threadEntryId(`turn-${n}`);
  await host.begin({
    commandId: threadCommandId(`begin:${n}`),
    threadId: id,
    entryId,
    descriptor: json({ prompt: `p${n}` }),
  });
  if (settle)
    await host.settle({
      commandId: threadCommandId(`settle:${n}`),
      threadId: id,
      entryId,
      payload: json({ reply: `r${n}` }),
    });
}

// The host's own tests (from Stories' stories-threads), over the memory rows; the backends run the conformance.
describe('thread host: durable', () => {
  it('keeps a thread across hosts on the same rows, loading it only when used', async () => {
    const memory = createMemoryRows<ThreadTable>();
    const first = createThreadHost({ rows: memory });
    await create(first);
    for (let n = 0; n < 5; n++) await turn(first, n);
    await turn(first, 5, false); // interrupted: begun, never settled

    let loads = 0;
    const counting: ThreadRows = {
      ...memory,
      load: (scope) => {
        loads++;
        return memory.load(scope);
      },
    };
    const second = createThreadHost({ rows: counting });
    expect(loads).toBe(0);
    const window = read(await second.window({ threadId: id, limit: 3, direction: 'backward' }));
    expect(loads).toBe(1);
    expect(window?.entries.map((e) => [e.position, e.descriptor.value, e.state.kind])).toEqual([
      [3, { prompt: 'p3' }, 'settled'],
      [4, { prompt: 'p4' }, 'settled'],
      [5, { prompt: 'p5' }, 'streaming'],
    ]);
    // The interrupted turn settles in the new host, and a replayed command is recognized.
    const settle = {
      commandId: threadCommandId('settle:5'),
      threadId: id,
      entryId: threadEntryId('turn-5'),
      payload: json({ reply: 'r5' }),
    };
    const settled = await second.settle(settle);
    expect(settled.status).toBe('accepted');
    expect(await second.settle(settle)).toEqual(settled);
    const again = createThreadHost({ rows: memory });
    const last = read(await again.window({ threadId: id, limit: 1, direction: 'backward' }));
    expect(last?.entries[0]?.state).toMatchObject({ kind: 'settled', payload: { value: { reply: 'r5' } } });
  });

  it('changes nothing when the write fails', async () => {
    const rows = createMemoryRows<ThreadTable>();
    let fail = false;
    const flaky: ThreadRows = {
      ...rows,
      apply: (...a) => (fail ? Promise.reject(new Error('disk full')) : rows.apply(...a)),
    };
    const host = createThreadHost({ rows: flaky });
    await create(host);
    await turn(host, 0);
    fail = true;
    await expect(turn(host, 1)).rejects.toThrow('disk full');
    fail = false;
    const window = read(await host.window({ threadId: id, limit: 8, direction: 'backward' }));
    expect(window?.entries.map((e) => e.entryId)).toEqual(['turn-0']);
    await turn(host, 1);
    expect(read(await host.window({ threadId: id, limit: 8, direction: 'backward' }))?.entries).toHaveLength(2);
  });

  it('forgets a thread with its rows', async () => {
    const rows = createMemoryRows<ThreadTable>();
    const host = createThreadHost({ rows });
    await create(host);
    await turn(host, 0);
    await host.forget(id);
    expect(await rows.load('thread:story-1')).toBeUndefined();
    expect(await host.window({ threadId: id, limit: 8, direction: 'backward' })).toEqual({
      status: 'not-found',
      reason: 'thread-missing',
    });
  });

  it('tells a subscriber which thread a command changed, and a watcher which one another handle changed', async () => {
    const changed: string[] = [];
    const heard: string[] = [];
    let announce: ((scope: string) => void) | undefined;
    let watching = 0;
    const rows: ThreadRows = {
      ...createMemoryRows<ThreadTable>(),
      watch(listener) {
        watching++;
        announce = listener;
        return () => {
          watching--;
          announce = undefined;
        };
      },
    };
    const host = createThreadHost({ rows });
    const stop = host.subscribe((thread) => changed.push(thread));
    expect(watching).toBe(0); // nothing is watched until someone watches
    await create(host);
    const unwatch = host.watch((thread) => heard.push(thread));
    await Promise.resolve();
    await Promise.resolve();
    expect(watching).toBe(1);
    announce?.('thread:story-1');
    announce?.('somebody-else');
    unwatch();
    expect(watching).toBe(0);
    stop();
    await turn(host, 0);
    expect({ changed, heard }).toEqual({ changed: ['story-1'], heard: ['story-1'] });
  });
});
