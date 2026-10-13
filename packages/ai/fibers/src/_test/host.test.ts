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
  type FiberHost,
  type FiberPayloadEnvelope,
  type FiberRows,
  type FiberTable,
  fiberCommandId,
  fiberEntryId,
  fiberId,
  fiberOwnerId,
} from '../index.js';

const id = fiberId('story-1');
const json = (value: unknown): FiberPayloadEnvelope => ({ codec: 'json', version: 1, value });
const create = (host: FiberHost) =>
  host.create({
    commandId: fiberCommandId('create'),
    fiberId: id,
    ownerId: fiberOwnerId('me'),
    limits: FIBER_LIMIT_HARD_BOUNDS,
  });

async function turn(host: FiberHost, n: number, settle = true) {
  const entryId = fiberEntryId(`turn-${n}`);
  await host.begin({
    commandId: fiberCommandId(`begin:${n}`),
    fiberId: id,
    entryId,
    descriptor: json({ prompt: `p${n}` }),
  });
  if (settle)
    await host.settle({
      commandId: fiberCommandId(`settle:${n}`),
      fiberId: id,
      entryId,
      payload: json({ reply: `r${n}` }),
    });
}

// The host's own tests (from Stories' stories-fibers), over the memory rows; the backends run the conformance.
describe('fiber host: durable', () => {
  it('keeps a fiber across hosts on the same rows, loading it only when used', async () => {
    const memory = createMemoryRows<FiberTable>();
    const first = createFiberHost({ rows: memory });
    await create(first);
    for (let n = 0; n < 5; n++) await turn(first, n);
    await turn(first, 5, false); // interrupted: begun, never settled

    let loads = 0;
    const counting: FiberRows = {
      ...memory,
      load: (scope) => {
        loads++;
        return memory.load(scope);
      },
    };
    const second = createFiberHost({ rows: counting });
    expect(loads).toBe(0);
    const window = await second.window({ fiberId: id, limit: 3, direction: 'backward' });
    expect(loads).toBe(1);
    expect(window?.entries.map((e) => [e.position, e.descriptor.value, e.state.kind])).toEqual([
      [3, { prompt: 'p3' }, 'settled'],
      [4, { prompt: 'p4' }, 'settled'],
      [5, { prompt: 'p5' }, 'streaming'],
    ]);
    // The interrupted turn settles in the new host, and a replayed command is recognized.
    const settle = {
      commandId: fiberCommandId('settle:5'),
      fiberId: id,
      entryId: fiberEntryId('turn-5'),
      payload: json({ reply: 'r5' }),
    };
    const settled = await second.settle(settle);
    expect(settled.status).toBe('accepted');
    expect(await second.settle(settle)).toEqual(settled);
    const again = createFiberHost({ rows: memory });
    const last = await again.window({ fiberId: id, limit: 1, direction: 'backward' });
    expect(last?.entries[0]?.state).toMatchObject({ kind: 'settled', payload: { value: { reply: 'r5' } } });
  });

  it('changes nothing when the write fails', async () => {
    const rows = createMemoryRows<FiberTable>();
    let fail = false;
    const flaky: FiberRows = {
      ...rows,
      apply: (...a) => (fail ? Promise.reject(new Error('disk full')) : rows.apply(...a)),
    };
    const host = createFiberHost({ rows: flaky });
    await create(host);
    await turn(host, 0);
    fail = true;
    await expect(turn(host, 1)).rejects.toThrow('disk full');
    fail = false;
    const window = await host.window({ fiberId: id, limit: 8, direction: 'backward' });
    expect(window?.entries.map((e) => e.entryId)).toEqual(['turn-0']);
    await turn(host, 1);
    expect((await host.window({ fiberId: id, limit: 8, direction: 'backward' }))?.entries).toHaveLength(2);
  });

  it('forgets a fiber with its rows', async () => {
    const rows = createMemoryRows<FiberTable>();
    const host = createFiberHost({ rows });
    await create(host);
    await turn(host, 0);
    await host.forget(id);
    expect(await rows.load('fiber:story-1')).toBeUndefined();
    expect(await host.window({ fiberId: id, limit: 8, direction: 'backward' })).toBeUndefined();
  });

  it('tells a subscriber which fiber a command changed, and a watcher which one another handle changed', async () => {
    const changed: string[] = [];
    const heard: string[] = [];
    let announce: ((scope: string) => void) | undefined;
    let watching = 0;
    const rows: FiberRows = {
      ...createMemoryRows<FiberTable>(),
      watch(listener) {
        watching++;
        announce = listener;
        return () => {
          watching--;
          announce = undefined;
        };
      },
    };
    const host = createFiberHost({ rows });
    const stop = host.subscribe((fiber) => changed.push(fiber));
    expect(watching).toBe(0); // nothing is watched until someone watches
    await create(host);
    const unwatch = host.watch((fiber) => heard.push(fiber));
    await Promise.resolve();
    await Promise.resolve();
    expect(watching).toBe(1);
    announce?.('fiber:story-1');
    announce?.('somebody-else');
    unwatch();
    expect(watching).toBe(0);
    stop();
    await turn(host, 0);
    expect({ changed, heard }).toEqual({ changed: ['story-1'], heard: ['story-1'] });
  });
});
