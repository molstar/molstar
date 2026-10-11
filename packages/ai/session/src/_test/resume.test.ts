/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createMemoryStorage } from '@molstar/ai-components';
import type { Clock, HarnessStorage } from '@molstar/ai-harness';
import { createToolkit } from '@molstar/ai-tools';
import type { TurnEvent } from '../index.js';
import { callOf, fixture, realClock, tool, toy } from './helpers.js';

/** Storage that stops answering, like a tab that just died: nothing it does reaches the records any more. */
function killable(inner: HarnessStorage) {
  let dead = false;
  const never = new Promise<never>(() => {});
  return {
    kill: () => {
      dead = true;
    },
    dead: () => dead,
    storage: {
      transact: (h, run) => (dead ? never : inner.transact(h, run)),
      putBlobs: (h, blobs) => (dead ? never : inner.putBlobs(h, blobs)),
      getBlobs: (h, refs) => (dead ? never : inner.getBlobs(h, refs)),
      watch: (h) => inner.watch(h),
      drop: (h) => (dead ? never : inner.drop(h)),
    } satisfies HarnessStorage,
  };
}

describe('resume after the tab dies mid-turn', () => {
  it('a new runtime reads the harness and finishes the interrupted turn', async () => {
    const records = createMemoryStorage();
    // Real time, shifted when the test moves past the old lease; in the dead tab, time no longer passes for it.
    let offset = 0;
    const tab1 = killable(records);
    const clock = (alive: () => boolean): Clock => ({
      now: () => realClock.now() + offset,
      sleep: (ms, signal) => (alive() ? realClock.sleep(ms, signal) : new Promise<never>(() => {})),
    });
    const host = toy();

    // Tab 1: the tool never returns; the tab dies while it waits.
    let started = false;
    const hanging = createToolkit<typeof host>({
      tools: [
        tool({
          id: 'slow_read',
          description: 'Reads, slowly.',
          input: { kind: 'object', fields: {} },
          modules: ['log'],
          read: true,
          run: () => {
            started = true;
            return new Promise(() => {});
          },
        }),
      ],
      core: ['slow_read'],
      groups: [],
    });
    const one = fixture(
      [{ calls: [callOf('slow_read')] }],
      {
        toolkit: hanging,
        clock: clock(() => !tab1.dead()),
        leaseMs: 1_000,
        harnessId: 'h-1',
        storage: tab1.storage,
      },
      {},
      host,
    );
    void one.session.turn({ prompt: 'Read it.' });
    for (let i = 0; i < 100 && !started; i++) await new Promise((r) => setTimeout(r, 10));
    expect(started).toBe(true);
    tab1.kill();

    // Tab 2, after the old lease expired: a fresh runtime over the same records.
    offset += 3_000;
    let calls = 0;
    const working = createToolkit<typeof host>({
      tools: [
        tool({
          id: 'slow_read',
          description: 'Reads, slowly.',
          input: { kind: 'object', fields: {} },
          modules: ['log'],
          read: true,
          run: ({ call }) => {
            calls += 1;
            host.log.lines.push(`read ${call.prior}`);
            return { ok: true, did: 'Read.' };
          },
        }),
      ],
      core: ['slow_read'],
      groups: [],
    });
    const two = fixture(
      [{ text: 'Read it.' }],
      { toolkit: working, clock: clock(() => true), leaseMs: 1_000, harnessId: 'h-1', storage: records },
      {},
      host,
    );
    expect(await two.session.interrupted()).toBe(true);
    const turns: TurnEvent[] = [];
    const result = await two.session.resume({ onTurn: (e) => turns.push(e) });
    expect(result?.status).toBe('completed');
    expect(result?.text).toBe('Read it.');
    expect(calls).toBe(1);
    // The call is told an earlier attempt may have run.
    expect(host.log.lines).toEqual(['read possible']);
    // The host is told of the turn it takes up.
    expect(turns).toEqual([{ kind: 'begun', turnId: expect.any(String), resumed: true }]);

    // Nothing left to resume, and the next turn in the same harness is not a first turn.
    expect(await two.session.resume()).toBeUndefined();
    expect(await two.session.isFresh()).toBe(false);
    one.session.dispose();
    two.session.dispose();
  });
});

describe('claims before the first await', () => {
  it('a resume holds the session while it reads the harness: a second resume, a turn and a forget wait', async () => {
    const { session } = fixture();
    const first = session.resume();
    expect(session.busy).toBe(true);
    await expect(session.resume()).resolves.toBeUndefined();
    await expect(session.turn({ prompt: 'go' })).rejects.toThrow('already running');
    await expect(session.forget()).rejects.toThrow('stop it first');
    await expect(first).resolves.toBeUndefined();
    expect(session.busy).toBe(false);
  });

  it('a forget holds the session too', async () => {
    const { session } = fixture();
    const forgetting = session.forget();
    expect(session.busy).toBe(true);
    await expect(session.resume()).resolves.toBeUndefined();
    await forgetting;
    expect(session.busy).toBe(false);
  });
});
