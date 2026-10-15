/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { LaunchReport, Runner } from '@molstar/ai-harness';
import { createSession, type Session } from '@molstar/ai-conformance';
import { settle } from '@molstar/ai-conformance/testing';
import { createMemoryStorage } from '../storage/memory.js';
import { createInProcessSupervisor, type InProcessSupervisorOptions } from '../supervisor/in-process.js';

const fixture = { name: 'memory', storage: () => createMemoryStorage() };
const WINDOWS = { leaseMs: 30_000, launchMs: 10_000 };

/** Lets real work (WebCrypto included) run, and virtual time move, until `done` holds. */
async function until(s: Session, done: () => boolean | Promise<boolean>, limit = 400): Promise<void> {
  for (let i = 0; i < limit; i += 1) {
    await settle();
    if (await done()) return;
    await s.clock.advance(250);
  }
  throw new Error('until: the condition never held');
}

const statuses = async (s: Session) => (await s.turns()).map((turn) => turn.status);

function supervised(s: Session, overrides: Partial<InProcessSupervisorOptions> = {}) {
  const runners: Runner[] = [];
  const supervisor = createInProcessSupervisor({
    client: s.client,
    clock: s.clock,
    startRunner: () => {
      const runner = s.runner();
      runners.push(runner);
      return runner;
    },
    ...overrides,
  });
  return { supervisor, runners };
}

describe('in-process supervisor (harness.md §4)', () => {
  it('starts a runner per launch, once per key; the report ends picked-up', async () => {
    const s = await createSession(fixture, { host: WINDOWS, model: [{ text: 'hi' }] });
    const { supervisor, runners } = supervised(s);
    const { launch } = await s.send('go', { run: false });
    const reports: LaunchReport[] = [];
    let completed = false;
    supervisor.ensure(s.harnessId, launch as string).subscribe(
      (r) => reports.push(r),
      undefined,
      () => {
        completed = true;
      },
    );
    supervisor.ensure(s.harnessId, launch as string);
    await until(s, async () => (await statuses(s))[0] === 'completed' && completed);
    expect(runners).toHaveLength(1);
    expect(reports.map((r) => r.phase)).toEqual(['launching', 'picked-up']);
    supervisor.close();
  });

  it('relaunches after a crash: the lease runs out, it checks, and a new runner finishes the turn', async () => {
    const s = await createSession(fixture, {
      host: WINDOWS,
      model: [{ text: 'lost', waitMs: 600_000 }, { text: 'found' }],
    });
    const { supervisor, runners } = supervised(s);
    const { launch } = await s.send('go', { run: false });
    supervisor.ensure(s.harnessId, launch as string);
    await until(s, () => s.model.invocations.length === 1);
    runners[0]?.abort('crash');
    await until(s, async () => (await statuses(s))[0] === 'completed');
    expect(runners).toHaveLength(2);
    supervisor.close();
  });

  it('a change while it reads (nothing to do, it found) still wakes it: no lost wake-up', async () => {
    // The first read sees an empty harness; the item arrives while that read is in flight. Waiting on
    // a wake that already happened would leave nobody to launch, or to relaunch after a lost lease.
    const s = await createSession(fixture, { host: WINDOWS, model: [{ text: 'found' }] });
    let release: () => void = () => undefined;
    const gate = new Promise<void>((r) => {
      release = r;
    });
    let first = true;
    const client = {
      ...s.client,
      execution: async () => {
        const execution = await s.client.execution();
        if (first) {
          first = false;
          await gate;
        }
        return execution;
      },
    } as typeof s.client;
    const { supervisor, runners } = supervised(s, { client });
    await settle();
    await s.send('go', { run: false });
    await settle();
    release();
    await until(s, async () => (await statuses(s))[0] === 'completed');
    expect(runners).toHaveLength(1);
    supervisor.close();
  });

  it('relaunches after a failed launch (the runner could not start)', async () => {
    const s = await createSession(fixture, { host: WINDOWS, model: [{ text: 'ok' }] });
    let attempts = 0;
    const supervisor = createInProcessSupervisor({
      client: s.client,
      clock: s.clock,
      startRunner: () => {
        attempts += 1;
        if (attempts === 1) throw new Error('no worker');
        return s.runner();
      },
    });
    const { launch } = await s.send('go', { run: false });
    const reports: LaunchReport[] = [];
    supervisor.ensure(s.harnessId, launch as string).subscribe((r) => reports.push(r));
    await until(s, async () => (await statuses(s))[0] === 'completed');
    expect(reports.map((r) => r.phase)).toEqual(['launching', 'failed']);
    expect(attempts).toBe(2);
    supervisor.close();
  });

  it('holds queued items when shouldStart says so', async () => {
    const s = await createSession(fixture, { host: WINDOWS, model: [{ text: 'ok' }] });
    const { supervisor, runners } = supervised(s, { shouldStart: () => false });
    const { launch } = await s.send('go', { run: false });
    const reports: LaunchReport[] = [];
    supervisor.ensure(s.harnessId, launch as string).subscribe((r) => reports.push(r));
    for (let i = 0; i < 20; i += 1) {
      await settle();
      await s.clock.advance(5_000);
    }
    expect(runners).toHaveLength(0);
    expect(reports.at(-1)).toMatchObject({ phase: 'failed', error: { reason: 'held' } });
    expect((await s.client.execution()).inbox).toHaveLength(1);
    supervisor.close();
  });

  it('calls onStalled once per stall, and again only after progress moved', async () => {
    let clock: Session['clock'] | undefined;
    const s = await createSession(fixture, {
      host: WINDOWS,
      tools: {
        hang: async (context) => {
          await clock?.sleep(200_000, context.signal);
          return { text: 'late' };
        },
      },
      model: [{ calls: [{ tool: 'hang', arguments: {} }] }, { text: 'final' }],
    });
    clock = s.clock;
    const stalls: number[] = [];
    const { supervisor } = supervised(s, {
      stalledAfterMs: (phase) => (phase === 'tools' ? 60_000 : undefined),
      onStalled: (state) => stalls.push(state.progress?.since ?? -1),
    });
    const { launch } = await s.send('go', { run: false });
    supervisor.ensure(s.harnessId, launch as string);
    await until(s, async () => (await statuses(s))[0] === 'completed', 2_000);
    expect(stalls).toHaveLength(1);
    supervisor.close();
  });
});
