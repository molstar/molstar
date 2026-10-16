/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { callOf, fixture } from './helpers.js';

describe('a session with nothing parked', () => {
  it('resume() answers undefined at once and leaves the session free', async () => {
    const { session } = fixture([{ text: 'Done.' }]);
    await expect(session.resume()).resolves.toBeUndefined();
    expect(session.busy).toBe(false);
    expect((await session.turn({ prompt: 'go' })).status).toBe('completed');
    await expect(session.resume()).resolves.toBeUndefined();
  });

  it('resume() during a turn answers undefined and leaves the turn alone', async () => {
    const { session } = fixture([{ calls: [callOf('get_state')] }, { text: 'Done.' }]);
    const running = session.turn({ prompt: 'go' });
    expect(session.busy).toBe(true);
    await expect(session.resume()).resolves.toBeUndefined();
    await expect(session.forget()).rejects.toThrow('stop it first');
    expect((await running).status).toBe('completed');
    expect(session.busy).toBe(false);
  });

  it('a stop while a turn is composed: the turn does not begin', async () => {
    const { session, model } = fixture([{ text: 'Done.' }]);
    const running = session.turn({ prompt: 'go' });
    await session.stop();
    expect((await running).status).toBe('stopped');
    expect(model.invocations).toHaveLength(0);
  });
});

describe('a resumed turn (state after a resume)', () => {
  it('reports the parked step’s attempts as tools and shows the calls that never ran as skipped', async () => {
    const f = fixture(
      [{ calls: [callOf('ask_service'), callOf('fetch_write'), callOf('set_value')] }, { text: 'Done.' }],
      {
        toolRetries: 0,
      },
    );
    f.api.flaky.failures = 1;
    const parked = await f.session.turn({ prompt: 'go' });
    expect(parked.status).toBe('awaiting-approval');
    const statuses: string[] = [];
    const call = parked.approval?.[0]?.callId as string;
    const r = await f.session.resume({
      decisions: { [call]: 'allow' },
      hooks: { onStatus: (e) => statuses.push(e.kind) },
    });
    expect(r?.status).toBe('completed');
    // ask_service and fetch_write ran as tools, then the model's reply.
    expect(statuses.filter((s) => s === 'tools' || s === 'model')).toEqual(['tools', 'tools', 'model']);
    expect(r?.tools.map((t) => `${t.name}:${t.status}`)).toEqual([
      'ask_service:succeeded',
      'fetch_write:failed',
      'set_value:skipped',
    ]);
    expect(r?.text).toBe('Done.');
  });

  it('shows a repair made before the park once, and keeps it in the result', async () => {
    const f = fixture([{ calls: [callOf('ghost'), callOf('ask_service')] }, { text: 'Done.' }]);
    const parked = await f.session.turn({ prompt: 'go' });
    expect(parked.tools.map((t) => `${t.name}:${t.status}`)).toEqual(['ghost:failed']);
    const shown: string[] = [];
    const call = parked.approval?.[0]?.callId as string;
    const r = await f.session.resume({
      decisions: { [call]: 'allow' },
      hooks: { onTool: (e) => shown.push(`${e.name}:${e.status}`) },
    });
    expect(shown).toEqual(['ask_service:started', 'ask_service:succeeded']);
    expect(r?.tools.map((t) => `${t.name}:${t.status}`)).toEqual(['ghost:failed', 'ask_service:succeeded']);
  });

  it('shows a denial once when the turn parks again for another call', async () => {
    const f = fixture([{ calls: [callOf('ask_service'), callOf('ask_service', { again: 1 })] }, { text: 'Done.' }]);
    const parked = await f.session.turn({ prompt: 'go' });
    const [first, second] = parked.approval ?? [];
    const again = await f.session.resume({ decisions: { [first?.callId as string]: 'deny' } });
    expect(again?.status).toBe('awaiting-approval');
    const shown: string[] = [];
    const r = await f.session.resume({
      decisions: { [second?.callId as string]: 'allow' },
      hooks: { onTool: (e) => e.status !== 'started' && shown.push(`${e.name}:${e.status}`) },
    });
    expect(shown).toEqual(['ask_service:succeeded']);
    expect(r?.tools.map((t) => `${t.name}:${t.status}`)).toEqual(['ask_service:failed', 'ask_service:succeeded']);
  });
});
