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
