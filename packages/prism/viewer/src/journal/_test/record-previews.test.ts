/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ThreadHost } from '@molstar/ai';
import { describe, expect, it, vi } from 'vitest';
import { createSessionRecord, previewsHeld, recordLines, type TurnEntry } from '../record.js';
import { fakeApi } from './runs.js';

// A turn's streamed reply read by several readers at once (journal.md §2): History and the assistant both read while a
// turn streams, and an export can race a refresh.

const hosts: ThreadHost[] = [];
vi.mock('@molstar/ai', async (original) => {
  const ai = await original<typeof import('@molstar/ai')>();
  return {
    ...ai,
    createThreadHost: (...args: Parameters<typeof ai.createThreadHost>) => {
      const host = ai.createThreadHost(...args);
      hosts.push(host);
      return host;
    },
  };
});

function manual() {
  const timers = new Set<() => void>();
  return {
    after: (fn: () => void) => {
      timers.add(fn);
      return () => timers.delete(fn);
    },
    fire: () => {
      for (const fn of [...timers]) {
        timers.delete(fn);
        fn();
      }
    },
  };
}

const replyOf = (page: { items: { entry: unknown }[] }) => (page.items[0].entry as TurnEntry).reply;

describe("a streaming reply's preview", () => {
  it('reads the same text however many readers read at once, and after', async () => {
    const { api } = fakeApi();
    const m = manual();
    const record = await createSessionRecord({ api, after: m.after });
    const turn = record.beginTurn('t1', 'hi');
    turn.reply('Hel');
    m.fire();
    await record.flush();
    turn.reply('Hello');
    m.fire();
    await record.flush();
    const [a, b, c] = await Promise.all([record.read(10), record.read(10), record.readAfter(0, 10)]);
    expect([replyOf(a), replyOf(b), replyOf(c)]).toEqual(['Hello', 'Hello', 'Hello']);
    expect(replyOf(await record.read(10))).toBe('Hello');

    // More arrives; an export racing a refresh still folds each delta once.
    turn.reply('Hello there');
    m.fire();
    await record.flush();
    const [lines, page] = await Promise.all([recordLines(record), record.read(10)]);
    expect(replyOf(page)).toBe('Hello there');
    expect(JSON.parse(lines.split('\n').at(-1)!).value.reply).toBe('Hello there');
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    await record.dispose();
  });

  it('is let go once its turn settles', async () => {
    hosts.length = 0;
    const { api } = fakeApi();
    const m = manual();
    const record = await createSessionRecord({ api, after: m.after });
    const host = hosts[0];
    const turn = record.beginTurn('t1', 'hi');
    turn.reply('Hello');
    m.fire();
    await record.flush();
    expect(replyOf(await record.read(10))).toBe('Hello');
    expect(previewsHeld(host)).toBe(1);
    await turn.settle({ reply: 'Hello, done.', outcome: 'done' });
    expect(previewsHeld(host)).toBe(0);
    expect(replyOf(await record.read(10))).toBe('Hello, done.');
    expect(previewsHeld(host)).toBe(0);
    await record.dispose();
  });
});
