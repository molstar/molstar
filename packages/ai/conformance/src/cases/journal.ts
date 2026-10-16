/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { hashMessage, type Message } from '@molstar/ai-agent';
import { type ConformanceFixture, createSession } from '../session.js';
import { settle } from '../testing/clock.js';
import { toolDefinition } from '../testing/messages.js';

/** Polls (letting the runner and the simulated clock's promises settle) until the condition holds. */
async function until(condition: () => boolean, rounds = 200): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (condition()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

/** Every way to change one part of a value: a leaf, a dropped or added element, a reordering. */
function variants(value: unknown): unknown[] {
  if (typeof value === 'string') return [`${value}!`];
  if (typeof value === 'number') return [value + 1];
  if (typeof value === 'boolean') return [!value];
  if (value === null) return ['not null'];
  if (Array.isArray(value)) {
    const reversed = [...value].reverse();
    return [
      ...value.flatMap((item, i) => variants(item).map((changed) => value.with(i, changed))),
      value.slice(0, -1),
      [...value, 'extra'],
      ...(JSON.stringify(reversed) === JSON.stringify(value) ? [] : [reversed]),
    ];
  }
  const record = value as Record<string, unknown>;
  return [
    ...Object.entries(record).flatMap(([key, item]) =>
      variants(item).map((changed) => ({ ...record, [key]: changed })),
    ),
    ...Object.keys(record).map((key) => Object.fromEntries(Object.entries(record).filter(([k]) => k !== key))),
    { ...record, extra: 1 },
  ];
}

/** The same value with every object's keys in reverse order. */
const reorder = (value: unknown): unknown =>
  Array.isArray(value)
    ? value.map(reorder)
    : value !== null && typeof value === 'object'
      ? Object.fromEntries(
          Object.entries(value)
            .reverse()
            .map(([key, item]) => [key, reorder(item)]),
        )
      : value;

/** J1–J4 (Journal). */
export function journalCases(fixture: ConformanceFixture): void {
  describe('Journal', () => {
    it('J3 stop-during-tools: Stop while tools run leaves main as it was before the round', async () => {
      let clock: { sleep(ms: number, signal?: AbortSignal): Promise<void> } | undefined;
      const s = await createSession(fixture, {
        tools: {
          slow: async (context) => {
            await clock?.sleep(60_000, context.signal);
            return { text: 'done' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.tools.callsOf('slow').length === 1);
      const before = await s.main();
      expect(before.some((message) => message.role === 'assistant')).toBe(false);

      await s.client.stop({ key: 'stop-1', turnId: (await s.client.execution()).turn?.id as string });
      await s.clock.runUntil(running);

      expect(await s.main()).toEqual(before);
      expect((await s.turns())[0]?.status).toBe('stopped');
    });

    it('J4 hash-stable: hashMessage ignores key order and changes with any part', async () => {
      const messages: Message[] = [
        {
          role: 'system',
          parts: [
            { kind: 'text', text: 'be brief' },
            {
              kind: 'tool-definition',
              definition: toolDefinition('lookup', { q: { type: 'string' } }),
              deferred: false,
            },
          ],
        },
        {
          role: 'assistant',
          parts: [
            { kind: 'text', text: 'looking' },
            {
              kind: 'tool-call',
              requestId: 'req-1',
              tool: { id: 'lookup' },
              arguments: { q: 'x', options: { a: 1, b: [true, null, 's'] } },
              providerOptions: { signature: 'abc' },
            },
          ],
        },
        {
          role: 'user',
          parts: [
            {
              kind: 'tool-result',
              requestId: 'req-1',
              tool: { id: 'lookup' },
              status: 'succeeded',
              parts: [{ kind: 'text', text: 'found' }],
            },
            { kind: 'reference', reference: { uri: 'mem://a', mediaType: 'text/plain', byteLength: 3 } },
          ],
        },
      ];

      for (const message of messages) {
        const hash = await hashMessage(message);
        expect(await hashMessage(structuredClone(message))).toBe(hash);
        expect(await hashMessage(reorder(message) as Message)).toBe(hash);

        const changed = variants(message);
        expect(changed.length).toBeGreaterThan(10);
        for (const other of changed) expect(await hashMessage(other as Message)).not.toBe(hash);
      }
    });
  });
}
