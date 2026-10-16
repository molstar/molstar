/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { hashing } from '../../messages/index.js';
import { createSession } from '../../testing/session.js';

/** Counts hashes (`hashMessage` and `digest` both go through `hashing.bytes`) while `run` works. */
async function withDigestSpy<T>(run: (digests: { calls(): number }) => Promise<T>): Promise<T> {
  const spy = vi.spyOn(hashing, 'bytes');
  try {
    return await run({ calls: () => spy.mock.calls.length });
  } finally {
    spy.mockRestore();
  }
}

/** P11 (Performance): hashing is done once. The one cost bound kept: hashing is real CPU in memory too. */
describe('Performance', () => {
  it('P11 hash-once: no round hashes more than the first, however long the head', async () => {
    await withDigestSpy(async (digests) => {
      const rounds = 16;
      const calls: number[] = [];
      const s = createSession({
        tools: {
          work: (context) => {
            calls.push(digests.calls());
            return { text: `worked ${JSON.stringify(context.arguments)}` };
          },
        },
        model: [
          ...Array.from({ length: rounds }, (_, round) => ({ calls: [{ tool: 'work', arguments: { round } }] })),
          { text: 'done' },
        ],
      });
      await s.send('go');
      const perRound = calls.slice(1).map((count, i) => count - (calls[i] as number));
      expect(digests.calls()).toBeGreaterThan(0);
      for (const count of perRound) expect(count).toBeLessThanOrEqual((perRound[0] as number) + 2);
    });
  });
});
