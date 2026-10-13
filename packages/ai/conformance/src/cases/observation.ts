/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ToolFrame } from '@molstar/ai-harness';
import { describe, expect, it } from '@jest/globals';
import { type ConformanceFixture, createSession } from '../session.js';

/** O1 (Observation). */
export function observationCases(fixture: ConformanceFixture): void {
  describe('Observation', () => {
    it('O1 throwing-subscriber: a subscriber that throws changes no event', async () => {
      const script = [{ text: 'hello' }];
      const quiet = await createSession(fixture, { model: script });
      await quiet.send('go');
      const loud = await createSession(fixture, { model: script });
      const runner = loud.runner();
      runner.observe().subscribe(() => {
        throw new Error('subscriber failure');
      });
      await loud.send('go', { runner });
      const shape = (events: Awaited<ReturnType<typeof quiet.events>>) =>
        events.map((event) => (event.kind === 'step' ? `${event.stepKind}:${event.ending}` : event.control));
      expect(shape(await loud.events())).toEqual(shape(await quiet.events()));
    });

    it("O1 tool-progress-observed: a handler's progress frames arrive as tool.progress observations, in order (defect 17)", async () => {
      async function* frames(): AsyncGenerator<ToolFrame> {
        for (const value of [1, 2, 3]) yield { kind: 'progress', value };
        yield { kind: 'result', status: 'succeeded', parts: [{ kind: 'text', text: 'done' }] };
      }
      const s = await createSession(fixture, {
        tools: { emit: () => frames() },
        model: [{ calls: [{ tool: 'emit', arguments: {} }] }, { text: 'final' }],
      });
      const runner = s.runner();
      await s.send('go', { runner });
      const [call] = s.tools.callsOf('emit');
      const progress = runner.observations.flatMap((o) => (o.kind === 'tool.progress' ? [o] : []));
      expect(progress.map((o) => o.value)).toEqual([1, 2, 3]);
      for (const o of progress) {
        expect(o.callId).toBe(call?.callId);
        expect(o.path).toEqual([]);
      }
    });
  });
}
