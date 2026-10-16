/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createSession } from '../session.js';

/** O1 (Observation). */
export function observationCases(): void {
  describe('Observation', () => {
    it('O1 throwing-subscriber: a subscriber that throws changes no event', async () => {
      const script = [{ text: 'hello' }];
      const quiet = createSession({ model: script });
      await quiet.send('go');
      const loud = createSession({ model: script });
      loud.harness.observe().subscribe(() => {
        throw new Error('subscriber failure');
      });
      await loud.send('go');
      const shape = (events: Awaited<ReturnType<typeof quiet.events>>) =>
        events.map((event) => (event.kind === 'step' ? `${event.stepKind}:${event.ending}` : event.control));
      expect(shape(await loud.events())).toEqual(shape(await quiet.events()));
    });
  });
}
