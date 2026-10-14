/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories' documents of the smoke test, a share a file: smoke.test.ts and smoke-more.test.ts each open half of the
// registry's stories (every other one), so the two halves run on two workers at once. The stories with the largest
// documents are the slow tier's (`SLOW_STORIES`); every other runs in the fast tier.

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { SMOKE, type SmokeEnv } from '../src/_stories/support/seqstar-smoke.js';
import { smoke } from './smoke-harness.js';

export const root = join(dirname(fileURLToPath(import.meta.url)), '../../../..');
const env: SmokeEnv = { readFile: (path) => new Uint8Array(readFileSync(join(root, path))) };

/** The stories whose documents take most of the run (a genome, a structure, case 9's and case 17's sizes). */
export const SLOW_STORIES = new Set([
  'large-documents.stories.tsx#GenomeScale',
  'molstar.stories.tsx#P53',
  'viewer.stories.tsx#Case17',
  'viewer.stories.tsx#Case9',
]);

/** The documents of the stories at `part` of every `parts` (by their index among the stories with documents). */
export function storyDocuments(part: number, parts: number) {
  describe("the stories' documents", () => {
    const stories = SMOKE.filter((s) => s.cases);
    for (let k = part; k < stories.length; k += parts) {
      const story = stories[k];
      const make = story.cases;
      if (!make) continue;
      it(story.key, { timeout: 120_000, tags: SLOW_STORIES.has(story.key) ? ['slow'] : [] }, async () => {
        const cases = await make(env);
        expect(cases.length, 'no document').toBeGreaterThan(0);
        for (const c of cases) {
          const findings = await smoke(c);
          expect(
            findings.map((f) => `[${c.label}] ${f.step}: ${f.problem}`),
            story.key,
          ).toEqual([]);
        }
      });
    }
  });
}
