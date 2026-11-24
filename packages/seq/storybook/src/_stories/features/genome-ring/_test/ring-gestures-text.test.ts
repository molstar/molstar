/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every ring story's text names the gestures that open a chromosome (the genome ring's interaction matrix, a block): the
// zoom, the double-click on its label or arc, the label's "Show only", the square; the partner keys on their story.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';

const source = readFileSync(new URL('../genome-ring.stories.ts', import.meta.url), 'utf8');

/** Each story's doc comment, by its export's name. */
function docs(): Map<string, string> {
  const out = new Map<string, string>();
  for (const m of source.matchAll(/\/\*\*((?:(?!\*\/)[\s\S])*)\*\/\s*export const (\w+): StoryObj/g))
    out.set(m[2], m[1].replace(/\n\s*\*\s?/g, ' '));
  return out;
}

describe('the ring stories name the gestures that open a chromosome', () => {
  const all = docs();

  it('every story of the ring says them', () => {
    expect(all.size).toBeGreaterThan(20);
    for (const [name, text] of all) {
      const opens = /double-click on its label or arc/.test(text) && /Show only/.test(text) && /square/.test(text);
      const back = /wheel out/.test(text);
      expect(opens || back, name).toBe(true);
      expect(/Cmd or Ctrl\+wheel/.test(text), name).toBe(true);
    }
  });

  it('the partner keys story names the keys', () => {
    expect(all.get('ThePartnerKeys')).toMatch(/`\.` and `,`/);
  });
});
