// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every ring story's rendered text (the story page's summary, not the docs-mode comment) names the gestures that open a
// chromosome (the genome ring's interaction matrix, a block): the zoom, the double-click on its label or arc, the label's
// "Show only", the square; the partner keys on their story. Read through the story's own render, the element its model
// carries.

import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { describe, expect, it } from 'vitest';
import { RingViews } from '@molstar/seqx-examples/genomes/ring/views';
import { exampleStory } from '../../../examples/example-root.js';
import { KindsViews, RingStoryViews, ringStoryView } from '../views.js';

const ids = [...new Set([...RingViews, ...RingStoryViews, ...KindsViews].map((v) => v.id))];

/** The text a story renders for a ring view: its model's summary, as the page prints it. */
const rendered = (id: string, keys = false) =>
  (
    exampleStory(ringStoryView(id, keys), { stacked: true, extensions: [rondoViewer] }) as unknown as {
      model: { summary: string };
    }
  ).model.summary;

describe('the ring stories’ rendered text names the gestures that open a chromosome', () => {
  it('every ring view says the double-click, and the zoom, beside how to get back or the label’s menu', () => {
    expect(ids.length).toBeGreaterThan(8);
    for (const id of ids) {
      const text = rendered(id);
      expect(text, id).toMatch(/Cmd or Ctrl\+wheel/);
      const opens = /double-click on its label or arc/.test(text) && /Show only/.test(text) && /square/.test(text);
      expect(opens || /shows every chromosome again/.test(text), id).toBe(true);
    }
  });

  it('the genome at the fit says all of them; the partner keys story names the keys', () => {
    expect(rendered('ring')).toMatch(/double-click on its label or arc, its label’s menu “Show only”, or the square/);
    expect(rendered('ring', true)).toMatch(/`\.` and `,`/);
  });
});
