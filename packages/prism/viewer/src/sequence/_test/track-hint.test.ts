/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { SecondaryTrack } from '../document.js';
import { createRoot } from '@molstar/prism-ui/dom/scope';
import { mountTrackGutter } from '../track-gutter.js';
import { observableState } from '@molstar/prism-ui';
import { headlessViewers, hostStructure } from './fixtures.js';

// The hidden tracks' hint (sequence.md §8): nothing while every track shows; "N hidden" at the foot of the tracks while
// some are hidden, its menu showing them again (Show all first, then one by one).

afterEach(() => {
  document.body.innerHTML = '';
});

async function until(check: () => boolean, ms = 5000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

describe('the hidden tracks hint', () => {
  it('shows only while tracks are hidden; its menu shows them again', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await until(() => h.viewer.tracks.value.length === 3);
    const root = document.createElement('div');
    root.className = 'msp-root';
    const surface = document.createElement('div');
    root.append(surface);
    document.body.append(root);
    // In a scope, as in the sequence control (its listeners are scoped: dom/scoped.ts).
    const unmount = createRoot((dispose) => {
      const off = mountTrackGutter(h.viewer, surface, root, observableState('main').asReadonly);
      return () => {
        off();
        dispose();
      };
    });
    const hint = () => [...root.querySelectorAll('button')].find((b) => b.textContent?.endsWith('hidden'));
    expect(hint()).toBeUndefined();

    await h.viewer.setTrackHidden(SecondaryTrack, true);
    await h.viewer.seq.graph.update();
    h.scheduler.advance(50);
    await until(() => !!hint());
    expect(hint()?.textContent).toBe('1 hidden');
    expect(hint()?.getAttribute('aria-label')).toBe('1 hidden track: show it');

    hint()!.click();
    const rows = [...document.querySelectorAll('[role="menuitem"]')].map((b) => b.textContent);
    expect(rows).toEqual(['Show all tracks', 'Show Secondary structure']);
    (document.querySelector('[role="menuitem"]') as HTMLElement).click();
    await until(() => h.viewer.tracks.value.every((t) => !t.hidden));
    await until(() => !hint());
    unmount();
    h.viewers.dispose();
  });
});
