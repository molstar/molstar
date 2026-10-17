// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { flush } from '@molstar/prism-ui/dom/scope';
import { shownInOrder } from '@molstar/prism-ui/testing/virtual-list';
import { afterEach, describe, expect, it } from 'vitest';
import { makeScene, type SceneFixture } from '../../../_stories/performance-fixtures.js';
import { SceneTree } from '../scene-tree.js';

// The scene tree follows the scene's change stream (`SceneModel.changes`, view-models.md §3): each changed node's
// subtree is flattened again in place, the rest kept. After any run of changes, its rows (order, level, place among
// siblings, open state) are those a tree mounted on the scene as it is now makes from scratch.

(globalThis as { ResizeObserver?: unknown }).ResizeObserver ??= class {
  observe() {}
  unobserve() {}
  disconnect() {}
};

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(fx: SceneFixture) {
  const box = document.createElement('div');
  document.body.append(box);
  const dispose = render(() => <SceneTree scene={fx.scene} />, box);
  flush();
  const rows = () =>
    shownInOrder(box, '[role=treeitem]').map((e) =>
      [
        e.dataset.treeKey,
        e.getAttribute('aria-level'),
        e.getAttribute('aria-posinset'),
        e.getAttribute('aria-setsize'),
        e.getAttribute('aria-expanded'),
      ].join(' '),
    );
  return { rows, dispose };
}

/** A small random generator (the same runs every time). */
function random(seed: number) {
  let s = seed;
  return () => {
    s = (s * 1103515245 + 12345) & 0x7fffffff;
    return s / 0x7fffffff;
  };
}

describe('the scene tree, updated by its changes', () => {
  it('matches a tree made from scratch after any run of openings, closings and changed children', () => {
    // Few enough rows that all are on screen (jsdom: a screen of 800 px).
    const fx = makeScene(5, 2);
    const t = mount(fx);
    const rnd = random(7);
    const all = [fx.structure, ...fx.components];
    const removed = new Map<number, (typeof fx.components)[number]['node'][]>();
    for (let step = 0; step < 120; step++) {
      // One to three changes a task (a burst flattens once).
      const n = 1 + Math.floor(rnd() * 3);
      for (let k = 0; k < n; k++) {
        const r = rnd();
        const target = all[Math.floor(rnd() * all.length)];
        if (r < 0.6) target.expanded.set(!target.expanded.value);
        else if (target !== fx.structure) {
          // A component's representations: one removed, or all back.
          const i = fx.components.indexOf(target);
          const kids = target.children.value;
          if (kids.length > 0 && rnd() < 0.5) {
            removed.set(i, [...(removed.get(i) ?? []), kids[kids.length - 1]]);
            target.children.set(kids.slice(0, -1));
          } else {
            target.children.set([...kids, ...(removed.get(i) ?? [])]);
            removed.delete(i);
          }
        }
      }
      flush();
      const fresh = mount(fx);
      expect(t.rows(), `step ${step}`).toEqual(fresh.rows());
      fresh.dispose();
    }
    t.dispose();
  });
});
