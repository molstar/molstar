/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The scene tree's hot paths in jsdom (1,000 components open: 4,000 rows; slots, `RecycledList`), shared by their quick
// cases (scene-tree.quick.ts) and heap guards (scene-tree.heap.ts): a pick moving between two rows, a scroll a row at a
// time, a jump far away. jsdom has no layout: the viewport's scroll position is a plain property here.

import { render } from '@molstar/prism-ui/dom/render';
import { flush } from '@molstar/prism-ui/dom/scope';
import { makeScene } from '../../../_stories/performance-fixtures.js';
import { SceneTree } from '../scene-tree.js';

export interface Case {
  step: (i: number) => void;
  sink: () => number;
}

function mountTree() {
  const fx = makeScene(1000, 3);
  for (const c of fx.components) c.expanded.set(true);
  const box = document.createElement('div');
  document.body.append(box);
  render(() => <SceneTree scene={fx.scene} />, box);
  flush();
  const viewport = box.querySelector<HTMLElement>('[data-axis] > div')!;
  let top = 0;
  Object.defineProperty(viewport, 'scrollTop', { get: () => top, set: (v: number) => (top = v), configurable: true });
  const scrollTo = (y: number) => {
    viewport.scrollTop = y;
    viewport.dispatchEvent(new Event('scroll'));
    flush();
  };
  return { fx, box, scrollTo, rows: () => box.querySelectorAll('[role=treeitem]').length };
}

/** The picked node moving between two shown rows. */
export function pick(): Case {
  const t = mountTree();
  const a = t.fx.components[3].node;
  const b = t.fx.components[5].node;
  return {
    step: (i) => {
      t.fx.scene.pick(i & 1 ? a : b);
      flush();
    },
    sink: () => t.rows(),
  };
}

/**
 * A row's step down, over 600 rows (then back to the top: a jump every 600 steps). Not more: jsdom's style parser keeps
 * every distinct value it parsed (a row's `translateY`), so more positions read as kept memory a browser does not keep;
 * the warm-up sees all of these.
 */
export function scrollStep(): Case {
  const t = mountTree();
  return { step: (i) => t.scrollTo((i % 600) * 26), sink: () => t.rows() };
}

/** A jump far away (a thumb drag): a window of new rows a step. */
export function jump(): Case {
  const t = mountTree();
  const span = 4000 * 26;
  return { step: (i) => t.scrollTo(Math.round((((i * 37) % 97) / 97) * span)), sink: () => t.rows() };
}
