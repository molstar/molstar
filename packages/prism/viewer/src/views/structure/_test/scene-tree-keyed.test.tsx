/** @vitest-environment jsdom */
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
import { ReactiveCounts } from '@molstar/prism-ui/reactive/state';
import { afterEach, describe, expect, it } from 'vitest';
import { makeScene } from '../../../_stories/performance-fixtures.js';
import { SceneTree } from '../scene-tree.js';

// Rows do not react to what did not change them (view-models.md §3): a row follows its own keys of the tree's shared
// states (the picked node by `keyedSelector`, its place by `keyedState`), so a pick renders the two rows it concerns
// and a reflow below the window renders none of the mounted rows.

afterEach(() => {
  document.body.innerHTML = '';
});

function mount() {
  const fx = makeScene(300, 3);
  const box = document.createElement('div');
  document.body.append(box);
  const dispose = render(() => <SceneTree scene={fx.scene} />, box);
  flush();
  const rows = () => box.querySelectorAll('[role=treeitem]').length;
  const renders = (change: () => void) => {
    ReactiveCounts.render = 0;
    change();
    flush();
    return ReactiveCounts.render;
  };
  return { fx, box, dispose, rows, renders };
}

describe('the scene tree, keyed', () => {
  it('a pick renders the two rows it moves between', () => {
    const t = mount();
    expect(t.rows()).toBeGreaterThan(20);
    expect(t.renders(() => t.fx.scene.pick(t.fx.components[3].node))).toBe(1);
    expect(t.renders(() => t.fx.scene.pick(t.fx.components[5].node))).toBe(2);
    expect(t.box.querySelector('[data-node-id="c5"]')!.closest('[role=treeitem]')!.getAttribute('aria-selected')).toBe(
      'true',
    );
    expect(t.box.querySelector('[data-node-id="c3"]')!.closest('[role=treeitem]')!.getAttribute('aria-selected')).toBe(
      'false',
    );
    expect(t.renders(() => t.fx.scene.pick(undefined))).toBe(1);
    t.dispose();
  });

  it('a reflow below the window renders none of the mounted rows; one above them renders those that moved', () => {
    const t = mount();
    const mounted = t.rows();
    expect(t.renders(() => t.fx.components[250].expanded.set(true))).toBe(0);
    expect(t.renders(() => t.fx.components[250].expanded.set(false))).toBe(0);
    expect(t.rows()).toBe(mounted);
    // A sibling opened on screen: its own row (open now) and the rows after it in the window (moved down), not those
    // before it.
    const n = t.renders(() => t.fx.components[2].expanded.set(true));
    expect(n).toBeGreaterThan(0);
    expect(n).toBeLessThan(mounted);
    t.dispose();
  });
});
