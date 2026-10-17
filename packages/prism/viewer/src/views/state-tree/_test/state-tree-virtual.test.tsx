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
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { afterEach, describe, expect, it } from 'vitest';
import type { StateTreeModel, StateTreeNode } from '../../../model/state-tree.js';
import { StateTreeTab } from '../state-tree-tab.js';

// The state tree virtualized (view-models.md §3): 5,000 cells, only the window in the DOM, ARIA's level, place and
// open state on each row, the keyboard across cells that are not rendered.

afterEach(() => {
  document.body.innerHTML = '';
});

/** 25 structures × 4 models × 49 cells (+ the models and structures): 5,025 cells. */
function makeTree(): StateTreeNode[] {
  return Array.from({ length: 25 }, (_, s) => ({
    ref: `s${s}`,
    label: `Structure ${s + 1}`,
    icon: 'structure',
    transform: 'Structure',
    children: Array.from({ length: 4 }, (_, m) => ({
      ref: `s${s}m${m}`,
      label: `Model ${m + 1}`,
      icon: 'layers',
      transform: 'Model',
      children: Array.from({ length: 49 }, (_, c) => ({
        ref: `s${s}m${m}c${c}`,
        label: `Component ${c}`,
        icon: 'cartoon',
        transform: 'Component',
        visible: c % 7 !== 0,
      })),
    })),
  })) as StateTreeNode[];
}

function mount() {
  const roots = observableState<readonly StateTreeNode[]>(makeTree());
  const selected = observableState<string | undefined>(undefined);
  const model: StateTreeModel = {
    roots: roots.asReadonly,
    selected: selected.asReadonly,
    select: (ref) => selected.set(ref),
    inspect: () => undefined,
    toggleVisible: () => {},
    remove: () => {},
  };
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <StateTreeTab model={model} />, container);
  const tree = container.querySelector<HTMLElement>('[role=tree]')!;
  const item = (ref: string) =>
    [...tree.querySelectorAll<HTMLElement>('[role=treeitem]')].find((i) => i.dataset.treeKey === ref);
  const key = (k: string) =>
    (document.activeElement as HTMLElement).dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const focused = () => (document.activeElement as HTMLElement).dataset.treeKey;
  return { roots, selected, tree, item, key, focused, dispose };
}

describe('the state tree, virtualized', () => {
  it('renders only the rows on screen of 5,025 cells', () => {
    const t = mount();
    const shown = t.tree.querySelectorAll('[role=treeitem]').length;
    expect(shown).toBeGreaterThan(20);
    expect(shown).toBeLessThan(100);
    expect(t.item('s24m3c48')).toBeUndefined();
    t.dispose();
  });

  it('gives each row its level, place among its siblings and open state', () => {
    const t = mount();
    const s0 = t.item('s0')!;
    expect([s0.getAttribute('aria-level'), s0.getAttribute('aria-posinset'), s0.getAttribute('aria-setsize')]).toEqual([
      '1',
      '1',
      '25',
    ]);
    expect(s0.getAttribute('aria-expanded')).toBe('true');
    const c2 = t.item('s0m0c2')!;
    expect([c2.getAttribute('aria-level'), c2.getAttribute('aria-posinset'), c2.getAttribute('aria-setsize')]).toEqual([
      '3',
      '3',
      '49',
    ]);
    // A leaf has nothing to open.
    expect(c2.hasAttribute('aria-expanded')).toBe(false);
    // Nested groups are gone: rows are siblings in one list.
    expect(t.tree.querySelector('[role=group]')).toBeNull();
    t.dispose();
  });

  it('moves the focus to cells not in the DOM, closes and opens, and selects on Enter', () => {
    const t = mount();
    t.tree.focus();
    expect(t.focused()).toBe('s0');
    t.key('End');
    expect(t.focused()).toBe('s24m3c48');
    expect(t.item('s0')).toBeUndefined();
    // Up to its model (open: Left closes it), then its structure.
    t.key('ArrowLeft');
    expect(t.focused()).toBe('s24m3');
    expect(t.item('s24m3')!.getAttribute('aria-expanded')).toBe('true');
    t.key('ArrowLeft');
    expect(t.item('s24m3')!.getAttribute('aria-expanded')).toBe('false');
    t.key('ArrowRight');
    expect(t.item('s24m3')!.getAttribute('aria-expanded')).toBe('true');
    t.key('ArrowRight');
    expect(t.focused()).toBe('s24m3c0');
    t.key('ArrowLeft');
    t.key('ArrowLeft');
    t.key('ArrowLeft');
    expect(t.focused()).toBe('s24');
    t.key('ArrowLeft');
    expect(t.item('s24')!.getAttribute('aria-expanded')).toBe('false');
    t.key('Home');
    expect(t.focused()).toBe('s0');
    t.key('ArrowDown');
    t.key('ArrowDown');
    expect(t.focused()).toBe('s0m0c0');
    t.key('Enter');
    expect(t.selected.value).toBe('s0m0c0');
    expect(t.item('s0m0c0')!.getAttribute('aria-selected')).toBe('true');
    t.dispose();
  });

  it('keeps the rows of cells that did not change when the tree changes', () => {
    const t = mount();
    const before = t.item('s0m0c1')!;
    const [s0, ...rest] = t.roots.value;
    // One cell relabelled: it and its ancestors are new records (as the source makes them), its siblings are not.
    const m0 = s0.children![0];
    const c0 = { ...m0.children![0], label: 'Renamed' };
    const next = { ...s0, children: [{ ...m0, children: [c0, ...m0.children!.slice(1)] }, ...s0.children!.slice(1)] };
    t.roots.set([next, ...rest]);
    expect(t.item('s0m0c0')!.textContent).toContain('Renamed');
    expect(t.item('s0m0c1')).toBe(before);
    t.dispose();
  });

  it('keeps the rows of cells opened again when the tree is flattened afresh', () => {
    const t = mount();
    t.tree.focus();
    // s0 closed and opened (its rows made by the incremental pass), then the same cells emitted again (the full pass).
    t.key('ArrowLeft');
    t.key('ArrowRight');
    const made = t.item('s0m0')!;
    expect(made).toBeDefined();
    t.roots.set([...t.roots.value]);
    flush();
    expect(t.item('s0m0')).toBe(made);
    t.dispose();
  });
});
