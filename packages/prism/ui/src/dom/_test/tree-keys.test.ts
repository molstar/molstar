/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

import { describe, expect, it, vi } from 'vitest';
import { treeFocus, treeKeys } from '../tree-keys.js';

// The tree's keyboard (accessibility.md): arrows move, open, close, enter and leave; Enter does what a click does.

function tree() {
  const root = document.createElement('div');
  root.setAttribute('role', 'tree');
  root.tabIndex = 0;
  root.innerHTML = `
    <div role="treeitem" aria-expanded="true" tabindex="-1" id="a">
      <div data-tree-row id="a-row"><button data-tree-toggle id="a-toggle"></button></div>
      <div role="group">
        <div role="treeitem" aria-expanded="false" tabindex="-1" id="b"><div data-tree-row><button data-tree-toggle id="b-toggle"></button></div></div>
      </div>
    </div>
    <div role="treeitem" aria-expanded="false" tabindex="-1" id="c" aria-selected="true"><div data-tree-row id="c-row"></div></div>`;
  root.addEventListener('keydown', treeKeys);
  root.addEventListener('focusin', treeFocus);
  document.body.replaceChildren(root);
  const el = (id: string) => document.getElementById(id) as HTMLElement;
  const key = (k: string) =>
    (document.activeElement as HTMLElement).dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  return { root, el, key };
}

describe('the tree keyboard', () => {
  it('enters at the picked item, moves with the arrows, Home and End', () => {
    const { root, el, key } = tree();
    root.focus();
    expect(document.activeElement).toBe(el('c'));
    key('Home');
    expect(document.activeElement).toBe(el('a'));
    key('ArrowDown');
    expect(document.activeElement).toBe(el('b'));
    expect([el('a').tabIndex, el('b').tabIndex]).toEqual([-1, 0]);
    key('End');
    expect(document.activeElement).toBe(el('c'));
    key('ArrowUp');
    expect(document.activeElement).toBe(el('b'));
  });

  it('opens and closes with Right and Left, enters and leaves, and clicks the row on Enter', () => {
    const { root, el, key } = tree();
    root.focus();
    key('Home');
    key('ArrowRight');
    expect(document.activeElement).toBe(el('b'));
    const open = vi.fn();
    el('b-toggle').addEventListener('click', open);
    key('ArrowRight');
    expect(open).toHaveBeenCalledOnce();
    key('ArrowLeft');
    expect(document.activeElement).toBe(el('a'));
    const close = vi.fn();
    el('a-toggle').addEventListener('click', close);
    key('ArrowLeft');
    expect(close).toHaveBeenCalledOnce();
    const pick = vi.fn();
    el('a-row').addEventListener('click', pick);
    key('Enter');
    expect(pick).toHaveBeenCalledOnce();
  });
});
