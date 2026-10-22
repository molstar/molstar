/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  type FlatNode,
  flattenTree,
  keepRows,
  VirtualList,
  type VirtualListHandle,
} from '../../components/virtual-list/virtual-list.js';
import { observableState } from '../../reactive/observable-state.js';
import { render } from '../render.js';
import { TreeNavigation } from '../tree-keys.js';

import { derive } from '../../reactive/state.js';
// The tree's keyboard (accessibility.md) over a virtualized tree: the arrows, Home and End move across rows that are
// not in the DOM (scrolled into view, then focused), Right and Left open, close, enter and leave, Enter does what a
// click does; focus roves, and the tree keeps it when its row scrolls away.

afterEach(() => {
  document.body.innerHTML = '';
});

interface N {
  id: string;
  children?: N[];
}

/** 50 groups × 99 leaves (+ the groups): 5,000 rows, all open. */
function makeTree(): N[] {
  return Array.from({ length: 50 }, (_, g) => ({
    id: `g${g}`,
    children: Array.from({ length: 99 }, (_, c) => ({ id: `g${g}c${c}` })),
  }));
}

function mount(selected?: string) {
  const roots = makeTree();
  const collapsed = observableState<ReadonlySet<string>>(new Set());
  const rows = observableState<readonly FlatNode<N>[]>([]);
  // As the views do: a row that stays the same node at the same depth is the same item (its nodes kept).
  const keep = keepRows<FlatNode<N>>(
    (r) => r.node.id,
    (a, b) => a.node === b.node && a.depth === b.depth,
  );
  const reflow = () =>
    rows.set(
      keep(
        flattenTree(
          roots,
          (n) => n.children,
          (n) => !collapsed.value.has(n.id),
        ),
      ),
    );
  reflow();
  const toggle = (id: string) => {
    const next = new Set(collapsed.value);
    if (next.has(id)) next.delete(id);
    else next.add(id);
    collapsed.set(next);
    reflow();
  };
  const picks: string[] = [];
  let list: VirtualListHandle | undefined;
  const lookups = { count: 0 };
  const nav = new TreeNavigation<FlatNode<N>>({
    rows: () => rows.value,
    key: (r) => {
      lookups.count++;
      return r.node.id;
    },
    level: (r) => r.depth + 1,
    scrollTo: (i) => list?.scrollToIndex(i),
    selected: () => selected,
  });
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <VirtualList
        items={rows.asReadonly}
        itemKey={(r) => r.node.id}
        rowHeight={26}
        role="tree"
        aria-label="Tree"
        contentRef={nav.tree}
        ref={(l) => {
          list = l;
        }}
      >
        {(r) => (
          <div
            role="treeitem"
            id={r.node.id}
            aria-level={r.depth + 1}
            aria-expanded={r.node.children ? derive(collapsed.asReadonly, (c) => !c.has(r.node.id)) : undefined}
            tabIndex={-1}
            ref={nav.row(r.node.id)}
          >
            <div data-tree-row="" onClick={() => picks.push(r.node.id)}>
              <button
                type="button"
                data-tree-toggle=""
                onClick={(e) => {
                  e.stopPropagation();
                  toggle(r.node.id);
                }}
              />
            </div>
          </div>
        )}
      </VirtualList>
    ),
    container,
  );
  const tree = container.querySelector<HTMLElement>('[role=tree]')!;
  const el = (id: string) => document.getElementById(id);
  const key = (k: string) =>
    (document.activeElement as HTMLElement).dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const focused = () => (document.activeElement as HTMLElement).id;
  return { tree, el, key, focused, picks, rows, reflow, roots, lookups, list: () => list!, dispose };
}

describe('VirtualList', () => {
  it('renders only the window of a 5,000-row tree, and a row scrolled to', () => {
    const t = mount();
    expect(t.rows.value).toHaveLength(5000);
    const shown = t.tree.querySelectorAll('[role=treeitem]').length;
    // A screenful (800px before the first measure) and the overscan.
    expect(shown).toBeGreaterThan(20);
    expect(shown).toBeLessThan(100);
    expect(t.el('g49c98')).toBeNull();
    t.dispose();
  });
});

describe('the tree keyboard', () => {
  it('enters at the picked row (rendered for it), moves with the arrows, Home and End across rows not in the DOM', () => {
    const t = mount('g30c5');
    expect(t.el('g30c5')).toBeNull();
    t.tree.focus();
    expect(t.focused()).toBe('g30c5');
    // One tab stop: the active row, not the tree.
    expect(t.el('g30c5')!.tabIndex).toBe(0);
    expect(t.tree.tabIndex).toBe(-1);
    t.key('Home');
    expect(t.focused()).toBe('g0');
    expect(t.el('g30c5')).toBeNull();
    t.key('ArrowDown');
    expect(t.focused()).toBe('g0c0');
    expect([t.el('g0')!.tabIndex, t.el('g0c0')!.tabIndex]).toEqual([-1, 0]);
    t.key('End');
    expect(t.focused()).toBe('g49c98');
    t.key('ArrowUp');
    expect(t.focused()).toBe('g49c97');
    expect(t.tree.querySelectorAll('[role=treeitem]').length).toBeLessThan(100);
    t.dispose();
  });

  it('finds the active row where it was, not by a scan of every row (perf audit 11/P9)', () => {
    const t = mount();
    t.tree.focus();
    t.key('End');
    expect(t.focused()).toBe('g49c98');
    t.lookups.count = 0;
    for (let i = 0; i < 10; i++) t.key('ArrowUp');
    expect(t.focused()).toBe('g49c88');
    // About five lookups a key (the row and its neighbor, on the key and on the focus), not a scan of 5,000 rows.
    expect(t.lookups.count).toBeLessThanOrEqual(10 * 6);
    t.dispose();
  });

  it('closes and opens with Left and Right, goes to the parent and the first child, clicks the row on Enter', () => {
    const t = mount();
    t.tree.focus();
    expect(t.focused()).toBe('g0');
    t.key('End');
    // A leaf: Left goes to its parent, 99 rows up.
    t.key('ArrowLeft');
    expect(t.focused()).toBe('g49');
    t.key('ArrowLeft');
    expect(t.el('g49')!.getAttribute('aria-expanded')).toBe('false');
    expect(t.rows.value).toHaveLength(5000 - 99);
    t.key('ArrowRight');
    expect(t.el('g49')!.getAttribute('aria-expanded')).toBe('true');
    t.key('ArrowRight');
    expect(t.focused()).toBe('g49c0');
    t.key('Enter');
    t.key(' ');
    expect(t.picks).toEqual(['g49c0', 'g49c0']);
    t.dispose();
  });

  it('keeps the focus in the tree when its row scrolls away, and moves on from it', () => {
    const t = mount();
    t.tree.focus();
    t.key('ArrowDown');
    expect(t.focused()).toBe('g0c0');
    // Scrolled far away (as the wheel would): the row leaves the DOM, the tree takes the focus and is the tab stop.
    t.list().scrollToIndex(3000);
    expect(t.el('g0c0')).toBeNull();
    expect(document.activeElement).toBe(t.tree);
    expect(t.tree.tabIndex).toBe(0);
    t.key('ArrowDown');
    expect(t.focused()).toBe('g0c1');
    t.dispose();
  });

  it('gives the focus back to a row rendered anew, and to the row in the place of one that went', async () => {
    const t = mount();
    t.tree.focus();
    t.key('ArrowDown');
    t.key('ArrowDown');
    expect(t.focused()).toBe('g0c1');
    // Its record replaced: a new element, focused again once the list has settled.
    t.roots[0].children![1] = { id: 'g0c1' };
    t.reflow();
    await Promise.resolve();
    expect(t.focused()).toBe('g0c1');
    // Gone: the row now in its place.
    t.roots[0].children!.splice(1, 1);
    t.reflow();
    await Promise.resolve();
    expect(t.focused()).toBe('g0c2');
    t.dispose();
  });

  it('ignores keys on a button in a row', () => {
    const t = mount();
    t.tree.focus();
    const toggle = t.el('g0')!.querySelector('button')!;
    const click = vi.fn();
    toggle.addEventListener('click', click);
    toggle.focus();
    t.key('ArrowDown');
    t.key('Enter');
    expect(document.activeElement).toBe(toggle);
    expect(click).not.toHaveBeenCalled();
    t.dispose();
  });
});
