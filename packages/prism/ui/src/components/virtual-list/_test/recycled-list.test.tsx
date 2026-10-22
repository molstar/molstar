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
import { render } from '../../../dom/render.js';
import { flush, flushFrame, watch } from '../../../dom/scope.js';
import { TreeNavigation } from '../../../dom/tree-keys.js';
import { setAttr, setText } from '../../../dom/writes.js';
import { observableState } from '../../../reactive/observable-state.js';
import { shownInOrder } from '../../../testing/virtual-list.js';
import { RecycledList } from '../recycled-list.js';
import type { VirtualListHandle } from '../virtual-list.js';
import { MaxScreensAhead, rowsAhead } from '../window.js';

// Recycled rows (view-models.md §3): a pool of slots, each made once, pointed at the items a moving window brings in;
// a window that moves creates, removes and reorders no node; the keyboard and the focus follow the rows, not the
// slots; keyed rows (what `recycle` does not take) live beside them.

afterEach(() => {
  document.body.innerHTML = '';
});

interface Item {
  readonly id: string;
  readonly label: string;
  readonly special?: boolean;
}

function mount(options: { count?: number; keyed?: boolean; reveal?: boolean } = {}) {
  const items = observableState<readonly Item[]>(
    Array.from({ length: options.count ?? 5000 }, (_, i) => ({
      id: `r${i}`,
      label: `Row ${i}`,
      special: options.keyed && i % 10 === 5,
    })),
  );
  const selected = observableState<string | undefined>(undefined);
  let made = 0;
  let keyedMade = 0;
  let list: VirtualListHandle | undefined;
  const nav = new TreeNavigation<Item>({
    rows: () => items.value,
    key: (r) => r.id,
    level: () => 1,
    scrollTo: (i) => list?.scrollToIndex(i),
  });
  const slot = () => {
    made++;
    let item: Item | undefined;
    const text = document.createTextNode('');
    const el = (
      <div role="treeitem" tabIndex={-1}>
        <div data-tree-row="" onClick={() => item && selected.set(item.id)}>
          {text}
        </div>
      </div>
    ) as HTMLElement;
    const point = nav.slot(el);
    const draw = () => {
      if (!item) return;
      setText(text, item.label);
      setAttr(el, 'aria-selected', String(selected.value === item.id));
    };
    watch([selected.asReadonly], draw);
    return {
      el,
      show(next: Item | undefined) {
        item = next;
        point(next?.id);
        draw();
      },
    };
  };
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <RecycledList
        items={items.asReadonly}
        itemKey={(r) => r.id}
        rowHeight={26}
        role="tree"
        aria-label="List"
        contentRef={nav.tree}
        ref={(l) => {
          list = l;
        }}
        slot={slot}
        recycle={(r) => !r.special}
        reveal={options.reveal}
      >
        {(r) => {
          keyedMade++;
          return (
            <div role="treeitem" tabIndex={-1} class="special" ref={nav.row(r.id)}>
              {r.label}
            </div>
          );
        }}
      </RecycledList>
    ),
    container,
  );
  const tree = container.querySelector<HTMLElement>('[role=tree]')!;
  const shown = () => shownInOrder(tree, '[role=treeitem]').map((e) => e.textContent);
  const key = (k: string) =>
    (document.activeElement as HTMLElement).dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const focused = () => (document.activeElement as HTMLElement).dataset.treeKey;
  return {
    items,
    selected,
    tree,
    shown,
    key,
    focused,
    list: () => list!,
    made: () => made,
    keyedMade: () => keyedMade,
    dispose,
  };
}

describe('the margin ahead of the scroll (window.ts)', () => {
  it('grows with the speed, two frames of it, up to three screens', () => {
    // 31 rows a screen of 26 px rows.
    expect(rowsAhead(0, 26, 31, 16)).toBe(16);
    expect(rowsAhead(260, 26, 31, 16)).toBe(20);
    expect(rowsAhead(1200, 26, 31, 16)).toBe(MaxScreensAhead * 31);
  });

  it("a wheel's large steps land on rows made from the second step on", () => {
    const t = mount();
    const viewport = t.tree.closest('[data-axis]')!.querySelector<HTMLElement>(':scope > div')!;
    let top = 0;
    Object.defineProperty(viewport, 'scrollTop', { get: () => top, set: (v: number) => (top = v), configurable: true });
    // jsdom lays nothing out: the list's screen is its default 800 px (31 rows).
    const madeRows = () => new Set(t.shown().map((l) => Number(l!.slice(4))));
    const step = (y: number) => {
      // A wheel's: the compositor's scroll (a scroll alone, a thumb's drag, has no margin for speed).
      viewport.dispatchEvent(new Event('wheel'));
      top = y;
      viewport.dispatchEvent(new Event('scroll'));
      flushFrame();
      flush();
    };
    // 1,200 px a frame (46 rows): the first step jumps; from the second, the next step's screen is made already.
    step(1200);
    step(2400);
    const made = madeRows();
    const next = Math.floor(3600 / 26);
    for (let i = next; i < next + 31; i++) expect(made.has(i), `row ${i}`).toBe(true);
    t.dispose();
  });
});

describe('RecycledList', () => {
  it('points its slots at the rows a moved window shows, and makes, removes and moves no node', () => {
    const t = mount();
    const first = t.shown();
    expect(first[0]).toBe('Row 0');
    expect(first.length).toBeGreaterThan(20);
    const pool = t.made();
    const changes: MutationRecord[] = [];
    const observer = new MutationObserver((records) => changes.push(...records.filter((r) => r.type === 'childList')));
    observer.observe(t.tree, { childList: true, subtree: true });
    // Far away, then a step, then back.
    t.list().scrollToIndex(3000);
    t.list().scrollToIndex(3040);
    t.list().scrollToIndex(10);
    t.list().scrollToIndex(4999);
    observer.takeRecords();
    observer.disconnect();
    expect(t.made()).toBe(pool);
    expect(changes).toEqual([]);
    const last = t.shown();
    expect(last[last.length - 1]).toBe('Row 4999');
    // Each shown row is the row at its place: consecutive labels.
    const n = last.map((l) => Number(l!.slice(4)));
    for (let i = 1; i < n.length; i++) expect(n[i]).toBe(n[i - 1] + 1);
    t.dispose();
  });

  it("follows the rows' changes and the list's: a new record, a selection, the rows replaced", () => {
    const t = mount({ count: 100 });
    t.selected.set('r2');
    // `watch` renders in a microtask (rendering.md §2.3): read after a flush.
    flush();
    const row2 = () => t.tree.querySelector<HTMLElement>('[data-tree-key="r2"]')!;
    expect(row2().getAttribute('aria-selected')).toBe('true');
    const items = t.items.value.slice();
    items[2] = { ...items[2], label: 'Renamed' };
    const before = row2();
    t.items.set(items);
    expect(row2()).toBe(before);
    expect(row2().textContent).toBe('Renamed');
    // A row inserted at the top: every row moves down one, no node is made.
    const pool = t.made();
    t.items.set([{ id: 'new', label: 'New' }, ...items]);
    expect(t.shown().slice(0, 4)).toEqual(['New', 'Row 0', 'Row 1', 'Renamed']);
    expect(t.made()).toBe(pool);
    // Fewer rows than slots: the rest are put away (hidden, showing nothing).
    t.items.set(items.slice(0, 3));
    expect(t.shown()).toEqual(['Row 0', 'Row 1', 'Renamed']);
    expect(t.tree.querySelectorAll('[data-tree-key]')).toHaveLength(3);
    t.dispose();
  });

  it('keeps what it does not recycle as keyed rows beside the slots', () => {
    const t = mount({ keyed: true });
    expect(t.shown().slice(0, 7)).toEqual(['Row 0', 'Row 1', 'Row 2', 'Row 3', 'Row 4', 'Row 5', 'Row 6']);
    expect(t.tree.querySelector('.special')!.textContent).toBe('Row 5');
    const made = t.keyedMade();
    t.list().scrollToIndex(2000);
    expect(t.keyedMade()).toBeGreaterThan(made);
    expect(t.shown()).toContain('Row 1995');
    expect([...t.tree.querySelectorAll('.special')].map((e) => e.textContent)).toContain('Row 1995');
    t.dispose();
  });

  it('moves the focus across re-pointed rows, and to the tree when its row is scrolled away', async () => {
    const t = mount();
    t.tree.focus();
    expect(t.focused()).toBe('r0');
    t.key('End');
    expect(t.focused()).toBe('r4999');
    expect(t.tree.querySelector('[data-tree-key="r0"]')).toBeNull();
    t.key('ArrowUp');
    t.key('ArrowUp');
    expect(t.focused()).toBe('r4997');
    t.key('Home');
    expect(t.focused()).toBe('r0');
    t.key('Enter');
    expect(t.selected.value).toBe('r0');
    // Its slot pointed at another row: the focus does not stay with the slot, it goes to the tree.
    const slot = document.activeElement as HTMLElement;
    t.list().scrollToIndex(3000);
    expect(slot.dataset.treeKey).not.toBe('r0');
    await Promise.resolve();
    expect(document.activeElement).toBe(t.tree);
    // And from there, the keys go on from the row it had (scrolled back into view).
    t.key('ArrowDown');
    expect(t.focused()).toBe('r1');
    t.dispose();
  });

  it('at rest, the rows are in the DOM in screen order; a focused row keeps its focus', async () => {
    vi.useFakeTimers();
    try {
      const t = mount({ count: 2000 });
      const domOrder = () =>
        [...t.tree.querySelectorAll<HTMLElement>('[role=treeitem]')]
          .filter((el) => !el.closest('[hidden]'))
          .map((el) => el.textContent);
      // Jumps away and back re-point slots out of order.
      for (const i of [900, 40, 1500, 300]) t.list().scrollToIndex(i);
      expect(domOrder()).not.toEqual(t.shown());
      // A shown row, not the first in the DOM (the pass moves rows before it).
      const rows = [...t.tree.querySelectorAll<HTMLElement>('[role=treeitem]')].filter((el) => !el.closest('[hidden]'));
      const row = rows[rows.length - 1];
      row.focus();
      expect(document.activeElement).toBe(row);
      vi.advanceTimersByTime(200);
      expect(domOrder()).toEqual(t.shown());
      expect(document.activeElement).toBe(row);
      t.dispose();
    } finally {
      vi.useRealTimers();
    }
  });
});

describe('a reveal (a node opened)', () => {
  it('shows no skeleton behind the rows fading in, and the skeleton again once they are in', async () => {
    vi.useFakeTimers();
    try {
      const m = mount({ count: 10, reveal: true });
      await flushFrame();
      const revealing = () => !!document.querySelector('[data-revealing]');
      expect(revealing()).toBe(false);
      const now = m.items.value;
      const opened = [0, 1, 2].map((i) => ({ id: `c${i}`, label: `Child ${i}` }));
      m.items.set([...now.slice(0, 3), ...opened, ...now.slice(3)]);
      flush();
      await flushFrame();
      expect(revealing()).toBe(true);
      vi.advanceTimersByTime(1000);
      expect(revealing()).toBe(false);
      m.dispose();
    } finally {
      vi.useRealTimers();
    }
  });
});
