/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush, flushFrame, render, ViewScope } from '../../../dom/index.js';
import { state } from '../../../state.js';
import type { VirtualListHandle } from '../list-window.js';
import { RecycledList, RecycledListStats } from '../recycled-list.js';
import { VirtualList } from '../virtual-list.js';

// What both lists share (list-window.ts): rows that say their place in the whole list, a rest that waits while the list
// is paused (dom §6), and (RecycledList) a row height measured from a row that is shown.

afterEach(() => {
  document.body.innerHTML = '';
  vi.useRealTimers();
});

const items = (n: number) => Array.from({ length: n }, (_, i) => ({ id: `r${i}`, label: `Row ${i}` }));
type Item = { id: string; label: string };

function mountVirtual(n: number, parent?: ViewScope) {
  const list = state<readonly Item[]>(items(n));
  let handle: VirtualListHandle | undefined;
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <VirtualList items={list.asReadonly} itemKey={(r) => r.id} role="listbox" positions ref={(h) => (handle = h)}>
        {(r) => <div role="option">{r.label}</div>}
      </VirtualList>
    ),
    container,
    parent,
  );
  return { list, container, dispose, handle: () => handle! };
}

describe('the list window', () => {
  it("a VirtualList's rows say their place in the whole list, after a jump too", async () => {
    const m = mountVirtual(1000);
    await flushFrame();
    flush();
    const first = m.container.querySelector('[role=option]')!;
    expect(first.getAttribute('aria-posinset')).toBe('1');
    expect(first.getAttribute('aria-setsize')).toBe('1000');
    m.handle().scrollToIndex(600);
    const row = [...m.container.querySelectorAll('[role=option]')].find((r) => r.textContent === 'Row 600')!;
    expect(row.getAttribute('aria-posinset')).toBe('601');
    m.list.set(items(800));
    flush();
    await flushFrame();
    flush();
    expect(m.container.querySelector('[role=option]')!.getAttribute('aria-setsize')).toBe('800');
    m.dispose();
  });

  it("a RecycledList's rows (slots and keyed) say their place in the whole list", async () => {
    const list = state<readonly Item[]>(items(500));
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => (
        <RecycledList
          items={list.asReadonly}
          itemKey={(r) => r.id}
          role="listbox"
          positions
          slot={() => {
            const el = (<div role="option" />) as HTMLElement;
            return {
              el,
              show(r: Item | undefined) {
                el.textContent = r?.label ?? '';
              },
            };
          }}
          recycle={(r) => r.id !== 'r3'}
        >
          {(r) => <div role="option">{r.label}</div>}
        </RecycledList>
      ),
      container,
    );
    await flushFrame();
    flush();
    const byText = (t: string) => [...container.querySelectorAll('[role=option]')].find((r) => r.textContent === t)!;
    expect(byText('Row 2').getAttribute('aria-posinset')).toBe('3');
    expect(byText('Row 3').getAttribute('aria-posinset')).toBe('4');
    expect(byText('Row 3').getAttribute('aria-setsize')).toBe('500');
    dispose();
  });

  it('waits for its rest while paused: no reordering of a hidden list; the rest comes when it shows again', async () => {
    vi.useFakeTimers();
    const parent = new ViewScope();
    const list = state<readonly Item[]>(items(5000));
    let handle: VirtualListHandle | undefined;
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => (
        <RecycledList
          items={list.asReadonly}
          itemKey={(r) => r.id}
          ref={(h) => (handle = h)}
          slot={() => {
            const el = (<div />) as HTMLElement;
            return { el, show: (r: Item | undefined) => void (el.textContent = r?.label ?? '') };
          }}
        />
      ),
      container,
      parent,
    );
    await flushFrame();
    flush();
    const before = RecycledListStats.settles;
    handle!.scrollToIndex(3000);
    parent.setActive(false);
    vi.advanceTimersByTime(1000);
    expect(RecycledListStats.settles).toBe(before);
    parent.setActive(true);
    vi.advanceTimersByTime(1000);
    expect(RecycledListStats.settles).toBe(before + 1);
    dispose();
  });
});
