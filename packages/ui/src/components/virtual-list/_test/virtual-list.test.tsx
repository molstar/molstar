/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { VirtualList } from '../virtual-list.js';

afterEach(() => {
  document.body.innerHTML = '';
  vi.restoreAllMocks();
});

interface Item {
  readonly id: number;
  readonly label: string;
}

describe('VirtualList', () => {
  it('keys its rows by itemKey (two items of one key are one row)', () => {
    // The JSX runtime keeps `key` for itself, so a `key` prop never arrived and rows were keyed by the item itself.
    vi.spyOn(console, 'warn').mockImplementation(() => {});
    const items = state<readonly Item[]>([
      { id: 1, label: 'a' },
      { id: 1, label: 'b' },
    ]);
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => (
        <VirtualList items={items.asReadonly} itemKey={(r) => r.id}>
          {(r) => <div class="row">{r.label}</div>}
        </VirtualList>
      ),
      container,
    );
    expect(container.querySelectorAll('.row').length).toBe(1);
    dispose();
  });

  it('measures its rows from the rows made, not the rows planned', () => {
    // 200 rows planned (a screenful of 800 px at 24 px is 34 rows, plus overscan), but `For` has made only a few: the
    // window's height over the planned count would be far too small a row.
    const observers: { cb: (entries: any[]) => void; target?: Element }[] = [];
    vi.stubGlobal(
      'ResizeObserver',
      class {
        cb: (entries: any[]) => void;
        constructor(cb: (entries: any[]) => void) {
          this.cb = cb;
          observers.push({ cb });
        }
        observe(target: Element) {
          observers[observers.length - 1] = { cb: this.cb, target };
          observers.find((o) => o.cb === this.cb)!.target = target;
        }
        unobserve() {}
        disconnect() {}
      },
    );
    const items = state<readonly Item[]>(Array.from({ length: 500 }, (_, i) => ({ id: i, label: `r${i}` })));
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => (
        <VirtualList items={items.asReadonly} itemKey={(r) => r.id}>
          {(r) => <div class="row">{r.label}</div>}
        </VirtualList>
      ),
      container,
    );
    const win = container.querySelector('[class*="window"]') as HTMLElement;
    const rows = win.childElementCount;
    expect(rows).toBeGreaterThan(0);
    // Remove all but 3 rows (as if the slice ran out), then the observer reports 3 rows of 30 px.
    while (win.childElementCount > 3) win.lastElementChild!.remove();
    for (const o of observers) if (o.target === win) o.cb([{ contentRect: { height: 3 * 30 } }]);
    const content = container.querySelector('[class*="content"]') as HTMLElement;
    expect(content.style.height).toBe(`${500 * 30}px`);
    dispose();
    vi.unstubAllGlobals();
  });
});
