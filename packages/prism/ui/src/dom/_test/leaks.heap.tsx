// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CommandManager } from '../../commands/manager.js';
import { openContextMenu } from '../../components/context-menu/context-menu.js';
import { RecycledList, type VirtualSlot } from '../../components/virtual-list/recycled-list.js';
import { VirtualList } from '../../components/virtual-list/virtual-list.js';
import { observableState } from '../../reactive/observable-state.js';
import { mockBlankMode, ShellStory } from '../../testing/shell.js';
import { CommandMenu } from '../../views/command-menu/command-menu.js';
import { describeScaling, heapUsedOutsideCode, isRetainedFlat } from '../../_test/perf.js';
import { render } from '../render.js';
import { flush, onCleanup } from '../scope.js';

// What a view attaches goes with it (rendering.md §2.4, the scoped helpers): representative views mounted and
// disposed again and again keep no heap (what stays after a collection over n and 4n cycles is flat) and leave no
// listener on the window, the document or its body (counted through a patched `addEventListener`). The positive
// control, a view that leaves a document listener, is caught.

// jsdom has none: the views' observers get one that observes nothing (as the viewer's shell tests do).
class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}
(globalThis as { ResizeObserver?: unknown }).ResizeObserver ??= NoResize;

/** Listeners on the shared targets now (window, document, body), net of removals. */
const shared = { count: 0 };
const proto = EventTarget.prototype;
const add = proto.addEventListener;
const remove = proto.removeEventListener;
const isShared = (t: EventTarget) => t === window || t === document || t === document.body;
const keyOf = (type: string, fn: unknown, o: unknown) =>
  `${type}|${typeof o === 'boolean' ? o : !!(o as { capture?: boolean } | undefined)?.capture}`;
const live = new Map<EventTarget, Map<unknown, Set<string>>>();
/** Where each listener on a shared target was added (to name a leak). */
const stacks = new Map<unknown, string>();
const leftOver = () => {
  const out: string[] = [];
  for (const byFn of live.values()) for (const [fn, keys] of byFn) if (keys.size) out.push(stacks.get(fn) ?? '?');
  return out;
};
proto.addEventListener = function (this: EventTarget, type: string, fn: unknown, o?: unknown) {
  if (isShared(this) && fn) {
    let byFn = live.get(this);
    if (!byFn) live.set(this, (byFn = new Map()));
    let keys = byFn.get(fn);
    if (!keys) byFn.set(fn, (keys = new Set()));
    const k = keyOf(type, fn, o);
    if (!keys.has(k)) {
      keys.add(k);
      shared.count++;
      stacks.set(fn, `${type}: ${(new Error().stack ?? '').split('\n').slice(2, 6).join(' / ')}`);
    }
  }
  return add.call(this, type, fn as EventListener, o as AddEventListenerOptions);
} as typeof proto.addEventListener;
proto.removeEventListener = function (this: EventTarget, type: string, fn: unknown, o?: unknown) {
  const byFn = isShared(this) ? live.get(this) : undefined;
  const keys = byFn?.get(fn);
  if (keys?.delete(keyOf(type, fn, o))) {
    shared.count--;
    // Forgotten once removed: the count must not keep a removed listener (and the view its closure holds) alive.
    if (keys.size === 0) {
      byFn!.delete(fn);
      stacks.delete(fn);
    }
  }
  return remove.call(this, type, fn as EventListener, o as EventListenerOptions);
} as typeof proto.removeEventListener;

/** A view's cycle: mounted, given a frame of work, disposed. */
type Cycle = (i: number) => void;

function mountCycle(view: () => Node): Cycle {
  return () => {
    const box = document.createElement('div');
    document.body.append(box);
    const dispose = render(view, box);
    flush();
    dispose();
    box.remove();
  };
}

const rows = Array.from({ length: 2000 }, (_, i) => ({ id: i, label: `Row ${i}` }));
const items = observableState<readonly { id: number; label: string }[]>(rows);

const cycles: Record<string, Cycle> = {
  // jsdom's own: 2,000 plain nodes made and removed (what jsdom keeps is not the runtime's).
  'plain DOM (jsdom alone)': () => {
    const box = document.createElement('div');
    document.body.append(box);
    for (let i = 0; i < 2000; i++) {
      const d = document.createElement('div');
      d.textContent = `Row ${i}`;
      box.append(d);
    }
    box.remove();
  },
  'a keyed list with rows': mountCycle(() => (
    <VirtualList items={items.asReadonly} itemKey={(r) => r.id} rowHeight={26}>
      {(r) => <div>{r.label}</div>}
    </VirtualList>
  )),
  'a recycled list with slots': mountCycle(() => (
    <RecycledList
      items={items.asReadonly}
      itemKey={(r) => r.id}
      rowHeight={26}
      slot={(): VirtualSlot<{ id: number; label: string }> => {
        const el = document.createElement('div');
        return {
          el,
          show: (r) => {
            el.textContent = r?.label ?? '';
          },
        };
      }}
    />
  )),
  'the command menu, opened and closed': mountCycle(() => {
    const manager = new CommandManager<object>({
      context: observableState<object>({}).asReadonly,
      formats: { get: () => undefined },
    });
    manager.register([{ id: 'a', title: 'Alpha', category: 'C', run: () => ({ message: '' }) }]);
    return <CommandMenu manager={manager} onClose={() => {}} />;
  }),
  'the shell, a tab shown': mountCycle(() => <ShellStory mode={mockBlankMode()} loaded openTab="help" />),
  'a context menu, opened and closed': () => {
    const root = document.createElement('div');
    root.className = 'msp-root';
    document.body.append(root);
    const close = openContextMenu(root, { x: 10, y: 10 }, [{ id: 'a', label: 'A', run: () => {} }]);
    close();
    root.remove();
  },
};

const gc = () => {
  (globalThis as { gc?: () => void }).gc!();
  (globalThis as { gc?: () => void }).gc!();
};
/** A frame between cycles: the microtasks a mount queued, and the animation frame its frame states asked for, run. */
const aFrame = () => new Promise((r) => setTimeout(r, 20));

/** `retainedScaling`, a frame between cycles. */
async function retained(n: number, cycle: Cycle) {
  const run = async (k: number) => {
    for (let i = 0; i < k; i++) {
      cycle(i);
      await aFrame();
    }
  };
  gc();
  const h0 = heapUsedOutsideCode();
  await run(n);
  gc();
  const h1 = heapUsedOutsideCode();
  await run(4 * n);
  gc();
  const h2 = heapUsedOutsideCode();
  return { short: h1 - h0, long: h2 - h1, scaled: h2 - h1 - (h1 - h0) };
}

describe('views leave nothing behind', () => {
  for (const [name, cycle] of Object.entries(cycles)) {
    it(name, async () => {
      for (let i = 0; i < 10; i++) {
        cycle(i);
        await aFrame();
      }
      const before = shared.count;
      const r = await retained(20, cycle);
      expect(
        shared.count - before,
        `listeners left on window, document or body: ${leftOver().slice(-3).join(' | ')}`,
      ).toBe(0);
      expect(isRetainedFlat(r, 20), describeScaling(r, 20)).toBe(true);
    });
  }

  it('a view that leaves a document listener is caught (positive control)', async () => {
    const leaky = mountCycle(() => {
      const el = document.createElement('div');
      // Not `on` (dom/scoped.ts): nothing removes it.
      document.addEventListener('keydown', () => el.remove());
      onCleanup(() => {});
      return el;
    });
    for (let i = 0; i < 5; i++) leaky(i);
    const before = shared.count;
    await retained(5, leaky);
    expect(shared.count - before).toBeGreaterThan(0);
  });
});
