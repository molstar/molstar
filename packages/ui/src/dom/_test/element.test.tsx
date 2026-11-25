/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ReactiveModel, State } from '@molstar/core/reactive';
import {
  type Child,
  defineElements,
  domView,
  type ElementRoot,
  elementView,
  focusables,
  menuKeys,
  rememberFocus,
  runInScope,
  trapTab,
  TreeNavigation,
  ViewElement,
  ViewScope,
} from '../index.js';

// Elements, views, focus and tree keys (dom §7).

class Counter extends ReactiveModel {
  readonly count = new State(0);
  mounted = 0;
  protected override onMount(): void {
    this.mounted++;
  }
}

const roots: string[] = [];

class CounterElement extends ViewElement<Counter> {
  static readonly tag = 'core-test-counter';
  protected render(model: Counter): Child {
    return <b>{model.count}</b>;
  }
  protected override createModel(): Counter {
    return new Counter();
  }
  protected override onConnected(root: ElementRoot): void {
    roots.push(`connected ${root.nodeName}`);
  }
  protected override onReleased(root: ElementRoot): void {
    roots.push(`released ${root.nodeName}`);
  }
}

defineElements([CounterElement]);

const microtask = () => new Promise<void>((r) => queueMicrotask(r));

describe('ViewElement', () => {
  it('renders its model when connected, keeps it across a move, disposes it after a disconnect', async () => {
    roots.length = 0;
    const el = document.createElement(CounterElement.tag) as CounterElement;
    document.body.appendChild(el);
    const model = el.model!;
    expect(model.mounted).toBe(1);
    model.count.set(2);
    expect(el.textContent).toBe('2');
    const other = document.createElement('div');
    document.body.appendChild(other);
    other.appendChild(el);
    await microtask();
    expect(el.textContent).toBe('2');
    expect(model.count.observed).toBe(true);
    el.remove();
    await microtask();
    await microtask();
    expect(el.textContent).toBe('');
    expect(model.count.observed).toBe(false);
    expect(roots).toEqual(['connected #document', 'connected #document', 'released #document']);
  });

  it('renders a model given before the element was defined, and again for another model', () => {
    const el = document.createElement('div');
    document.body.appendChild(el);
    const first = new Counter();
    const view = document.createElement(CounterElement.tag) as CounterElement;
    view.model = first;
    el.appendChild(view);
    expect(first.mounted).toBe(0);
    const second = new Counter();
    second.count.set(7);
    view.model = second;
    expect(view.textContent).toBe('7');
    expect(first.count.observed).toBe(false);
    el.remove();
  });
});

describe('views', () => {
  it('domView renders under the host scope; elementView sets the context', () => {
    const host = document.createElement('div');
    const scope = new ViewScope();
    const label = new State('a');
    const dispose = domView<{ label: State<string> }>(({ ctx }) => <i>{ctx.label}</i>).mount(host, { label }, scope);
    scope.setActive(false);
    label.set('b');
    expect(host.textContent).toBe('a');
    scope.setActive(true);
    expect(host.textContent).toBe('b');
    dispose();
    const disposeEl = elementView<number>('x-thing').mount(host, 5, scope);
    expect((host.querySelector('x-thing') as any).ctx).toBe(5);
    disposeEl();
    expect(host.innerHTML).toBe('');
  });
});

describe('focus', () => {
  it('trapTab keeps Tab inside; rememberFocus gives it back; menuKeys moves between items', () => {
    const outside = document.createElement('button');
    const dialog = document.createElement('div');
    dialog.innerHTML = '<button id="a">a</button><button id="b">b</button><button hidden>c</button>';
    document.body.append(outside, dialog);
    outside.focus();
    const restore = rememberFocus();
    const [a, b] = focusables(dialog);
    expect(focusables(dialog).length).toBe(2);
    b.focus();
    const tab = new KeyboardEvent('keydown', { key: 'Tab', cancelable: true });
    trapTab(tab, dialog);
    expect(document.activeElement).toBe(a);
    restore(dialog);
    expect(document.activeElement).toBe(outside);
    const menu = document.createElement('div');
    menu.innerHTML = '<div role="menuitem" tabindex="-1">1</div><div role="menuitem" tabindex="-1">2</div>';
    document.body.append(menu);
    const items = menu.querySelectorAll<HTMLElement>('[role="menuitem"]');
    items[0].focus();
    const down = new KeyboardEvent('keydown', { key: 'ArrowDown', cancelable: true });
    Object.defineProperty(down, 'currentTarget', { value: menu });
    menuKeys(down);
    expect(document.activeElement).toBe(items[1]);
    outside.remove();
    dialog.remove();
    menu.remove();
  });

  it('TreeNavigation moves the roving focus by the rows model', () => {
    const rows = ['r1', 'r2', 'r3'];
    const scope = new ViewScope();
    const tree = document.createElement('div');
    tree.setAttribute('role', 'tree');
    document.body.append(tree);
    const nav = new TreeNavigation<string>({
      rows: () => rows,
      key: (r) => r,
      level: () => 1,
      scrollTo: () => {},
    });
    const els = rows.map((r) => {
      const el = document.createElement('div');
      el.setAttribute('role', 'treeitem');
      el.textContent = r;
      tree.append(el);
      return el;
    });
    runInScope(scope, () => {
      nav.tree(tree);
      rows.forEach((r, i) => nav.row(r)(els[i]));
    });
    nav.focusIndex(0);
    expect(document.activeElement).toBe(els[0]);
    els[0].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
    expect(document.activeElement).toBe(els[1]);
    expect(els[1].tabIndex).toBe(0);
    expect(els[0].tabIndex).toBe(-1);
    els[1].dispatchEvent(new KeyboardEvent('keydown', { key: 'End', bubbles: true }));
    expect(nav.active).toBe('r3');
    scope.dispose();
    tree.remove();
  });
});
