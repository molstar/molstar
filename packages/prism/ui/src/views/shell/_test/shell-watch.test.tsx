/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { flush, RuntimeProbe } from '../../../dom/scope.js';
import { domView } from '../../../dom/view.js';
import { NotificationsManager } from '../../../managers/notifications.js';
import { ShellModel, type ShellViewContext } from '../../../shell/shell.js';
import { ShellView } from '../shell-view.js';

// The shell frame as a function of the layout (rendering.md §2.6): its variables and attributes are written by watches
// that follow the dock and shell states, each renders once for a burst, and the tab buttons follow the open tab by key.

let rendered = 0;
beforeEach(() => {
  rendered = 0;
  RuntimeProbe.on = (kind) => {
    if (kind === 'render') rendered++;
  };
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      disconnect() {}
      unobserve() {}
    },
  );
});
afterEach(() => {
  RuntimeProbe.on = undefined;
  vi.unstubAllGlobals();
  document.body.innerHTML = '';
});

function mount() {
  const shell = new ShellModel();
  shell.registerTabs([
    { id: 'structure', title: 'Structure', icon: 'house', render: domView<ShellViewContext>(() => <i>panel</i>) },
    { id: 'settings', title: 'Settings', icon: 'settings', render: domView<ShellViewContext>(() => <i>settings</i>) },
  ]);
  shell.show('structure');
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ShellView shell={shell} status={new NotificationsManager()} canvas={<canvas />} />, container);
  flush();
  const root = container.firstElementChild as HTMLElement;
  const viewport = root.querySelector<HTMLElement>('[data-msp-viewport]')!;
  const tab = (id: string) => container.querySelector<HTMLElement>(`[role=tab][data-tab="${id}"]`)!;
  return { shell, root, viewport, tab };
}

describe('the shell frame as a view of the layout', () => {
  it('writes the mode, the insets and the layers’ marks as the layout changes', () => {
    const v = mount();
    expect(v.root.getAttribute('data-mode')).toBe(v.shell.mode.value);
    const pinnedPx = () => v.viewport.style.getPropertyValue('--msp-canvas-pinned');
    const marked = () => v.viewport.querySelector('[data-pinned]') !== null;
    const start = v.shell.dock('right').state.value.pinned;
    expect(marked()).toBe(start);

    v.shell.setPinned('right', !start);
    flush();
    expect(marked()).toBe(!start);
    expect(pinnedPx() === '0px').toBe(start);

    v.shell.setPinned('right', start);
    flush();
    expect(marked()).toBe(start);
    expect(pinnedPx() === '0px').toBe(!start);
  });

  it('a burst of layout changes in one task renders each watch once', () => {
    const v = mount();
    flush();
    const start = rendered;
    v.shell.setSize('right', 400);
    v.shell.setSize('right', 410);
    v.shell.setSize('right', 420);
    flush();
    const burst = rendered - start;
    const again = rendered;
    v.shell.setSize('right', 430);
    flush();
    // Three sets in a task cost what one does.
    expect(burst).toBe(rendered - again);
  });

  it('the tab buttons follow the open tab by key', () => {
    const v = mount();
    expect(v.tab('structure').getAttribute('aria-selected')).toBe('true');
    expect(v.tab('settings').getAttribute('aria-selected')).toBe('false');
    v.shell.show('settings');
    flush();
    expect(v.tab('structure').getAttribute('aria-selected')).toBe('false');
    expect(v.tab('settings').getAttribute('aria-selected')).toBe('true');
  });
});
