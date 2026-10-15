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
import { domView } from '../../../dom/view.js';
import { render } from '../../../dom/render.js';
import { ShellModel, type ShellViewContext } from '../../../shell/shell.js';
import { StatusModel } from '../../../shell/status.js';
import { ShellView } from '../shell-view.js';

// Expanded (`view.expand`, navigation.expanded): the canvas fills the viewer; the rail, the panels and the band are
// left out, the canvas stack stays; Escape or opening a panel brings them back.

beforeEach(() => {
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      disconnect() {}
      unobserve() {}
    },
  );
});
afterEach(() => vi.unstubAllGlobals());

function mount() {
  const shell = new ShellModel();
  shell.registerTabs([
    { id: 'structure', title: 'Structure', icon: 'house', render: domView<ShellViewContext>(() => <i>panel</i>) },
    {
      id: 'sequence',
      title: 'Sequence',
      icon: 'sequence',
      dock: 'top',
      frame: 'none',
      render: domView<ShellViewContext>(() => <s>band</s>),
    },
  ]);
  shell.registerToolbar([{ id: 'expand', place: 'canvas', group: 'canvas', icon: 'maximize-2', title: 'Expand' }]);
  shell.show('sequence');
  shell.show('structure');
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <ShellView shell={shell} status={new StatusModel()} canvas={<canvas />} />, container);
  const root = container.firstElementChild as HTMLElement;
  return {
    shell,
    root,
    has: (selector: string) => container.querySelector(selector) !== null,
    dispose: () => {
      dispose();
      container.remove();
    },
  };
}

describe('the expanded shell', () => {
  it('leaves out the rail, the panel and the band, and keeps the canvas stack', () => {
    const v = mount();
    expect(v.has('i')).toBe(true);
    expect(v.has('s')).toBe(true);
    expect(v.has('[aria-label="Expand"]')).toBe(true);

    v.shell.setExpanded(true);
    expect(v.root.hasAttribute('data-expanded')).toBe(true);
    expect(v.has('i')).toBe(false);
    expect(v.has('s')).toBe(false);
    expect(v.has('[role="tablist"]')).toBe(false);
    // The way back stays.
    expect(v.has('[aria-label="Expand"]')).toBe(true);

    v.shell.setExpanded(false);
    expect(v.root.hasAttribute('data-expanded')).toBe(false);
    expect(v.has('i')).toBe(true);
    expect(v.has('s')).toBe(true);
    v.dispose();
  });

  it('Escape and opening a panel bring the chrome back', () => {
    const v = mount();
    v.shell.setExpanded(true);
    v.root.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(v.shell.expanded.value).toBe(false);
    // Escape again is the panel's as before (it closes).
    expect(v.has('i')).toBe(true);

    v.shell.setExpanded(true);
    v.shell.show('structure');
    expect(v.shell.expanded.value).toBe(false);
    expect(v.has('i')).toBe(true);
    v.dispose();
  });
});
