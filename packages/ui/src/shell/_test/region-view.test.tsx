/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush, Motion } from '../../dom/index.js';
import { axeViolations } from '../../testing/axe.js';
import { findView, mainNode } from '../layout.js';
import { key, mountShell, probeView } from './fixture.js';

// The main area's region in the view (shell §2.1): panes beside the app's content in one tree, a stage alone without a
// header, the active pane marked, a maximized pane alone (Escape back), the top overlay where the app places it, a stage's body
// kept as it moves between a pane and a window.

afterEach(() => {
  document.body.innerHTML = '';
});

const views = () => [
  probeView('home', { place: 'head', order: 0 }),
  probeView('stage', { title: 'Stage B', canvas: true, home: 'stage', instances: 'many' }),
  probeView('plot', { home: 'stage', region: { side: 'bottom' } }),
  probeView('seq', { title: 'Sequence', home: 'top' }),
];
const content = (m: ReturnType<typeof mountShell>) => m.q('[data-pane="content"]')!;
const pane = (m: ReturnType<typeof mountShell>, ref: string) =>
  m.q(`[data-msui-view="${ref}"]`)?.closest<HTMLElement>('[data-pane]') ?? null;

describe('the main area’s region', () => {
  it('one pane: its status is active, one of one; the shell draws no mark', () => {
    const m = mountShell(views());
    expect(m.shell.paneStatus('content').value).toMatchObject({ active: true, several: false });
    expect(content(m).hasAttribute('data-active')).toBe(false);
    m.dispose();
  });

  it('nothing docked: the content alone, no panes marked; a stage opens beside it, the content still the same node', () => {
    const m = mountShell(views());
    const app = m.q('[data-app-main]')!;
    expect(m.qa('[data-pane]').map((p) => p.dataset.pane)).toEqual(['content']);
    expect(m.q('[data-panes]')).toBeNull();
    m.shell.show('stage', { key: 'b' });
    flush();
    expect(m.qa('[data-pane]')).toHaveLength(2);
    // The app's content is the same element: a canvas in it is not remounted.
    expect(m.q('[data-app-main]')).toBe(app);
    expect(pane(m, 'stage#b')?.closest('[data-msui-shell-main]')).not.toBeNull();
    m.dispose();
  });

  it('a stage alone in its pane has no header; a panel beside it has one; ✕ takes it out', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.show('plot');
    flush();
    expect(m.frame('stage#b')!.querySelector('header')!.hidden).toBe(true);
    expect(m.frame('stage#b')!.getAttribute('aria-label')).toBe('Stage B');
    const header = m.frame('plot')!.querySelector('header')!;
    expect(header.hidden).toBe(false);
    header.querySelector<HTMLButtonElement>('button[data-tip="Close Plot"], button[aria-label="Close Plot"]')!.click();
    flush();
    expect(findView(m.shell.layout.value, 'plot')).toBeUndefined();
    m.dispose();
  });

  it('the pane pressed or focused is the active one, marked while there are several', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    flush();
    m.q('[data-input="stage"]')!.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    flush();
    expect(m.shell.activePane.value).toBe('stage#b');
    expect(pane(m, 'stage#b')!.hasAttribute('data-active')).toBe(true);
    expect(content(m).hasAttribute('data-active')).toBe(false);
    m.q('[data-app-main]')!.dispatchEvent(new FocusEvent('focusin', { bubbles: true }));
    flush();
    expect(m.shell.activePane.value).toBe('content');
    expect(content(m).hasAttribute('data-active')).toBe(true);
    // Each pane's status, for its content to react to: active, hovered, focused (the view feeds hover and focus).
    expect(m.shell.paneStatus('content').value).toMatchObject({ active: true, focused: true, several: true });
    expect(m.shell.paneStatus('stage#b').value).toMatchObject({ active: false, focused: false });
    m.q('[data-input="stage"]')!.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    expect(m.shell.paneStatus('stage#b').value.hovered).toBe(true);
    expect(m.shell.paneStatus('content').value.hovered).toBe(false);
    m.dispose();
  });

  it('panes show and hide at once (no slide: a panel’s motion, not a pane’s)', () => {
    const close = vi.fn(() => undefined);
    const was = Motion.close;
    Motion.close = close;
    try {
      const m = mountShell(views());
      m.shell.show('stage', { key: 'b' });
      m.shell.show('plot');
      flush();
      m.shell.maximize('stage#b');
      flush();
      m.shell.maximize(undefined);
      flush();
      expect(close).not.toHaveBeenCalled();
      expect(m.qa('[data-pane][data-reveal]')).toEqual([]);
      m.dispose();
    } finally {
      Motion.close = was;
    }
  });

  it('a splitter is the layout’s: pressing or focusing it changes no active pane', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.show('plot');
    flush();
    expect(m.shell.activePane.value).toBe('content');
    const splitter = pane(m, 'plot')!.querySelector<HTMLElement>(':scope > [role="separator"]')!;
    splitter.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    splitter.dispatchEvent(new FocusEvent('focusin', { bubbles: true }));
    expect(m.shell.activePane.value).toBe('content');
    m.dispose();
  });

  it('a maximized pane shows alone, the others away and paused; Escape brings them back', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.show('plot');
    flush();
    m.shell.maximize('stage#b');
    flush();
    expect(content(m).getAttribute('data-shown')).toBe('false');
    expect(pane(m, 'plot')!.getAttribute('data-shown')).toBe('false');
    expect(m.probes.plot.scope?.active).toBe(false);
    expect(pane(m, 'stage#b')!.getAttribute('data-shown')).toBeNull();
    expect(pane(m, 'stage#b')!.hasAttribute('data-flex')).toBe(true);
    // Alone on screen: no mark of being the active one.
    expect(m.q('[data-panes]')).toBeNull();
    key(m.q('[data-input="stage"]')!, 'Escape');
    expect(m.shell.maximized.value).toBeUndefined();
    expect(content(m).getAttribute('data-shown')).toBeNull();
    expect(m.probes.plot.scope?.active).toBe(true);
    m.dispose();
  });

  it('expand is the shell’s, maximize the pane’s: expand keeps every pane and hides the rail; maximize keeps the rail', () => {
    const m = mountShell(views());
    m.shell.show('home');
    m.shell.show('stage', { key: 'b' });
    flush();
    // A pane's own button: maximize. The pane alone, the shell's chrome (the rail) stays.
    m.shell.panePort('stage#b').toggleMaximized();
    flush();
    expect(content(m).getAttribute('data-shown')).toBe('false');
    expect(m.rail()).not.toBeNull();
    expect(m.shell.expanded.value).toBe(false);
    key(m.q('[data-input="stage"]')!, 'Escape');
    expect(m.shell.maximized.value).toBeUndefined();
    // Global: the shell's chrome away, the main area with all its panes.
    m.shell.toggleExpanded(true);
    flush();
    expect(m.rail()).toBeNull();
    expect(content(m).getAttribute('data-shown')).toBeNull();
    expect(pane(m, 'stage#b')!.getAttribute('data-shown')).toBeNull();
    // Both, from the app: Escape steps back one at a time (expand first).
    m.shell.panePort('stage#b').toggleMaximized();
    key(m.q('[data-input="stage"]')!, 'Escape');
    expect(m.shell.expanded.value).toBe(false);
    expect(m.shell.maximized.value).toBeDefined();
    m.dispose();
  });

  it('maximizes a window’s pane: out over the main area, its body kept; back in its window after', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.float('stage#b');
    flush();
    m.q<HTMLInputElement>('[data-input="stage"]')!.value = 'kept';
    m.shell.panePort('stage#b').toggleMaximized();
    flush();
    expect(m.shell.paneStatus('stage#b').value.maximized).toBe(true);
    const host = m.q('[data-maximized-pane="stage#b"]')!;
    expect(host.querySelector('[data-input="stage"]')).not.toBeNull();
    expect(m.rail()).not.toBeNull();
    m.shell.maximize(undefined);
    flush();
    expect(m.q('[data-maximized-pane]')).toBeNull();
    expect(m.win('stage#b')).not.toBeNull();
    expect(m.q<HTMLInputElement>('[data-input="stage"]')!.value).toBe('kept');
    expect(m.probes.stage.mounts).toBe(1);
    m.dispose();
  });

  it('a lined-up grid: a cell’s splitter lights its line across (the other column’s with it)', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.show('stage', { key: 'c' });
    m.shell.show('stage', { key: 'd' });
    m.shell.setArea(1000, 800);
    m.shell.arrangeGrid();
    flush();
    // The horizontal splitters of the region's two columns (one a line: on the second cell of each).
    const lines = m
      .qa('[data-msui-shell-main] [role="separator"][aria-orientation="horizontal"]')
      .filter((x) => !x.hidden);
    expect(lines).toHaveLength(2);
    lines[0].dispatchEvent(new PointerEvent('pointerenter'));
    expect(lines[1].hasAttribute('data-linked')).toBe(true);
    lines[0].dispatchEvent(new PointerEvent('pointerleave'));
    expect(lines[1].hasAttribute('data-linked')).toBe(false);
    m.dispose();
  });

  it('a press in a window makes its view the active pane', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.float('stage#b');
    flush();
    m.q('[data-input="stage"]')!.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect(m.shell.activePane.value).toBe('stage#b');
    expect(m.shell.paneStatus('stage#b').value.active).toBe(true);
    m.dispose();
  });

  it('a stage’s body is kept as it floats out and docks back beside the content', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    flush();
    m.q<HTMLInputElement>('[data-input="stage"]')!.value = 'kept';
    m.shell.float('stage#b');
    flush();
    expect(m.win('stage#b')).not.toBeNull();
    m.shell.split('stage#b', 'left');
    flush();
    expect(pane(m, 'stage#b')).not.toBeNull();
    expect(m.q<HTMLInputElement>('[data-input="stage"]')!.value).toBe('kept');
    expect(m.probes.stage.mounts).toBe(1);
    m.dispose();
  });

  it('the top overlay: the view’s body where the app places it (a band slot), no header', () => {
    const m = mountShell(
      views(),
      { top: 'overlay' },
      {
        main: ({ top }) => (
          <div data-app-main="">
            <div data-band-slot="">{top()}</div>
            scene
          </div>
        ),
      },
    );
    m.shell.show('seq');
    flush();
    expect(mainNode(m.shell.layout.value).top?.views).toEqual(['seq']);
    const frame = m.frame('seq')!;
    expect(frame.closest('[data-band-slot]')).not.toBeNull();
    expect(frame.querySelector('header')!.hidden).toBe(true);
    m.dispose();
  });

  it('the pane menu: maximize, float, close; a view elsewhere offers the main view', () => {
    const m = mountShell(views());
    m.shell.show('home');
    m.shell.show('stage', { key: 'b' });
    flush();
    const right = m.frame('home')!.querySelector('header')!;
    right.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true, clientX: 5, clientY: 5 }));
    flush();
    const labels = () => [...document.querySelectorAll('[role="menu"] [role^="menuitem"]')].map((r) => r.textContent);
    expect(labels()).toContain('Beside the main view');
    [...document.querySelectorAll<HTMLElement>('[role="menu"] [role^="menuitem"]')]
      .find((r) => r.textContent === 'Beside the main view')!
      .click();
    flush();
    expect(pane(m, 'home')).not.toBeNull();
    m.frame('home')!
      .querySelector('header')!
      .dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true, clientX: 5, clientY: 5 }));
    flush();
    expect(labels()).toEqual(expect.arrayContaining(['Maximize', 'Float', 'Close']));
    // Closed before the shell goes: an open menu's listeners would outlive the test (files share a worker).
    key(document.activeElement!, 'Escape');
    expect(document.querySelector('[role="menu"]')).toBeNull();
    m.dispose();
  });

  it('is accessible (axe): panes, a stage without a header, a panel', async () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.show('plot');
    flush();
    expect(await axeViolations(m.root, ['scrollable-region-focusable'])).toBe('');
    m.dispose();
  });
});
