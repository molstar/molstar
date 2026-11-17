/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { domView, flush } from '../../dom/index.js';
import { findView } from '../layout.js';
import { idPart } from '../view-frame.js';
import { key, mountShell, probeView } from './fixture.js';

// The panes in the view (shell §2.1, §4, §6): the active pane's mark across rebuilds, a window's pane status, a tab's
// menu by where it is, names that follow and ids that do not clash, the shell's words, F6 by the shell's own marks, a
// drop a view took, and a window's move that leaves the docked trees alone.

afterEach(() => {
  document.body.innerHTML = '';
});

const views = () => [
  probeView('home', { place: 'head', order: 0 }),
  probeView('stage', { title: 'Stage B', stage: true, home: 'stage', instances: 'many' }),
  probeView('plot', { home: 'stage', stageSide: 'bottom' }),
  probeView('panel'),
];
const pane = (m: ReturnType<typeof mountShell>, id: string) => m.q(`[data-pane="${id}"]`);
const closeMenu = () => {
  document
    .querySelector<HTMLElement>('[role="menu"]')
    ?.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
  flush();
};

describe('the active pane’s mark', () => {
  it('stays on the active pane after the region is rebuilt (a pane shown beside two)', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    flush();
    expect(pane(m, 'content')!.hasAttribute('data-active')).toBe(true);
    m.shell.show('plot');
    flush();
    expect(m.qa('[data-pane]')).toHaveLength(3);
    expect(pane(m, 'content')!.hasAttribute('data-active')).toBe(true);
    expect(m.qa('[data-pane][data-active]')).toHaveLength(1);
    m.dispose();
  });
});

describe('a window’s pane', () => {
  it('is hovered and focused as a region’s pane is', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.float('stage#b');
    flush();
    const input = m.win('stage#b')!.querySelector<HTMLElement>('[data-input="stage"]')!;
    input.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    expect(m.shell.paneStatus('stage#b').value.hovered).toBe(true);
    m.win('stage#b')!.dispatchEvent(new PointerEvent('pointerleave'));
    expect(m.shell.paneStatus('stage#b').value.hovered).toBe(false);
    input.dispatchEvent(new FocusEvent('focusin', { bubbles: true }));
    flush();
    expect(m.shell.paneStatus('stage#b').value).toMatchObject({ focused: true, active: true });
    input.dispatchEvent(new FocusEvent('focusout', { bubbles: true, relatedTarget: m.main() }));
    expect(m.shell.paneStatus('stage#b').value.focused).toBe(false);
    m.dispose();
  });

  it('is named by its pane’s title as it changes', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.float('stage#b');
    flush();
    expect(m.win('stage#b')!.getAttribute('aria-label')).toBe('Stage B');
    const undo = m.shell.panePort('stage#b').setTitle('Stage D');
    flush();
    expect(m.win('stage#b')!.getAttribute('aria-label')).toBe('Stage D');
    undo();
    flush();
    expect(m.win('stage#b')!.getAttribute('aria-label')).toBe('Stage B');
    m.dispose();
  });
});

describe('element ids', () => {
  it('two shells on a page do not share a title id; refs are encoded one to one', () => {
    const a = mountShell(views());
    const b = mountShell(views());
    a.shell.show('panel');
    b.shell.show('panel');
    flush();
    const idA = a.frame('panel')!.querySelector('h2')!.id;
    const idB = b.frame('panel')!.querySelector('h2')!.id;
    expect(idA).not.toBe(idB);
    expect(document.querySelectorAll(`[id="${idA}"]`)).toHaveLength(1);
    expect(idPart('a#b')).not.toBe(idPart('a_b'));
    expect(idPart('stage-1')).toBe('stage-1');
    a.dispose();
    b.dispose();
  });
});

describe('a tab’s menu', () => {
  it('is the pane’s in the main area, the view’s when docked', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.move('plot', { kind: 'group', id: findView(m.shell.layout.value, 'stage#b')!.group.id });
    flush();
    const tab = m.q<HTMLElement>('[data-pane] [role="tab"]')!;
    tab.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true }));
    const items = () => [...document.querySelectorAll('[role="menuitem"]')].map((i) => i.textContent);
    expect(items()).toContain('Maximize');
    closeMenu();
    expect(document.querySelector('[role="menu"]')).toBeNull();
    m.shell.show('panel');
    m.shell.splitBelow('panel');
    // A lower group of two views: its tabs.
    m.shell.move('plot', { kind: 'group', id: findView(m.shell.layout.value, 'panel')!.group.id });
    flush();
    const docked = m.qa<HTMLElement>('[data-group] [role="tab"]').find((t) => !t.closest('[data-pane]'))!;
    key(docked, 'F10', { shiftKey: true });
    expect(items()).toContain('Float');
    expect(items()).not.toContain('Maximize');
    closeMenu();
    m.dispose();
  });
});

describe('the shell’s words (shell §6)', () => {
  it('the view reads the model’s labels', () => {
    const m = mountShell(
      views(),
      {
        labels: {
          mainArea: 'Zone principale',
          resize: 'Redimensionner',
          rail: 'Vues et outils',
          close: (t) => `Fermer ${t}`,
          workspace: 'Espace',
        },
      },
      { onDropFiles: () => {} },
    );
    m.shell.show('plot');
    flush();
    expect(m.main().getAttribute('aria-label')).toBe('Zone principale');
    expect(m.q('nav')!.getAttribute('aria-label')).toBe('Vues et outils');
    expect(m.q('[role="separator"]')!.getAttribute('aria-label')).toBe('Redimensionner');
    expect(m.q('[role="application"]')!.getAttribute('aria-label')).toBe('Espace');
    expect(
      m.frame('plot')!.querySelector('button[aria-label="Fermer Plot"], button[title="Fermer Plot"]'),
    ).not.toBeNull();
    m.dispose();
  });

  it('F6 finds the shell’s own regions, not by a label nor an app’s nav in a view', () => {
    vi.spyOn(Element.prototype, 'getClientRects').mockImplementation(function (this: Element) {
      return (this.isConnected ? [{}] : []) as unknown as DOMRectList;
    });
    const withNav = probeView('nav', {
      render: domView(() => (
        <nav aria-label="App nav">
          <button type="button" data-app-nav="">
            link
          </button>
        </nav>
      )),
    });
    const m = mountShell([probeView('layers'), withNav], { labels: { rail: 'Vues' } });
    m.shell.show('layers');
    flush();
    m.main().focus();
    const seen = new Set<string>();
    for (let i = 0; i < 4; i++) {
      key(document.activeElement!, 'F6');
      const a = document.activeElement!;
      expect(a.hasAttribute('data-app-nav')).toBe(false);
      seen.add(a.closest('[data-group]') ? 'group' : a.closest('nav') ? 'rail' : 'main');
    }
    expect(seen).toEqual(new Set(['group', 'rail', 'main']));
    m.dispose();
  });
});

describe('a file drop', () => {
  it('a view that took the drop (prevented its default) keeps it from the shell', async () => {
    const dropped: File[][] = [];
    const m = mountShell(views(), {}, { onDropFiles: (f) => dropped.push(f) });
    await Promise.resolve();
    const files = [new File(['x'], 'a.pdb')];
    const drop = (handled: boolean) => {
      const e = new Event('drop', { bubbles: true, cancelable: true }) as DragEvent;
      Object.defineProperty(e, 'dataTransfer', { value: { types: ['Files'], files } });
      if (handled) e.preventDefault();
      m.main().dispatchEvent(e);
    };
    drop(true);
    expect(dropped).toHaveLength(0);
    drop(false);
    expect(dropped).toHaveLength(1);
    m.dispose();
  });
});

describe('a window’s move', () => {
  it('writes no docked tree (they read nothing a window’s rectangle changes)', () => {
    const m = mountShell(views());
    m.shell.show('stage', { key: 'b' });
    m.shell.show('panel');
    m.shell.float('plot');
    flush();
    const win = findView(m.shell.layout.value, 'plot')!.window!;
    const shows = vi.spyOn(m.shell, 'shows');
    m.shell.setWindowRect(win.group.id, { ...win.rect, x: win.rect.x + 10 });
    flush();
    expect(shows).not.toHaveBeenCalled();
    expect(m.win('plot')!.style.left).toBe(`${win.rect.x + 10}px`);
    m.dispose();
  });
});
