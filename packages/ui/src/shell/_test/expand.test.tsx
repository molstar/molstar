/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { flush } from '../../dom/index.js';
import { shellCommands } from '../shell-commands.js';
import { key, mountShell, probeView } from './fixture.js';

// Expand and Fullscreen (shell §3.9).

afterEach(() => {
  document.body.innerHTML = '';
});

const views = () => [
  probeView('home', { place: 'head', order: 0 }),
  probeView('layers', { order: 1 }),
  probeView('tool', { home: 'float' }),
];

describe('Expand', () => {
  it('shows only the main area: the rail, the docked column and the windows away and paused', () => {
    const m = mountShell(views());
    m.shell.show('home');
    m.shell.show('tool');
    flush();
    expect(m.rail()).not.toBeNull();
    m.shell.toggleExpanded(true);
    flush();
    expect(m.rail()).toBeNull();
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBe('false');
    expect(m.probes.home.scope?.active).toBe(false);
    expect(m.q('[data-window]')!.closest('[inert]')).not.toBeNull();
    expect(m.root.querySelector('[role="application"]')?.hasAttribute('data-expanded')).toBe(true);
    m.dispose();
  });

  it('Escape leaves it (after anything inside took it), and so does showing a docked view', () => {
    const m = mountShell(views());
    m.shell.show('home');
    m.shell.toggleExpanded(true);
    flush();
    key(m.main(), 'Escape');
    expect(m.shell.expanded.value).toBe(false);
    m.shell.toggleExpanded(true);
    m.shell.show('layers');
    flush();
    expect(m.shell.expanded.value).toBe(false);
    expect(m.frame('layers')).not.toBeNull();
    m.dispose();
  });

  it('`fill: page` puts the shell over the page while expanded', () => {
    const m = mountShell(views(), {}, { fill: 'page' });
    const root = m.root.querySelector('[role="application"]')!;
    expect(root.hasAttribute('data-fill')).toBe(false);
    m.shell.toggleExpanded(true);
    flush();
    expect(root.getAttribute('data-fill')).toBe('page');
    m.dispose();
  });
});

describe('Fullscreen', () => {
  it('is asked of the shell’s root; the result (and the browser’s own exit) is written back', async () => {
    const m = mountShell(views());
    const root = m.root.querySelector<HTMLElement>('[role="application"]')!;
    let current: Element | null = null;
    Object.defineProperty(document, 'fullscreenElement', { configurable: true, get: () => current });
    root.requestFullscreen = async () => {
      current = root;
    };
    document.exitFullscreen = async () => {
      current = null;
    };
    // The availability is decided as the view mounts: here, with `requestFullscreen` stubbed after, set it as it would.
    m.shell.setFullscreenAvailable(true);
    m.shell.requestFullscreen(true);
    await Promise.resolve();
    await Promise.resolve();
    expect(m.shell.fullscreen.value).toBe(true);
    current = null;
    document.dispatchEvent(new Event('fullscreenchange'));
    expect(m.shell.fullscreen.value).toBe(false);
    m.dispose();
    delete (document as any).fullscreenElement;
  });

  it('where it is not available, asking does nothing and its command says why', () => {
    const m = mountShell(views());
    m.shell.setFullscreenAvailable('Fullscreen is not allowed here');
    m.shell.requestFullscreen(true);
    expect(m.shell.fullscreenRequest.value).toBeUndefined();
    const cmd = shellCommands(m.shell).find((c) => c.id === 'view.fullscreen')!;
    expect(cmd.when?.(undefined)).toBe('Fullscreen is not allowed here');
    m.dispose();
  });
});
