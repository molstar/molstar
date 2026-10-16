/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { domView } from '../../../dom/view.js';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { ShellModel, type ShellViewContext } from '../../../shell/shell.js';
import { NotificationsManager } from '../../../managers/notifications.js';
import { SectionList } from '../section-list.js';
import { ShellView } from '../shell-view.js';

// Hidden UI is inert and does no UI updates (styling.md §5): what the shell keeps mounted while hiding it (a folded
// section, a panel behind an overlay card) follows no model until it shows again.

describe('hidden shell content', () => {
  it('a folded section holds no subscriptions and its body shows the latest value when it opens', () => {
    const text = observableState('a');
    const shell = new ShellModel();
    shell.registerSections([
      {
        slot: 'x',
        id: 'one',
        title: 'One',
        render: domView<ShellViewContext>(() => <i>{text}</i>),
      },
    ] as never);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <SectionList shell={shell} slot="x" collapsed />, container);
    const i = container.querySelector('i') as HTMLElement;
    const body = i.closest('[data-section] > div') as HTMLElement;
    expect(body.hidden).toBe(true);
    expect(body.querySelector('[inert]')).not.toBeNull();
    expect(i.textContent).toBe('');

    text.set('b');
    expect(i.textContent).toBe('');

    (container.querySelector('button[aria-expanded]') as HTMLElement).click();
    expect(body.hidden).toBe(false);
    expect(i.textContent).toBe('b');
    expect(body.querySelector('[inert]')).toBeNull();
    text.set('c');
    expect(i.textContent).toBe('c');
    container.remove();
  });

  it('a panel behind an overlay card is inert and paused, and back as it was when the card closes', () => {
    vi.stubGlobal(
      'ResizeObserver',
      class {
        observe() {}
        disconnect() {}
        unobserve() {}
      },
    );
    const text = observableState('a');
    const shell = new ShellModel();
    shell.registerTabs([
      {
        id: 'structure',
        title: 'Structure',
        icon: 'house',
        render: domView<ShellViewContext>(() => <i>{text}</i>),
      },
    ]);
    shell.registerOverlays([{ id: 'card', title: 'Card', render: domView<ShellViewContext>(() => <b>card</b>) }]);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ShellView shell={shell} status={new NotificationsManager()} canvas={<canvas />} />, container);
    const panel = () => container.querySelector('[role="tabpanel"]') as HTMLElement;
    const i = container.querySelector('i') as HTMLElement;
    expect(i.textContent).toBe('a');
    expect(panel().hasAttribute('inert')).toBe(false);

    shell.openOverlay('card');
    expect(panel().hasAttribute('inert')).toBe(true);
    text.set('b');
    expect(i.textContent).toBe('a');

    shell.closeOverlay();
    expect(panel().hasAttribute('inert')).toBe(false);
    expect(container.querySelector('i')).toBe(i);
    expect(i.textContent).toBe('b');

    // Closed: the panel is not there at all.
    shell.toggle('right');
    expect(container.querySelector('[role="tabpanel"]')).toBeNull();
    container.remove();
    vi.unstubAllGlobals();
  });
});
