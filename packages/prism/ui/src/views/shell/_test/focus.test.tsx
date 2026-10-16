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

// The keyboard in the shell (accessibility.md): the rail's tabs are a tab list, a closed panel gives the focus to its
// tab, an overlay card opened from the keyboard takes the focus and gives it back when it closes; the one live region
// says what is announced.

const frame = () => new Promise((r) => setTimeout(r, 32));

function mount() {
  const shell = new ShellModel();
  shell.registerTabs([
    { id: 'home', title: 'Home', icon: 'house', render: domView<ShellViewContext>(() => <i>home</i>) },
    {
      id: 'structure',
      title: 'Structure',
      icon: 'box',
      render: domView<ShellViewContext>(() => <input aria-label="Filter" />),
    },
  ]);
  shell.registerOverlays([
    { id: 'card', title: 'Card', render: domView<ShellViewContext>(() => <input aria-label="Width" />) },
  ]);
  const status = new StatusModel();
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <ShellView
        shell={shell}
        status={status}
        canvas={<canvas />}
        context={() => <button type="button">Clear selection</button>}
      />
    ),
    container,
  );
  return {
    shell,
    status,
    container,
    tab: (title: string) => container.querySelector<HTMLElement>(`[role=tab][aria-label="${title}"]`)!,
    done: () => {
      dispose();
      container.remove();
    },
  };
}

describe('shell focus', () => {
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

  it('tabs control their panel; closing the panel from inside it focuses its tab', () => {
    const t = mount();
    t.shell.show('structure');
    const panel = t.container.querySelector<HTMLElement>('[role=tabpanel]')!;
    const tab = t.tab('Structure');
    expect(tab.getAttribute('aria-controls')).toBe(panel.id);
    expect(panel.querySelector('h2')?.textContent).toBe('Structure');
    panel.querySelector('input')!.focus();
    (panel.querySelector('button[data-tip="Close panel"]') as HTMLElement).click();
    expect(t.container.querySelector('[role=tabpanel]')).toBeNull();
    expect(document.activeElement).toBe(tab);
    t.done();
  });

  it('a panel switched with the focus on the page leaves the focus there', () => {
    const t = mount();
    t.shell.show('structure');
    expect(document.activeElement).toBe(document.body);
    t.shell.show('home');
    expect(t.container.querySelector('[role=tabpanel] h2')?.textContent).toBe('Home');
    expect(document.activeElement).toBe(document.body);
    t.done();
  });

  it('a panel the focus left for the page (a click on the canvas) leaves it there when switched', async () => {
    const t = mount();
    t.shell.show('structure');
    const input = t.container.querySelector<HTMLElement>('[role=tabpanel] input')!;
    input.focus();
    input.blur();
    await Promise.resolve();
    t.shell.show('home');
    expect(document.activeElement).toBe(document.body);
    t.done();
  });

  it('an overlay card opened from the keyboard takes the focus, keeps Tab, and gives it back', async () => {
    const t = mount();
    const anchor = document.createElement('button');
    t.container.append(anchor);
    // From the keyboard: the focus ring shows on the button (jsdom's :focus-visible does not follow the keyboard).
    const matches = Element.prototype.matches;
    vi.spyOn(Element.prototype, 'matches').mockImplementation(function (this: Element, s: string) {
      return s === ':focus-visible' ? document.activeElement === this : matches.call(this, s);
    });
    anchor.focus();
    t.shell.openOverlay('card', anchor);
    const card = t.container.querySelector<HTMLElement>('[role=dialog]')!;
    expect(card.getAttribute('aria-modal')).toBe('true');
    expect(card.getAttribute('aria-labelledby')).toBe(card.querySelector('h2')!.id);
    await frame();
    expect(document.activeElement).toBe(card.querySelector('input'));
    t.shell.closeOverlay();
    expect(document.activeElement).toBe(anchor);
    vi.restoreAllMocks();
    t.done();
  });

  it('an overlay card stays open while the canvas is turned or its bottom bar used, and closes on a pointer down elsewhere', () => {
    const t = mount();
    t.shell.openOverlay('card');
    const down = (el: Element) =>
      el.dispatchEvent(new (window.PointerEvent ?? MouseEvent)('pointerdown', { bubbles: true, composed: true }));
    down(t.container.querySelector('[data-msp-canvas] canvas')!);
    expect(t.shell.overlay.value?.overlay.id).toBe('card');
    // The bar at the canvas's bottom acts on the scene too: the selection's Clear keeps the card (a PAE plot's pick).
    const clear = [...t.container.querySelectorAll('button')].find((b) => b.textContent === 'Clear selection')!;
    down(clear);
    expect(t.shell.overlay.value?.overlay.id).toBe('card');
    const outside = document.createElement('div');
    document.body.append(outside);
    down(outside);
    expect(t.shell.overlay.value).toBeUndefined();
    outside.remove();
    t.done();
  });

  it('says announced results in the live region, with no message on show', () => {
    const t = mount();
    const region = t.container.querySelector('[role=status]')!;
    t.status.announce('Selected 1 residue.');
    expect(region.textContent).toBe('Selected 1 residue.');
    expect(t.status.current.value).toBeUndefined();
    expect(t.status.history.value).toEqual([]);
    // The same words again are said again (the region's text changes).
    t.status.announce('Selected 1 residue.');
    expect(region.textContent).not.toBe('Selected 1 residue.');
    expect(region.textContent?.trim()).toBe('Selected 1 residue.');
    t.status.post('Loaded 1CBS.');
    expect(region.textContent).toBe('Loaded 1CBS.');
    t.done();
  });
});
