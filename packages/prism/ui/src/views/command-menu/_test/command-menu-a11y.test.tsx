/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { describe, expect, it } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import type { CommandDefinition } from '../../../commands/types.js';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandMenu } from '../command-menu.js';

const tick = () => new Promise((r) => setTimeout(r, 0));

function setup(extra = 0) {
  const ran: string[] = [];
  const command = (id: string, title: string, category: string): CommandDefinition<object> => ({
    id,
    title,
    category,
    run: () => {
      ran.push(id);
      return { message: title };
    },
  });
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  manager.register([
    command('view.reset', 'Reset the camera', 'View'),
    command('view.focus', 'Focus on selection', 'View'),
    command('look.cartoon', 'Show as cartoon', 'Look'),
    ...Array.from({ length: extra }, (_, i) => command(`more.c${i}`, `Command ${i}`, 'More')),
  ]);
  // What opened the menu (a button), then the menu, as the shell mounts it.
  const opener = document.createElement('button');
  opener.textContent = 'Commands';
  document.body.append(opener);
  const container = document.createElement('div');
  document.body.append(container);
  opener.focus();
  let dispose: (() => void) | undefined;
  const close = () => {
    dispose?.();
    dispose = undefined;
  };
  dispose = render(() => <CommandMenu manager={manager} onClose={close} />, container);
  const input = container.querySelector<HTMLInputElement>('input[role=combobox]')!;
  const key = (el: Element, k: string, init: KeyboardEventInit = {}) =>
    el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true, ...init }));
  const type = (text: string) => {
    input.value = text;
    input.dispatchEvent(new Event('input', { bubbles: true }));
  };
  const cleanup = () => {
    close();
    opener.remove();
    container.remove();
  };
  return { manager, ran, opener, container, input, key, type, close, cleanup, isOpen: () => !!dispose };
}

describe('CommandMenu accessibility (README §3.7)', () => {
  it('is a combobox whose input names the highlighted option; captions are not options', async () => {
    const t = setup();
    await tick();
    const { input, container } = t;
    const list = container.querySelector('[role=listbox]')!;
    expect(input.getAttribute('aria-controls')).toBe(list.id);
    expect(input.getAttribute('aria-autocomplete')).toBe('list');
    expect(input.getAttribute('aria-expanded')).toBe('true');
    // Every option has an id; the listbox holds options and presentational captions only.
    const options = [...list.querySelectorAll('[role=option]')];
    expect(options.length).toBe(3);
    expect(options.every((o) => o.id)).toBe(true);
    const captions = [...list.querySelectorAll('[role=presentation]')].map((c) => c.textContent);
    expect(captions).toEqual(['View', 'Look']);
    // The active descendant follows the highlight.
    expect(input.getAttribute('aria-activedescendant')).toBe(options[0].id);
    t.key(input, 'ArrowDown');
    expect(input.getAttribute('aria-activedescendant')).toBe(options[1].id);
    expect(options[1].getAttribute('aria-selected')).toBe('true');
    expect(options[0].getAttribute('aria-selected')).toBe('false');
    t.cleanup();
  });

  it('renders a screenful of a long list, and the highlighted option wherever it goes', async () => {
    const t = setup(100);
    await tick();
    const { input, container } = t;
    const options = () => [...container.querySelectorAll('[role=option]')];
    // 103 commands: a screenful is rendered, the rest is a spacer of their height.
    expect(options().length).toBeLessThan(40);
    // Up from the first wraps to the last: it is rendered, and the input names it.
    t.key(input, 'ArrowUp');
    const last = options().at(-1)!;
    expect(options().length).toBe(103);
    expect(last.textContent).toContain('Command 99');
    expect(input.getAttribute('aria-activedescendant')).toBe(last.id);
    expect(last.getAttribute('aria-selected')).toBe('true');
    // A search starts over from a screenful.
    t.type('command');
    expect(options().length).toBeLessThan(40);
    expect(input.getAttribute('aria-activedescendant')).toBe(options()[0].id);
    t.cleanup();
  });

  it('says what a search found, once the typing pauses, and when nothing matches', async () => {
    const t = setup();
    const status = t.container.querySelector('[role=status]')!;
    expect(status.getAttribute('aria-live')).toBe('polite');
    t.type('cartoon');
    await new Promise((r) => setTimeout(r, 450));
    expect(status.textContent).toBe('1 result');
    expect(t.input.getAttribute('aria-activedescendant')).toBe(t.container.querySelector('[role=option]')!.id);
    t.type('zzzz');
    await new Promise((r) => setTimeout(r, 450));
    expect(status.textContent).toBe('Nothing matches.');
    expect(t.input.hasAttribute('aria-activedescendant')).toBe(false);
    t.cleanup();
  });

  it('is modal: Tab stays inside, and closing gives the focus back to what opened it', async () => {
    const t = setup();
    await tick();
    const dialog = t.container.querySelector('[role=dialog]')!;
    expect(dialog.getAttribute('aria-modal')).toBe('true');
    expect(document.activeElement).toBe(t.input);
    // Tab and Shift+Tab stay inside (here the input is its one stop: the browser's default move is prevented).
    expect(t.key(t.input, 'Tab', { shiftKey: true })).toBe(false);
    expect(t.key(t.input, 'Tab')).toBe(false);
    expect(document.activeElement).toBe(t.input);
    expect(dialog.contains(document.activeElement)).toBe(true);
    // Escape closes it; the focus goes back to the button.
    t.key(document.activeElement!, 'Escape');
    expect(t.isOpen()).toBe(false);
    expect(document.activeElement).toBe(t.opener);
    t.cleanup();
  });

  it('gives the focus back after a command ran', async () => {
    const t = setup();
    await tick();
    t.key(t.input, 'Enter');
    await tick();
    expect(t.ran).toEqual(['view.reset']);
    expect(t.isOpen()).toBe(false);
    expect(document.activeElement).toBe(t.opener);
    t.cleanup();
  });

  it("offers a command's place in the UI: ⌥↵ (or its button) closes the menu and opens it", async () => {
    const t = setup();
    const opened: string[] = [];
    t.manager.register([
      {
        id: 'export.models',
        title: 'Export models…',
        category: 'Export',
        params: () => ({ name: PD.Text('') }),
        view: { title: 'Export › Models', open: () => opened.push('export') },
        run: () => ({ message: 'Exported.' }),
      },
    ]);
    t.manager.open('export.models');
    await tick();
    const form = t.container.querySelector<HTMLInputElement>('input[role=combobox]')!;
    expect(t.container.textContent).toContain('Export › Models');
    t.key(form, 'Enter', { altKey: true });
    await tick();
    expect(opened).toEqual(['export']);
    expect(t.isOpen()).toBe(false);
    t.cleanup();
  });
});
