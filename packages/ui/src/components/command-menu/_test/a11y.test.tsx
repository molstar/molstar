/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import type { CommandDefinition } from '../../../commands/types.js';
import { flush } from '../../../dom/index.js';
import { manager, mount, tick } from './fixture.js';

// ARIA (commands §5): a modal dialog; the input a combobox that keeps the focus and names the highlighted option; rows
// `option`s in a `listbox`, captions not; long lists a screenful at a time; a polite status for a search.

afterEach(() => {
  document.body.innerHTML = '';
});

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
  const m = manager([
    command('view.reset', 'Reset the camera', 'View'),
    command('view.focus', 'Focus on selection', 'View'),
    command('look.cartoon', 'Show as cartoon', 'Look'),
    ...Array.from({ length: extra }, (_, i) => command(`more.c${i}`, `Command ${i}`, 'More')),
  ]);
  const opener = document.createElement('button');
  document.body.append(opener);
  opener.focus();
  const t = mount({ manager: m, onClose: () => t.dispose() });
  return { m, ran, opener, t };
}

describe('the command menu’s accessibility', () => {
  it('is a combobox whose input names the highlighted option; captions are not options', async () => {
    const { t } = setup();
    await tick();
    const input = t.input();
    const list = t.container.querySelector('[role=listbox]')!;
    expect(input.getAttribute('aria-controls')).toBe(list.id);
    expect(input.getAttribute('aria-autocomplete')).toBe('list');
    expect(input.getAttribute('aria-expanded')).toBe('true');
    const options = [...list.querySelectorAll('[role=option]')];
    expect(options.length).toBe(3);
    expect(options.every((o) => o.id)).toBe(true);
    expect([...list.querySelectorAll('[role=presentation]')].map((c) => c.textContent)).toEqual(['View', 'Look']);
    expect(input.getAttribute('aria-activedescendant')).toBe(options[0].id);
    t.key('ArrowDown');
    expect(input.getAttribute('aria-activedescendant')).toBe(options[1].id);
    expect(options[1].getAttribute('aria-selected')).toBe('true');
    expect(options[0].getAttribute('aria-selected')).toBe('false');
  });

  it('renders a screenful of a long list, and the highlighted option wherever it goes', async () => {
    const { t } = setup(100);
    await tick();
    expect(t.options().length).toBeLessThan(40);
    t.key('ArrowUp');
    const last = t.options().at(-1)!;
    expect(t.options().length).toBe(103);
    expect(last.textContent).toContain('Command 99');
    expect(t.input().getAttribute('aria-activedescendant')).toBe(last.id);
    t.type('command');
    expect(t.options().length).toBeLessThan(40);
    expect(t.input().getAttribute('aria-activedescendant')).toBe(t.options()[0].id);
  });

  it('says what a search found, once the typing pauses, and when nothing matches', async () => {
    const { t } = setup();
    const status = t.container.querySelector('[role=status]')!;
    expect(status.getAttribute('aria-live')).toBe('polite');
    t.type('cartoon');
    await new Promise((r) => setTimeout(r, 450));
    expect(status.textContent).toBe('1 result');
    t.type('zzzz');
    await new Promise((r) => setTimeout(r, 450));
    expect(status.textContent).toBe('Nothing matches.');
    expect(t.input().hasAttribute('aria-activedescendant')).toBe(false);
  });

  it('is modal: Tab stays inside, and closing gives the focus back to what opened it', async () => {
    const { t, opener } = setup();
    await tick();
    const dialog = t.container.querySelector('[role=dialog]')!;
    expect(dialog.getAttribute('aria-modal')).toBe('true');
    expect(document.activeElement).toBe(t.input());
    expect(t.key('Tab', { shiftKey: true }).defaultPrevented).toBe(true);
    expect(dialog.contains(document.activeElement)).toBe(true);
    t.key('Escape', {}, document.activeElement!);
    expect(t.closed()).toBe(1);
    expect(document.activeElement).toBe(opener);
  });

  it('gives the focus back after a command ran', async () => {
    const { t, ran, opener } = setup();
    await tick();
    t.key('Enter');
    await tick();
    expect(ran).toEqual(['view.reset']);
    expect(t.closed()).toBe(1);
    expect(document.activeElement).toBe(opener);
  });

  it('names the sentence’s fields, and the input after the active one', async () => {
    const { t, m } = setup();
    m.register([
      {
        id: 'look.color',
        title: 'Color…',
        category: 'Look',
        params: () => ({ what: PD.Text(''), color: PD.Select('red', PD.arrayToOptions(['red', 'blue'])) }),
        fields: { color: { before: 'in' } },
        run: () => ({}),
      },
    ]);
    m.open('look.color');
    await tick();
    flush();
    const group = t.container.querySelector('[role=group]')!;
    expect(group.getAttribute('aria-label')).toBe('Color: arguments');
    expect(t.input().getAttribute('aria-label')).toBe('What');
    expect([...group.querySelectorAll('button[aria-label]')].map((b) => b.getAttribute('aria-label'))).toContain(
      'Color: empty',
    );
    expect(document.activeElement).toBe(t.input());
  });

  it('offers a command’s place in the UI: ⌥↵ closes the menu and opens it', async () => {
    const { t, m } = setup();
    const opened: string[] = [];
    m.register([
      {
        id: 'export.models',
        title: 'Export models…',
        category: 'Export',
        params: () => ({ name: PD.Text('') }),
        view: { title: 'Export › Models', open: () => opened.push('export') },
        run: () => ({ message: 'Exported.' }),
      },
    ]);
    m.open('export.models');
    await tick();
    flush();
    expect(t.container.textContent).toContain('Export › Models');
    t.key('Enter', { altKey: true });
    await tick();
    expect(opened).toEqual(['export']);
    expect(t.closed()).toBe(1);
  });
});
