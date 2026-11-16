/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { KeyPlatform } from '../../../commands/keys.js';
import { flush } from '../../../dom/index.js';
import type { CommandMenuAssistant } from '../command-menu.js';
import { manager, mount, tick } from './fixture.js';

// Prefixes, `?` help, the assistant's Ask row, and a menu opened for one command (commands §5.1–5.2, §5.5).

const mac = KeyPlatform.mac;
afterEach(() => {
  KeyPlatform.mac = mac;
  document.body.innerHTML = '';
});

describe('a prefix’s rows', () => {
  it('fills the text from a completion and keeps the menu open; a run closes it', async () => {
    const ran = vi.fn();
    const t = mount({
      manager: manager(),
      prefixes: [
        {
          char: 'say ',
          title: 'Say',
          items: (text) =>
            text.startsWith('hello')
              ? [{ id: 'run', label: `Run: ${text}`, run: () => void ran(text) }]
              : [{ id: 'fill', label: 'hello', run: () => ({ text: 'say hello world' }) }],
        },
      ],
    });
    t.type('say he');
    expect(t.rows()).toEqual(['hello']);
    t.key('Enter');
    expect(t.input().value).toBe('say hello world');
    expect(t.closed()).toBe(0);
    expect(t.rows()).toEqual(['Run: hello world']);
    t.key('Enter');
    expect(ran).toHaveBeenCalledWith('hello world');
    expect(t.closed()).toBe(1);
  });

  it('`?` lists the prefixes (picking one types it) and the commands with keys', () => {
    KeyPlatform.mac = true;
    const m = manager([
      { id: 'a', title: 'Reset', category: 'View', shortcut: { key: 'r', mod: true }, run: () => ({}) },
      { id: 'b', title: 'Quiet', category: 'View', run: () => ({}) },
    ]);
    const t = mount({ manager: m, prefixes: [{ char: '@', title: 'Scene', items: () => [] }] });
    t.type('?');
    expect(t.rows()).toEqual(['@  Scene', 'Reset⌘R']);
    t.key('Enter');
    expect(t.input().value).toBe('@');
    expect(t.input().placeholder).toBe('Scene');
  });
});

describe('the Ask row', () => {
  function withAssistant(ready = true) {
    const assistant: CommandMenuAssistant & { asked: string[]; setups: number } = {
      asked: [],
      setups: 0,
      ready: () => ready,
      propose: (text) => (text.includes('red') ? { command: 'look.color', tail: 'red' } : undefined),
      ask(text) {
        this.asked.push(text);
      },
      setup() {
        this.setups++;
      },
    };
    const m = manager([
      {
        id: 'look.color',
        title: 'Color',
        category: 'Look',
        alias: 'color',
        params: () => ({ color: PD.Select('blue', PD.arrayToOptions(['red', 'blue'])) }),
        run: () => ({}),
      },
    ]);
    return { m, assistant, t: mount({ manager: m, assistant }) };
  }

  it('↵ asks when nothing matches; a proposal opens its sentence, filled, and says where it came from', async () => {
    const { t, m, assistant } = withAssistant();
    t.type('paint it all red');
    const ask = t.container.querySelector<HTMLElement>('button[data-default]');
    expect(ask?.textContent).toContain('“paint it all red”');
    t.key('Enter');
    flush();
    expect(m.active.value?.state.value.values).toEqual({ color: 'red' });
    expect(t.container.textContent).toContain('Proposed from “paint it all red”');
    expect(assistant.asked).toEqual([]);
    m.close();
    flush();
    t.type('what is this');
    t.key('Enter', { metaKey: true, ctrlKey: true });
    expect(assistant.asked).toEqual(['what is this']);
    expect(t.closed()).toBe(1);
  });

  it('without a model, offers the setup', () => {
    const { t, assistant } = withAssistant(false);
    const row = t.container.querySelector<HTMLElement>('button[data-setup]')!;
    expect(row.textContent).toContain('Set up');
    row.click();
    expect(assistant.setups).toBe(1);
  });
});

describe('a menu opened for one command', () => {
  async function open(direct: boolean) {
    const m = manager([
      {
        id: 'look.show',
        title: 'Show as',
        category: 'Look',
        params: () => ({ type: PD.Text('cartoon') }),
        run: () => ({}),
      },
    ]);
    m.open('look.show', undefined, { type: 'surface' });
    const t = mount({ manager: m, direct });
    await tick();
    return { m, t };
  }

  it('closes with one Esc, and says so', async () => {
    const { t, m } = await open(true);
    expect(t.container.textContent).toContain('esc close');
    t.key('Escape');
    expect(m.active.value).toBeUndefined();
    expect(t.closed()).toBe(1);
  });

  it('opened from the list, its Esc goes back to the commands', async () => {
    const { t, m } = await open(false);
    expect(t.container.textContent).toContain('esc commands');
    t.key('Escape');
    expect(m.active.value).toBeUndefined();
    expect(t.closed()).toBe(0);
  });
});
