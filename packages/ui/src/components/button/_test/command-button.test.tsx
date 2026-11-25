/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it } from 'vitest';
import { KeyPlatform } from '../../../commands/keys.js';
import { CommandManager } from '../../../commands/manager.js';
import { flush, render } from '../../../dom/index.js';
import { Button } from '../button.js';
import { CommandButton } from '../command-button.js';

const mac = KeyPlatform.mac;
afterEach(() => {
  KeyPlatform.mac = mac;
  document.body.innerHTML = '';
});

function mount(content: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  render(content, container);
  return () => container.querySelector('button') as HTMLButtonElement;
}

describe('Button', () => {
  it('names an icon-only button by its title, unless it has an aria-label', () => {
    expect(mount(() => <Button icon={<span />} title="Reset" />)().getAttribute('aria-label')).toBe('Reset');
    const own = mount(() => <Button icon={<span />} title="Reset (Ctrl+R)" aria-label="Reset" />)();
    expect(own.getAttribute('aria-label')).toBe('Reset');
    expect(own.getAttribute('data-tip')).toBe('Reset (Ctrl+R)');
  });
});

describe('CommandButton', () => {
  it('an icon-only one is named by its title; the reason is its description, the key its shortcut', () => {
    KeyPlatform.mac = false;
    const ctx = new State({ loaded: false });
    const m = new CommandManager<{ loaded: boolean }>({ context: ctx.asReadonly, formats: { get: () => undefined } });
    m.register([
      {
        id: 'view.reset',
        title: 'Reset…',
        category: 'View',
        icon: 'rotate-ccw',
        shortcut: { key: 'r', mod: true, shift: true },
        when: (c) => (c.loaded ? true : 'Nothing is loaded.'),
        run: () => ({}),
      },
    ]);
    const button = mount(() => <CommandButton manager={m} command="view.reset" />);
    expect(button().getAttribute('aria-label')).toBe('Reset');
    expect(button().getAttribute('data-tip')).toBe('Nothing is loaded.');
    expect(button().getAttribute('aria-description')).toBe('Nothing is loaded.');
    expect(button().getAttribute('aria-keyshortcuts')).toBe('Control+Shift+R');
    ctx.set({ loaded: true });
    flush();
    expect(button().getAttribute('aria-label')).toBe('Reset');
    expect(button().getAttribute('data-tip')).toBe('Reset (Ctrl+Shift+R)');
    expect(button().hasAttribute('aria-description')).toBe(false);
  });

  it('a labelled one is named by its label', () => {
    const m = new CommandManager<object>({ context: new State({}).asReadonly, formats: { get: () => undefined } });
    m.register([{ id: 'a', title: 'Apply', category: 'View', icon: 'check', run: () => ({}) }]);
    const button = mount(() => <CommandButton manager={m} command="a" label="Go" />);
    expect(button().hasAttribute('aria-label')).toBe(false);
    expect(button().textContent).toBe('Go');
  });
});
