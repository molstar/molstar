/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import { flush, render } from '../../../dom/index.js';
import { CommandButton } from '../../button/command-button.js';
import { manager, mount, tick } from './fixture.js';

// What follows the app while it is shown: a row's reason it cannot run follows the context while the menu is open; a
// button for a lazy command shows it without registering inside its derived value (commands §4.2).

afterEach(() => {
  document.body.innerHTML = '';
});

describe('live rows and buttons', () => {
  it("a row's reason it cannot run follows the context while the menu is open", async () => {
    const ctx = new State<{ loaded: boolean }>({ loaded: false });
    const run = vi.fn(() => ({ message: 'ok' }));
    const m = manager(
      [
        {
          id: 'look.spin',
          title: 'Spin',
          category: 'Look',
          when: (c: any) => (c.loaded ? true : 'Nothing is loaded.'),
          run,
        },
      ],
      { context: ctx.asReadonly as any },
    );
    const t = mount({ manager: m });
    await tick();
    t.type('spin');
    const row = () => t.options()[0];
    expect(row().getAttribute('aria-disabled')).toBe('true');
    expect(row().textContent).toContain('Nothing is loaded.');
    ctx.set({ loaded: true });
    flush();
    expect(row().hasAttribute('aria-disabled')).toBe(false);
    t.key('Enter');
    await tick();
    expect(run).toHaveBeenCalledTimes(1);
  });

  it('a button for a lazy command shows it; making it does not run inside the derived lookup', async () => {
    const m = new CommandManager<object>({
      context: new State<object>({}).asReadonly,
      formats: { get: () => undefined },
    });
    const lazy = m.registerLazy({
      make: () => ({ commands: [{ id: 'lazy.one', title: 'One', category: 'Lazy', run: () => ({}) }], groups: [] }),
    });
    const get = vi.spyOn(m, 'get');
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(() => <CommandButton manager={m} command="lazy.one" label />, container);
    flush();
    expect(container.querySelector('button')?.textContent).toContain('One');
    expect(get).not.toHaveBeenCalled();
    // Stale and remade while shown: made in the microtask, the button still there.
    lazy.stale();
    await tick();
    flush();
    expect(container.querySelector('button')?.textContent).toContain('One');
    dispose();
  });
});
