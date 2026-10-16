/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import type { CommandDefinition, KeyBindings } from '../../../commands/types.js';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { HelpPanel } from '../help-panel.js';

// Help's keys (commands.md §6): the effective keys, conflicts, and remapping a command's key: recorded, cleared, reset.

const tick = () => new Promise((r) => setTimeout(r, 0));

afterEach(() => {
  document.body.innerHTML = '';
});

function setup() {
  const saved: KeyBindings[] = [];
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
    keymap: { app: { 'app.cleared': null }, store: { load: () => undefined, save: (b) => void saved.push(b) } },
  });
  const cmd = (
    id: string,
    title: string,
    extra: Partial<CommandDefinition<object>> = {},
  ): CommandDefinition<object> => ({
    id,
    title,
    category: 'Test',
    run: () => ({}),
    ...extra,
  });
  manager.register([
    cmd('next', 'Next snapshot', { shortcut: { key: 'ArrowRight' } }),
    cmd('reset', 'Reset camera', { shortcut: { key: 'r', alt: true } }),
    cmd('app.cleared', 'Cleared by the app', { shortcut: { key: 'x', alt: true } }),
    cmd('plain', 'No key at first'),
  ]);
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <HelpPanel commands={manager} bindings={[]} />, container);
  const commands = () => container.querySelector('[aria-label="Commands"]') as HTMLElement;
  const titles = () => [...commands().querySelectorAll('dt')].map((d) => d.textContent);
  const row = (title: string) =>
    [...commands().querySelectorAll('dt')].find((d) => d.textContent === title)?.parentElement as HTMLElement;
  const button = (scope: HTMLElement, label: string) =>
    [...scope.querySelectorAll('button')].find((b) => b.textContent === label) as HTMLButtonElement;
  return { manager, saved, container, commands, titles, row, button, dispose };
}

describe('Help’s keys', () => {
  it('lists the effective keys: the app’s removal is gone', () => {
    const p = setup();
    expect(p.titles()).toEqual(['Next snapshot', 'Reset camera']);
    p.dispose();
  });

  it('records a new key, says a conflict, clears and resets', async () => {
    const p = setup();
    p.button(p.commands(), 'Change keys').click();
    await tick();
    p.button(p.row('Reset camera'), 'Change').click();
    await tick();
    const recorder = p.row('Reset camera').querySelector('[role="button"]') as HTMLElement;
    expect(document.activeElement).toBe(recorder);
    // A modifier alone waits; the chord is taken.
    recorder.dispatchEvent(new KeyboardEvent('keydown', { key: 'Alt', altKey: true, bubbles: true }));
    expect(p.saved).toEqual([]);
    recorder.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    await tick();
    expect(p.saved.at(-1)).toEqual({ reset: { key: 'ArrowRight' } });
    // Two commands on → now: said in the conflicts and the live region.
    expect(p.container.querySelector('[aria-label="Shortcut conflicts"]')?.textContent).toContain(
      'Next snapshot and Reset camera',
    );
    expect(p.container.querySelector('[role="status"]')?.textContent).toContain('also used by Next snapshot');
    p.button(p.row('Reset camera'), 'Reset').click();
    await tick();
    expect(p.manager.shortcut(p.manager.get('reset')!)).toEqual({ key: 'r', alt: true });
    expect(p.container.querySelector('[aria-label="Shortcut conflicts"]')).toBeNull();
    p.button(p.row('Next snapshot'), 'Clear').click();
    await tick();
    expect(p.manager.shortcut(p.manager.get('next')!)).toBeUndefined();
    // Cleared by the person: still listed, to reset it.
    expect(p.row('Next snapshot').textContent).toContain('no key');
    p.dispose();
  });

  it('cancels a recording with Escape', async () => {
    const p = setup();
    p.button(p.commands(), 'Change keys').click();
    await tick();
    p.button(p.row('Reset camera'), 'Change').click();
    await tick();
    const recorder = p.row('Reset camera').querySelector('[role="button"]') as HTMLElement;
    const escape = new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, cancelable: true });
    recorder.dispatchEvent(escape);
    await tick();
    expect(escape.defaultPrevented).toBe(true);
    expect(p.row('Reset camera').querySelector('[role="button"]')).toBeNull();
    expect(p.saved).toEqual([]);
    p.dispose();
  });
});
