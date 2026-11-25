/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it } from 'vitest';
import { flush, render } from '../../dom/index.js';
import { axeViolations } from '../../testing/axe.js';
import { HelpView } from '../help-view.js';
import { KeysModel } from '../keys-model.js';
import { KeyPlatform } from '../keys.js';
import { CommandManager } from '../manager.js';
import type { CommandDefinition, KeyBindings } from '../types.js';

// Help's keys (commands §6.6): the bindings, the effective keys, conflicts, and a command's key recorded, cleared and
// reset from the view; the live region says it.

const tick = () => new Promise((r) => setTimeout(r, 0));
const mac = KeyPlatform.mac;

let dispose: (() => void) | undefined;
afterEach(() => {
  dispose?.();
  dispose = undefined;
  KeyPlatform.mac = mac;
  document.body.innerHTML = '';
});

function setup() {
  KeyPlatform.mac = true;
  const saved: KeyBindings[] = [];
  const manager = new CommandManager<object>({
    context: new State<object>({}).asReadonly,
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
    cmd('next', 'Next step', { shortcut: { key: 'ArrowRight' } }),
    cmd('reset', 'Reset camera', { shortcut: { key: 'r', alt: true } }),
    cmd('app.cleared', 'Cleared by the app', { shortcut: { key: 'x', alt: true } }),
    cmd('plain', 'No key at first'),
  ]);
  const container = document.createElement('div');
  document.body.append(container);
  dispose = render(
    () => (
      <HelpView
        keys={new KeysModel(manager)}
        bindings={[{ title: 'Camera', bindings: [{ action: 'Rotate', input: ['Drag'], or: ['One finger'] }] }]}
      />
    ),
    container,
  );
  const commands = () => container.querySelector('[aria-label="Commands"]') as HTMLElement;
  const titles = () => [...commands().querySelectorAll('dt')].map((d) => d.textContent);
  const row = (title: string) =>
    [...commands().querySelectorAll('dt')].find((d) => d.textContent === title)?.parentElement as HTMLElement;
  const button = (scope: HTMLElement, label: string) =>
    [...scope.querySelectorAll('button')].find((b) => b.textContent === label) as HTMLButtonElement;
  return { manager, saved, container, commands, titles, row, button };
}

describe('the help view', () => {
  it('lists the app’s bindings, the keys everywhere and the commands’ effective keys', () => {
    const p = setup();
    expect(p.container.querySelector('[aria-label="Camera"]')?.textContent).toContain('RotateDragorOne finger');
    expect(p.container.querySelector('[aria-label="Everywhere"]')?.textContent).toContain('Esc');
    expect(p.titles()).toEqual(['Next step', 'Reset camera']);
    expect(p.row('Reset camera').querySelectorAll('kbd kbd')).toHaveLength(2);
  });

  it('records a new key, says a conflict, resets and clears', async () => {
    const p = setup();
    p.button(p.commands(), 'Change keys').click();
    flush();
    p.button(p.row('Reset camera'), 'Change').click();
    flush();
    await tick();
    const recorder = p.row('Reset camera').querySelector('[role="button"]') as HTMLElement;
    expect(document.activeElement).toBe(recorder);
    // A modifier alone waits; the chord is taken.
    recorder.dispatchEvent(new KeyboardEvent('keydown', { key: 'Alt', altKey: true, bubbles: true }));
    expect(p.saved).toEqual([]);
    recorder.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    flush();
    expect(p.saved.at(-1)).toEqual({ reset: { key: 'ArrowRight' } });
    expect(p.container.querySelector('[aria-label="Shortcut conflicts"]')?.textContent).toContain(
      'Next step and Reset camera',
    );
    expect(p.container.querySelector('[role="status"]')?.textContent).toContain('also used by Next step');
    p.button(p.row('Reset camera'), 'Reset').click();
    flush();
    expect(p.manager.shortcut(p.manager.get('reset')!)).toEqual({ key: 'r', alt: true });
    expect(p.container.querySelector('[aria-label="Shortcut conflicts"]')).toBeNull();
    p.button(p.row('Next step'), 'Clear').click();
    flush();
    expect(p.manager.shortcut(p.manager.get('next')!)).toBeUndefined();
    expect(p.row('Next step').textContent).toContain('no key');
  });

  it('cancels a recording with Escape, keeping the key from what is behind', async () => {
    const p = setup();
    p.button(p.commands(), 'Change keys').click();
    flush();
    p.button(p.row('Reset camera'), 'Change').click();
    flush();
    await tick();
    const recorder = p.row('Reset camera').querySelector('[role="button"]') as HTMLElement;
    const escape = new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, cancelable: true });
    recorder.dispatchEvent(escape);
    flush();
    expect(escape.defaultPrevented).toBe(true);
    expect(p.row('Reset camera').querySelector('[role="button"]')).toBeNull();
    expect(p.saved).toEqual([]);
  });

  it('lets Tab through the recorder: it leaves, and the blur cancels', async () => {
    const p = setup();
    p.button(p.commands(), 'Change keys').click();
    flush();
    p.button(p.row('Reset camera'), 'Change').click();
    flush();
    await tick();
    const recorder = p.row('Reset camera').querySelector('[role="button"]') as HTMLElement;
    const tab = new KeyboardEvent('keydown', { key: 'Tab', bubbles: true, cancelable: true });
    recorder.dispatchEvent(tab);
    expect(tab.defaultPrevented).toBe(false);
    expect(p.saved).toEqual([]);
    recorder.blur();
    flush();
    expect(p.row('Reset camera').querySelector('[role="button"]')).toBeNull();
  });

  it('gives a key from a choice that shows its placeholder again after a pick and a cancel', async () => {
    const p = setup();
    p.button(p.commands(), 'Change keys').click();
    flush();
    const choice = () => p.container.querySelector('button[aria-label="Give a key to a command"]') as HTMLElement;
    expect(choice().textContent).toContain('Give a key to…');
    choice().click();
    flush();
    const option = [...document.querySelectorAll('[role="option"]')].find((o) =>
      o.textContent?.includes('No key at first'),
    ) as HTMLElement;
    option.click();
    flush();
    await tick();
    // Recording it: the choice names it; cancelled, the placeholder again.
    expect(choice().textContent).toContain('No key at first');
    const recorder = p.row('No key at first').querySelector('[role="button"]') as HTMLElement;
    recorder.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, cancelable: true }));
    flush();
    expect(choice().textContent).toContain('Give a key to…');
    expect(document.querySelector('[role="listbox"]')).toBeNull();
  });

  it('speaks the app’s words', () => {
    KeyPlatform.mac = true;
    const manager = new CommandManager<object>({
      context: new State<object>({}).asReadonly,
      formats: { get: () => undefined },
    });
    const container = document.createElement('div');
    document.body.append(container);
    dispose = render(
      () => (
        <HelpView
          keys={new KeysModel(manager)}
          labels={{ everywhere: 'Partout', closeAction: 'Fermer', commands: 'Commandes' }}
        />
      ),
      container,
    );
    expect(container.querySelector('[aria-label="Partout"]')?.textContent).toContain('Fermer');
    expect(container.querySelector('[aria-label="Commandes"]')).not.toBeNull();
  });

  it('says a key change in the words of the keys model it is given', () => {
    const manager = new CommandManager<object>({
      context: new State<object>({}).asReadonly,
      formats: { get: () => undefined },
    });
    manager.register([{ id: 'next', title: 'Suivant', category: 'Test', shortcut: { key: 'n' }, run: () => ({}) }]);
    const keys = new KeysModel(manager, { labels: { cleared: (c) => `${c} : aucune touche.` } });
    const container = document.createElement('div');
    document.body.append(container);
    dispose = render(() => <HelpView keys={keys} />, container);
    keys.clear('next');
    flush();
    expect(container.querySelector('[role="status"]')?.textContent).toBe('Suivant : aucune touche.');
  });

  it('is accessible (axe), changing keys too', async () => {
    const p = setup();
    p.button(p.commands(), 'Change keys').click();
    flush();
    expect(await axeViolations(p.container)).toBe('');
  });
});
