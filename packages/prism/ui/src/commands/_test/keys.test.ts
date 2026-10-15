/** @vitest-environment jsdom */
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { observableState } from '../../reactive/observable-state.js';
import { formatChord } from '../keys.js';
import { CommandManager } from '../manager.js';
import type { CommandDefinition } from '../types.js';

interface Ctx {
  log: string[];
}

let dispose: (() => void) | undefined;
afterEach(() => {
  dispose?.();
  dispose = undefined;
  document.body.innerHTML = '';
});

function setup() {
  const context = observableState<Ctx>({ log: [] });
  const command = (id: string, shortcut: CommandDefinition<Ctx>['shortcut']): CommandDefinition<Ctx> => ({
    id,
    title: id,
    category: 'Scene',
    ...(shortcut ? { shortcut } : {}),
    run: (_v, c) => {
      c.log.push(id);
      return {};
    },
  });
  const manager = new CommandManager<Ctx>({ context: context.asReadonly, formats: { get: () => undefined } });
  manager.register([command('scene.edit', { key: 'e', alt: true }), command('story.save', { key: 's', mod: true })]);
  const onKey = (e: KeyboardEvent) => {
    if ((e.altKey || e.metaKey || e.ctrlKey) && manager.handleKey(e, true)) e.preventDefault();
  };
  window.addEventListener('keydown', onKey);
  dispose = () => window.removeEventListener('keydown', onKey);
  return context;
}

function press(el: Element, init: KeyboardEventInit): KeyboardEvent {
  const e = new KeyboardEvent('keydown', { bubbles: true, cancelable: true, ...init });
  el.dispatchEvent(e);
  return e;
}

function field(tag: 'textarea' | 'input' | 'div'): HTMLElement {
  const el = document.createElement(tag);
  // jsdom has no isContentEditable; the attribute is what the manager reads there.
  if (tag === 'div') el.setAttribute('contenteditable', 'true');
  document.body.append(el);
  return el;
}

describe('keybindings', () => {
  it("⌥E typed in a textarea ('ę' on Polish Pro) types the letter; it does not run the command", () => {
    const context = setup();
    const e = press(field('textarea'), { key: 'ę', code: 'KeyE', altKey: true });
    expect(e.defaultPrevented).toBe(false);
    expect(context.value.log).toEqual([]);
  });

  it('the US acute dead key (⌥E) in an editable element is not taken', () => {
    const context = setup();
    const e = press(field('div'), { key: 'Dead', code: 'KeyE', altKey: true });
    expect(e.defaultPrevented).toBe(false);
    expect(context.value.log).toEqual([]);
  });

  it('outside a text field, ⌥E runs it, by its code when the key is not a letter', async () => {
    const context = setup();
    expect(press(document.body, { key: 'Dead', code: 'KeyE', altKey: true }).defaultPrevented).toBe(true);
    await Promise.resolve();
    expect(press(document.body, { key: '´', code: 'KeyE', altKey: true }).defaultPrevented).toBe(true);
    expect(context.value.log).toEqual(['scene.edit', 'scene.edit']);
  });

  it('a ⌘ chord, which types nothing, runs from a text field', () => {
    const context = setup();
    expect(press(field('input'), { key: 's', code: 'KeyS', metaKey: true }).defaultPrevented).toBe(true);
    expect(context.value.log).toEqual(['story.save']);
  });

  it("matches by key first: on AZERTY ⌘A is code KeyQ and does not run ⌘Q's command", () => {
    const context = setup();
    press(document.body, { key: 'a', code: 'KeyS', metaKey: true });
    expect(context.value.log).toEqual([]);
  });

  it('holding the chord runs it once; a key in an IME composition is not a chord', () => {
    const context = setup();
    press(document.body, { key: '´', code: 'KeyE', altKey: true });
    press(document.body, { key: '´', code: 'KeyE', altKey: true, repeat: true });
    press(document.body, { key: 'Process', code: 'KeyE', altKey: true });
    press(document.body, { key: '´', code: 'KeyE', altKey: true, isComposing: true });
    expect(context.value.log).toEqual(['scene.edit']);
  });
});

describe('chord display', () => {
  it('shows ⌘ on macOS and Ctrl elsewhere', () => {
    expect(formatChord({ key: 'k', mod: true }, true)).toBe('⌘K');
    expect(formatChord({ key: 'k', mod: true }, false)).toBe('Ctrl+K');
    expect(formatChord({ key: 'p', mod: true, shift: true }, true)).toBe('⇧⌘P');
    expect(formatChord({ key: 'p', mod: true, shift: true }, false)).toBe('Ctrl+Shift+P');
    expect(formatChord({ key: 'ArrowRight', alt: true }, false)).toBe('Alt+→');
  });
});
