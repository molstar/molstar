/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush, render } from '../../../dom/index.js';
import { TextInput } from '../text-input.js';

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(content: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(content, container);
  return { input: container.querySelector('input') as HTMLInputElement, dispose };
}

describe('TextInput', () => {
  it('follows its value and placeholder, and reports what is typed', () => {
    const value = new State('1tqn');
    const placeholder = new State('PDB ID');
    const onInput = vi.fn((v: string) => value.set(v));
    const { input } = mount(() => <TextInput value={value} placeholder={placeholder} onInput={onInput} />);
    expect(input.value).toBe('1tqn');
    expect(input.placeholder).toBe('PDB ID');
    value.set('7bv2');
    placeholder.set('');
    flush();
    expect(input.value).toBe('7bv2');
    expect(input.hasAttribute('placeholder')).toBe(false);
    input.value = '1cbs';
    input.dispatchEvent(new Event('input', { bubbles: true }));
    expect(onInput).toHaveBeenCalledWith('1cbs');
  });

  it('submits on Enter; a search field empties on Escape, and lets an empty one through', () => {
    const onSubmit = vi.fn();
    const onInput = vi.fn();
    // What it is in (a panel closing on Escape) sees only the Escape it lets through.
    const outer = vi.fn();
    const { input } = mount(() => (
      <div onKeyDown={outer}>
        <TextInput search value="abc" onInput={onInput} onSubmit={onSubmit} />
      </div>
    ));
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
    expect(onSubmit).toHaveBeenCalledWith('abc');
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(input.value).toBe('');
    expect(onInput).toHaveBeenCalledWith('');
    expect(outer).toHaveBeenCalledTimes(1); // the Enter
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(outer).toHaveBeenCalledTimes(2);
  });

  it('keeps what was typed in an uncontrolled field when the placeholder or disabled changes', () => {
    const placeholder = new State('PDB ID');
    const disabled = new State(false);
    const { input } = mount(() => <TextInput value="" placeholder={placeholder} disabled={disabled} />);
    input.value = '1tq';
    input.dispatchEvent(new Event('input', { bubbles: true }));
    placeholder.set('Search');
    flush();
    expect(input.value).toBe('1tq');
    disabled.set(true);
    flush();
    expect(input.value).toBe('1tq');
  });

  it('leaves Enter and Escape that end an IME composition to it', () => {
    const onSubmit = vi.fn();
    const onInput = vi.fn();
    const { input } = mount(() => <TextInput search value="かな" onInput={onInput} onSubmit={onSubmit} />);
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', isComposing: true, bubbles: true }));
    expect(onSubmit).not.toHaveBeenCalled();
    const escape = new KeyboardEvent('keydown', { key: 'Escape', isComposing: true, bubbles: true, cancelable: true });
    input.dispatchEvent(escape);
    expect(input.value).toBe('かな');
    expect(onInput).not.toHaveBeenCalled();
    expect(escape.defaultPrevented).toBe(false);
  });

  it('is disabled while its state says so, and a password hides its characters', () => {
    const disabled = new State(true);
    const { input } = mount(() => <TextInput type="password" disabled={disabled} />);
    expect(input.type).toBe('password');
    expect(input.disabled).toBe(true);
    disabled.set(false);
    flush();
    expect(input.disabled).toBe(false);
  });
});
