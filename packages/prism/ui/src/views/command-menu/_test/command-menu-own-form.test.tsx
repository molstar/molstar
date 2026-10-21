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
import { afterEach, describe, expect, it, vi } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import type { CommandArgs } from '../../../commands/types.js';
import { render } from '../../../dom/render.js';
import { flush } from '../../../dom/scope.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandMenu, type CommandFormView } from '../command-menu.js';

// A command's own form (commands.md §5.1): drawn under the sentence, the options of the fields it draws not listed, the
// keys of the input offered to it first, its form-only values (`setParam`) in the run's args.

afterEach(() => {
  document.body.innerHTML = '';
});

function open(form: CommandFormView) {
  const ran: CommandArgs[] = [];
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  manager.register([
    {
      id: 'look.show',
      title: 'Show as',
      category: 'Look',
      params: () => ({
        what: PD.Select('a', PD.arrayToOptions(['a', 'b'])),
        type: PD.Select('cartoon', PD.arrayToOptions(['cartoon', 'surface', 'spacefill'])),
        options: PD.Value<Record<string, unknown>>({}, { isHidden: true }),
      }),
      defaults: () => ({ what: 'a', type: 'cartoon' }),
      run: (args) => {
        ran.push(args);
        return {};
      },
    },
  ]);
  const opened = manager.open('look.show');
  if (opened.kind !== 'form') throw new Error(opened.kind);
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <CommandMenu manager={manager} forms={{ 'look.show': form }} onClose={vi.fn()} />, container);
  const input = () => container.querySelector<HTMLInputElement>('input[role=combobox]')!;
  const options = () => [...container.querySelectorAll<HTMLElement>('[role=option]')];
  const key = (k: string) => {
    const e = new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true });
    input().dispatchEvent(e);
    flush();
    return e;
  };
  return { manager, invocation: opened.invocation, container, input, options, key, ran };
}

describe('a command’s own form in the command menu', () => {
  it('draws under the sentence; its fields’ options are its own, the others listed', () => {
    const t = open({
      fields: ['type'],
      render: (p) => {
        p.invocation.focus('type');
        return <div data-own="">own</div>;
      },
    });
    flush();
    expect(t.container.querySelector('[data-own]')).not.toBeNull();
    // The type is active: its options are the form's, none listed.
    t.invocation.clear('type');
    flush();
    expect(t.options()).toEqual([]);
    // Another field: listed as before.
    t.invocation.focus('what');
    t.invocation.clear('what');
    flush();
    expect(t.options().map((o) => o.querySelector('span:nth-child(2)')?.textContent)).toEqual(['A', 'B']);
    // No "+ N more options": the hidden param is the form's.
    expect(t.container.textContent).not.toContain('more option');
  });

  it('sees the keys first; what it takes the menu does not act on', async () => {
    const seen: string[] = [];
    const t = open({
      render: (p) => {
        p.keys((e) => {
          seen.push(e.key);
          return e.key === 'ArrowDown';
        });
        return <div />;
      },
    });
    const pickBefore = t.invocation.state.value.pick;
    const down = t.key('ArrowDown');
    expect(down.defaultPrevented).toBe(true);
    expect(t.invocation.state.value.pick).toBe(pickBefore);
    t.key('Enter');
    await Promise.resolve();
    expect(seen).toEqual(['ArrowDown', 'Enter']);
    expect(t.ran.length).toBe(1);
  });

  it('Tab past the last field goes into its controls, and past them back to the input', () => {
    const t = open({
      fields: ['type'],
      render: () => (
        <div>
          <button type="button" tabIndex={-1} data-tile="">
            tile
          </button>
          <button type="button" data-size="">
            size
          </button>
          <button type="button" data-quality="">
            quality
          </button>
        </div>
      ),
    });
    flush();
    const size = t.container.querySelector<HTMLElement>('[data-size]')!;
    const quality = t.container.querySelector<HTMLElement>('[data-quality]')!;
    const tab = (shift = false) => {
      const e = new KeyboardEvent('keydown', { key: 'Tab', shiftKey: shift, bubbles: true, cancelable: true });
      (document.activeElement as HTMLElement).dispatchEvent(e);
      flush();
      return e;
    };
    t.input().focus();
    t.invocation.focus('what');
    flush();
    // Not the last field: the sentence's Tab.
    tab();
    expect(document.activeElement).toBe(t.input());
    expect(t.invocation.state.value.field).toBe('type');
    // The last field: into the form, past its tiles (not Tab stops), to its first control.
    expect(tab().defaultPrevented).toBe(true);
    expect(document.activeElement).toBe(size);
    // Between its controls: the browser's Tab (not taken).
    expect(tab().defaultPrevented).toBe(false);
    quality.focus();
    // Past the last: the input, on the first field.
    expect(tab().defaultPrevented).toBe(true);
    expect(document.activeElement).toBe(t.input());
    expect(t.invocation.state.value.field).toBe('what');
    // Shift+Tab on the first field: the form's last control; Shift+Tab before its first: the input, the last field.
    tab(true);
    expect(document.activeElement).toBe(quality);
    size.focus();
    tab(true);
    expect(document.activeElement).toBe(t.input());
    expect(t.invocation.state.value.field).toBe('type');
  });

  it('keeps Tab in the sentence when the form has no control to take the focus', () => {
    const t = open({ render: () => <div /> });
    t.input().focus();
    t.invocation.focus('type');
    flush();
    t.key('Tab');
    expect(document.activeElement).toBe(t.input());
    expect(t.invocation.state.value.field).toBe('what');
  });

  it('runs with the form-only values it set', async () => {
    const t = open({ render: () => <div /> });
    expect(t.invocation.setParam('options', { sizeFactor: 2 })).toBe(true);
    // Inline fields go through `set`; unknown names are refused.
    expect(t.invocation.setParam('type', 'surface')).toBe(false);
    expect(t.invocation.setParam('nope', 1)).toBe(false);
    expect(t.invocation.state.value.params).toEqual({ options: { sizeFactor: 2 } });
    await t.invocation.confirm();
    expect(t.ran[0]).toMatchObject({ what: 'a', type: 'cartoon', options: { sizeFactor: 2 } });
  });
});
