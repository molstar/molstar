/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import type { CommandArgs } from '../../../commands/types.js';
import { flush } from '../../../dom/index.js';
import type { CommandFieldView, CommandFormView } from '../command-menu.js';
import { manager, mount } from './fixture.js';

// Custom UIs (commands §5.4): a command's own form (its fields' choices its own, the keys first, Tab into its controls
// and back, `options: false`), and a field editor by format.

afterEach(() => {
  document.body.innerHTML = '';
});

function open(props: { form?: CommandFormView; fieldViews?: Record<string, CommandFieldView> }) {
  const ran: CommandArgs[] = [];
  const m = manager([
    {
      id: 'look.show',
      title: 'Show as',
      category: 'Look',
      params: () => ({
        what: PD.Select('a', PD.arrayToOptions(['a', 'b'])),
        type: PD.Select('cartoon', PD.arrayToOptions(['cartoon', 'surface', 'spacefill'])),
        tint: PD.Text('red'),
        options: PD.Value<Record<string, unknown>>({}, { isHidden: true }),
        size: PD.Color(0 as any),
      }),
      fields: { tint: { format: 'color' } },
      defaults: () => ({ what: 'a', type: 'cartoon', tint: 'red' }),
      run: (args) => {
        ran.push(args);
        return {};
      },
    },
  ]);
  const opened = m.open('look.show');
  if (opened.kind !== 'form') throw new Error(opened.kind);
  const forms = props.form ? { 'look.show': props.form } : undefined;
  const t = mount({
    manager: m,
    ...(forms ? { forms } : {}),
    ...(props.fieldViews ? { fieldViews: props.fieldViews } : {}),
  });
  flush();
  return { m, t, invocation: opened.invocation, ran };
}

describe('a command’s own form', () => {
  it('draws under the sentence; its fields’ choices are its own, the others listed', () => {
    const { t, invocation } = open({ form: { fields: ['type'], render: () => <div data-own="">own</div> } });
    expect(t.container.querySelector('[data-own]')).not.toBeNull();
    invocation.focus('type');
    flush();
    expect(t.options()).toEqual([]);
    expect(t.input().getAttribute('aria-expanded')).toBe('false');
    invocation.focus('what');
    flush();
    expect(t.rows()).toEqual(['A', 'B']);
    // The Options fold stays for the params it does not draw (unless it says `options: false`).
    expect(t.container.textContent).toContain('Options');
  });

  it('with `options: false`, no Options fold', () => {
    const { t } = open({ form: { options: false, render: () => <div /> } });
    expect(t.container.textContent).not.toContain('Options');
  });

  it('sees the keys first; what it takes the menu does not act on', async () => {
    const seen: string[] = [];
    const { t, invocation, ran } = open({
      form: {
        options: false,
        render: (p) => {
          p.keys((e) => {
            seen.push(e.key);
            return e.key === 'ArrowDown';
          });
          return <div />;
        },
      },
    });
    const before = invocation.state.value.pick;
    expect(t.key('ArrowDown').defaultPrevented).toBe(true);
    expect(invocation.state.value.pick).toBe(before);
    t.key('Enter');
    await Promise.resolve();
    expect(seen).toEqual(['ArrowDown', 'Enter']);
    expect(ran.length).toBe(1);
  });

  it('Tab past the last field goes into its controls, and past them back to the input', () => {
    const { t, invocation } = open({
      form: {
        options: false,
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
      },
    });
    const size = t.container.querySelector<HTMLElement>('[data-size]')!;
    const quality = t.container.querySelector<HTMLElement>('[data-quality]')!;
    const tab = (shift = false) => t.key('Tab', { shiftKey: shift }, document.activeElement!);
    t.input().focus();
    invocation.focus('type');
    flush();
    // Not the last field: the sentence's Tab.
    tab();
    expect(document.activeElement).toBe(t.input());
    expect(invocation.state.value.field).toBe('tint');
    // The last field: into the form, past its tiles (not Tab stops), to its first control.
    expect(tab().defaultPrevented).toBe(true);
    expect(document.activeElement).toBe(size);
    // Between its controls: the browser's Tab.
    expect(tab().defaultPrevented).toBe(false);
    quality.focus();
    // Past the last: the input, on the first field.
    expect(tab().defaultPrevented).toBe(true);
    expect(document.activeElement).toBe(t.input());
    expect(invocation.state.value.field).toBe('what');
    // Shift+Tab on the first field: the form's last control; before its first: the input, on the last field.
    tab(true);
    expect(document.activeElement).toBe(quality);
    size.focus();
    tab(true);
    expect(document.activeElement).toBe(t.input());
    expect(invocation.state.value.field).toBe('tint');
  });

  it('keeps Tab in the sentence when nothing below takes the focus', () => {
    const { t, invocation } = open({ form: { options: false, render: () => <div /> } });
    t.input().focus();
    invocation.focus('tint');
    flush();
    t.key('Tab');
    expect(document.activeElement).toBe(t.input());
    expect(invocation.state.value.field).toBe('what');
  });

  it('runs with the params it set through the invocation, hidden ones too', async () => {
    const { invocation, ran } = open({ form: { render: () => <div /> } });
    expect(invocation.set('options', { sizeFactor: 2 })).toBe(true);
    expect(invocation.set('nope', 1)).toBe(false);
    await invocation.confirm();
    expect(ran[0]).toMatchObject({ what: 'a', type: 'cartoon', options: { sizeFactor: 2 } });
  });
});

describe('a field editor by format', () => {
  const swatches: CommandFieldView = {
    render: (p) => (
      <div data-swatches="">
        {['red', 'green', 'blue'].map((c) => (
          <button type="button" data-color={c} onClick={() => p.set(c)}>
            {c}
          </button>
        ))}
      </div>
    ),
  };

  it('replaces the choices while a field of its format is active, and sets the value', () => {
    const { t, invocation } = open({ fieldViews: { color: swatches } });
    invocation.focus('what');
    flush();
    expect(t.container.querySelector('[data-swatches]')).toBeNull();
    invocation.focus('tint');
    flush();
    expect(t.container.querySelector('[data-swatches]')).not.toBeNull();
    expect(t.options()).toEqual([]);
    t.container.querySelector<HTMLElement>('[data-color=blue]')!.click();
    flush();
    expect(invocation.state.value.values.tint).toBe('blue');
    invocation.focus('what');
    flush();
    expect(t.container.querySelector('[data-swatches]')).toBeNull();
  });

  it('sees the keys first while it shows', () => {
    const seen: string[] = [];
    const { t, invocation } = open({
      fieldViews: {
        color: {
          render: (p) => {
            p.keys((e) => {
              seen.push(e.key);
              if (e.key !== 'ArrowRight') return false;
              p.set('green');
              return true;
            });
            return <div />;
          },
        },
      },
    });
    invocation.focus('tint');
    flush();
    t.key('ArrowRight');
    expect(invocation.state.value.values.tint).toBe('green');
    invocation.focus('what');
    flush();
    t.key('ArrowDown');
    expect(seen).toEqual(['ArrowRight']);
  });
});
