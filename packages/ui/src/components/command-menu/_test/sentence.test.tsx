/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import type { CommandInvocation } from '../../../commands/invocation.js';
import type { CommandArgs, CommandFormat } from '../../../commands/types.js';
import { DomProbe, flush } from '../../../dom/index.js';
import { manager, mount, tick } from './fixture.js';

// The sentence is the editor (commands §5.3): the caret in the active field's token, empty fields as they read, the
// list always the active field's, the Options fold's edits in the run, outcomes.

afterEach(() => {
  document.body.innerHTML = '';
});

const residues: CommandFormat<object> = {
  hint: () => 'A chain, a residue range: A 10-20',
  check: (text) =>
    /^[A-Z]( \d+(-\d+)?)?$/.test(text) ? { ok: true, value: text } : { ok: false, message: 'Not a selection.' },
  complete: (text) =>
    ['A', 'B', 'A 10-20']
      .filter((x) => x.toLowerCase().startsWith(text.toLowerCase()))
      .map((x) => ({ value: x, label: x, detail: 'chain' })),
};

function open(fail = false) {
  const ran: CommandArgs[] = [];
  const m = manager(
    [
      {
        id: 'look.show',
        title: 'Show as…',
        category: 'Look',
        params: () => ({
          what: PD.Text(''),
          type: PD.Select('cartoon', PD.arrayToOptions(['cartoon', 'surface', 'spacefill'])),
          size: PD.Optional(PD.Numeric(1, { min: 0, max: 5, step: 0.1 })) as PD.Any,
          quality: PD.Select('auto', PD.arrayToOptions(['auto', 'high', 'low'])),
          color: PD.Color(0xff0000 as any),
          smooth: PD.Boolean(true),
          style: PD.Group({ flat: PD.Boolean(false) }, { isFlat: true }),
        }),
        fields: {
          what: { format: 'residues' },
          type: { before: 'as' },
          size: { before: 'size', presets: [0.5, 1, 2] },
        },
        action: 'Apply',
        run: (args) => {
          if (fail) throw new Error('No structure.');
          ran.push(args);
          return { message: 'ok' };
        },
      },
    ],
    { formats: { get: (n) => (n === 'residues' ? residues : undefined) } },
  );
  const r = m.open('look.show');
  if (r.kind !== 'form') throw new Error(r.kind);
  const t = mount({ manager: m });
  flush();
  const token = (id: string) => t.container.querySelector<HTMLElement>(`[data-field="${id}"]`)!;
  return { m, t, ran, invocation: r.invocation as CommandInvocation<object>, token };
}

describe('the sentence', () => {
  it('puts the input in the active field’s token; others read as their value, name or default', () => {
    const { t, token } = open();
    // `what` is required and empty: active, the input in its token.
    expect(token('what').contains(t.input())).toBe(true);
    expect(t.input().getAttribute('aria-label')).toBe('What');
    // `type` (required, empty): its name, dashed; `size` (optional, empty): its default, muted; before-words quiet.
    const box = (id: string) => token(id).querySelector<HTMLElement>('[data-state]')!;
    expect(box('type').dataset.state).toBe('missing');
    expect(token('type').textContent).toBe('asType');
    expect(box('size').dataset.state).toBe('default');
    expect(token('size').textContent).toBe('size1');
    // Each token is a button named by its field and what it holds.
    expect(token('type').querySelector('button')?.getAttribute('aria-label')).toBe('Type: empty');
    // The fields past the sentence are not tokens.
    expect(token('color')).toBeNull();
  });

  it('typing inside the token filters its choices; Enter takes one and moves on, the input with it', () => {
    const { t, token, invocation } = open();
    // The format's completions and hint.
    expect(t.rows().map((r) => r.slice(0, 1))).toEqual(['A', 'B', 'A']);
    expect(t.container.textContent).toContain('A chain, a residue range');
    t.type('B');
    expect(t.input().value).toBe('B');
    expect(token('what').contains(t.input())).toBe(true);
    expect(t.input().parentElement?.getAttribute('data-text')).toBe('B');
    t.key('Enter');
    expect(invocation.state.value.values.what).toBe('B');
    // The next empty field is active: the input moved into its token, the list is its options.
    expect(token('type').contains(t.input())).toBe(true);
    expect(token('what').textContent).toBe('B');
    expect(t.rows().map((r) => r.trim())).toEqual(
      ['cartoon', 'surface', 'spacefill'].map((x) => x[0].toUpperCase() + x.slice(1)),
    );
    // A click on a choice takes it; the focus stays in the input.
    t.input().focus();
    t.options()[1].click();
    flush();
    expect(invocation.state.value.values.type).toBe('surface');
    expect(document.activeElement).toBe(t.input());
  });

  it('the list follows the active field: a number’s presets and its range, a select’s options', () => {
    const { t, invocation } = open();
    invocation.focus('size');
    flush();
    expect(t.rows()).toEqual(['0.5', '1', '2']);
    const range = t.container.querySelector<HTMLInputElement>('input[type=range]')!;
    expect(range).not.toBeNull();
    range.value = '3';
    range.dispatchEvent(new Event('input', { bubbles: true }));
    flush();
    expect(invocation.state.value.values.size).toBe(3);
    expect(invocation.state.value.field).toBe('size');
    // "60%" reads as 0.6.
    t.type('60%');
    t.key('Enter');
    expect(invocation.state.value.values.size).toBe(0.6);
    invocation.focus('quality');
    flush();
    expect(t.container.querySelector('input[type=range]')).toBeNull();
    expect(t.rows()).toEqual(['Auto', 'High', 'Low']);
  });

  it('a slider scrubs: each step sets the value, the same input keeps the focus', async () => {
    const { t, invocation } = open();
    invocation.focus('size');
    flush();
    const range = t.container.querySelector<HTMLInputElement>('input[type=range]')!;
    await tick();
    range.focus();
    for (const v of ['1.5', '2', '2.5', '3']) {
      range.value = v;
      range.dispatchEvent(new Event('input', { bubbles: true }));
      flush();
      await tick();
      expect(t.container.querySelector('input[type=range]')).toBe(range);
      expect(document.activeElement).toBe(range);
    }
    expect(invocation.state.value.values.size).toBe(3);
  });

  it('clicking a token, Tab, and ← / → at the text’s edge move between fields; ⌫ clears then steps back', () => {
    const { t, token, invocation } = open();
    token('quality').querySelector('button')!.click();
    flush();
    expect(invocation.state.value.field).toBe('quality');
    expect(token('quality').contains(t.input())).toBe(true);
    t.key('ArrowLeft');
    expect(invocation.state.value.field).toBe('size');
    t.key('ArrowRight');
    expect(invocation.state.value.field).toBe('quality');
    // Text typed: the arrows move the caret, not the field.
    t.type('hi');
    t.input().setSelectionRange(1, 1);
    t.key('ArrowLeft');
    expect(invocation.state.value.field).toBe('quality');
    t.type('');
    t.key('Tab', { shiftKey: true });
    expect(invocation.state.value.field).toBe('size');
    invocation.set('size', 2, { stay: true });
    flush();
    t.key('Backspace');
    expect(invocation.state.value.values.size).toBeUndefined();
    t.key('Backspace');
    expect(invocation.state.value.field).toBe('type');
  });

  it('a move of the pick renders the input and the two options it concerns, no more', () => {
    const { t, invocation } = open();
    invocation.focus('quality');
    flush();
    t.key('ArrowDown');
    let rendered = 0;
    DomProbe.on = (kind) => {
      if (kind === 'render') rendered++;
    };
    try {
      t.key('ArrowDown');
    } finally {
      DomProbe.on = undefined;
    }
    expect(rendered).toBeLessThanOrEqual(3);
  });

  it('⌫ at the very start goes back to the list; Esc too', () => {
    const { t, m } = open();
    t.key('Backspace');
    expect(m.active.value).toBeUndefined();
    flush();
    expect(t.input().getAttribute('aria-label')).toBe('Command');
  });

  it('a failure keeps the sentence open with its error; done closes the menu', async () => {
    const failing = open(true);
    failing.invocation.set('what', 'A');
    failing.invocation.set('type', 'surface');
    failing.invocation.set('quality', 'auto');
    failing.invocation.set('smooth', true);
    flush();
    failing.t.key('Enter');
    await tick();
    flush();
    expect(failing.m.active.value).toBe(failing.invocation);
    expect(failing.t.container.querySelector('[role=alert]')?.textContent).toBe('No structure.');
    expect(failing.t.closed()).toBe(0);
    document.body.innerHTML = '';

    const ok = open();
    ok.invocation.set('what', 'A');
    ok.invocation.set('type', 'surface');
    ok.invocation.set('quality', 'auto');
    ok.invocation.set('smooth', true);
    ok.t.container.querySelector<HTMLButtonElement>('button[class*=run]')!.click();
    await tick();
    expect(ok.ran.length).toBe(1);
    expect(ok.t.closed()).toBe(1);
  });
});

describe('the Options fold', () => {
  it('is closed, says how many; opened, its edits reach the run and its summary', async () => {
    const { t, invocation, ran } = open();
    const head = t.container.querySelector<HTMLButtonElement>('button[aria-controls^=msui-command-options]')!;
    expect(head.getAttribute('aria-expanded')).toBe('false');
    expect(head.textContent).toBe('Options2 at their defaults');
    head.click();
    flush();
    expect(head.getAttribute('aria-expanded')).toBe('true');
    // The form shows the options (color, style), not the sentence's fields.
    const body = t.container.querySelector(`#${head.getAttribute('aria-controls')}`)!;
    expect(body.textContent).toContain('Flat');
    expect(body.textContent).not.toContain('Quality');
    body.querySelector<HTMLElement>('[role=switch]')!.click();
    flush();
    expect(invocation.state.value.optionValues).toEqual({ style: { flat: true } });
    expect(head.textContent).toContain('style changed');
    invocation.set('what', 'A');
    invocation.set('type', 'surface');
    invocation.set('quality', 'low');
    invocation.set('smooth', false);
    await invocation.confirm();
    expect(ran[0]).toMatchObject({
      what: 'A',
      type: 'surface',
      quality: 'low',
      style: { flat: true },
      color: 0xff0000,
    });
  });

  it('follows an option set elsewhere', () => {
    const { t, invocation } = open();
    invocation.set('color', 0x00ff00);
    flush();
    expect(t.container.textContent).toContain('color #00ff00');
  });
});
