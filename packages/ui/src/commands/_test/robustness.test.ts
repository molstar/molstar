/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { describe, expect, it, vi } from 'vitest';
import type { CommandInvocation } from '../invocation.js';
import { CommandManager } from '../manager.js';
import type { CommandDefinition, CommandFormat } from '../types.js';

// A format's check that throws never strands an invocation (commands §2.2); typing is one change and lists the
// choices once; a set that changes nothing tells nothing; the slider's range is one object while it stands; a set of
// another value does not list the active field's choices again; a keydown checks staleness by a count (§6.1).

function open(command: CommandDefinition<object>, format?: CommandFormat<object>) {
  const m = new CommandManager<object>({
    context: new State<object>({}).asReadonly,
    formats: { get: (n) => (n === 'sel' ? format : undefined) },
  });
  m.register([command]);
  const r = m.open(command.id);
  if (r.kind !== 'form') throw new Error(r.kind);
  return { m, invocation: r.invocation as CommandInvocation<object> };
}

const show = (run = vi.fn(() => ({ message: 'ok' }))): CommandDefinition<object> => ({
  id: 'show',
  title: 'Show',
  category: 'Look',
  params: () => ({
    what: PD.Text(''),
    size: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }),
    color: PD.Color(0xff0000 as any),
  }),
  fields: { what: { format: 'sel' } },
  run,
});

describe('a format that throws', () => {
  it('reads as text it does not take, wherever it is checked; Enter leaves the invocation editable', async () => {
    let fail = false;
    const format: CommandFormat<object> = {
      check: (v) => {
        if (fail) throw new Error('The structure went away.');
        return { ok: true, value: v };
      },
      complete: () => [],
    };
    const run = vi.fn(() => ({ message: 'ok' }));
    const { invocation } = open(show(run), format);
    invocation.type('A');
    expect(invocation.accept()).toBe(true);
    fail = true;
    // Typing and the choices: no throw.
    invocation.focus('what');
    expect(() => invocation.type('B')).not.toThrow();
    expect(invocation.choices.value?.items).toEqual([]);
    expect(invocation.accept()).toBe(false);
    expect(invocation.state.value.error?.message).toBe('The structure went away.');
    // The run's last check: refused, and still editing.
    invocation.type('');
    invocation.set('size', 2, { stay: true });
    await invocation.confirm();
    expect(invocation.state.value.status).toBe('editing');
    expect(invocation.state.value.error).toEqual({ field: 'what', message: 'The structure went away.' });
    expect(run).not.toHaveBeenCalled();
    fail = false;
    await invocation.confirm();
    expect(run).toHaveBeenCalledTimes(1);
  });
});

describe('changes the views are told', () => {
  it('a keystroke is one change, and lists the choices once', () => {
    const complete = vi.fn((text: string) => [{ value: `${text}!`, label: `${text}!` }]);
    const format: CommandFormat<object> = { check: (v) => ({ ok: true, value: v }), complete };
    const { invocation } = open(show(), format);
    const seen = vi.fn();
    const unlisten = invocation.state.listen(seen);
    const choices = vi.fn();
    const unlistenChoices = invocation.choices.listen(choices);
    complete.mockClear();
    invocation.type('A');
    expect(seen).toHaveBeenCalledTimes(1);
    expect(invocation.state.value.pick).toBe(0);
    expect(complete).toHaveBeenCalledTimes(1);
    unlisten();
    unlistenChoices();
  });

  it('a set that changes nothing leaves the state as it was', () => {
    const { invocation } = open(show());
    invocation.set('size', 2, { stay: true });
    const before = invocation.state.value;
    invocation.set('size', 2, { stay: true });
    expect(invocation.state.value).toBe(before);
  });

  it('the slider range is one object while its field and bounds stand; another value lists no choices', () => {
    const complete = vi.fn(() => []);
    const format: CommandFormat<object> = { check: (v) => ({ ok: true, value: v }), complete };
    const { invocation } = open(show(), format);
    invocation.focus('size');
    const range = invocation.range.value;
    expect(range).toEqual({ field: 'size', min: 0, max: 5, step: 0.1 });
    for (const v of [1.5, 2, 2.5]) invocation.set('size', v, { stay: true });
    expect(invocation.range.value).toBe(range);
    // The active field is `what`; a set of an option (or another field) does not list its completions again.
    invocation.focus('what');
    const listed = invocation.choices.value;
    complete.mockClear();
    invocation.set('color', 0x00ff00);
    invocation.set('size', 3);
    expect(invocation.choices.value).toBe(listed);
    expect(complete).not.toHaveBeenCalled();
  });
});

describe('lazy commands and keys', () => {
  it('a keydown makes stale lazy commands only for a bound key; a made or disposed contribution is not stale', () => {
    const m = new CommandManager<object>({
      context: new State<object>({}).asReadonly,
      formats: { get: () => undefined },
      keymap: { app: { 'lazy.one': { key: 'l' } } },
    });
    const ran = vi.fn(() => ({}));
    const make = vi.fn(() => ({
      commands: [{ id: 'lazy.one', title: 'One', category: 'Lazy', run: ran }],
      groups: [],
    }));
    const lazy = m.registerLazy({ make });
    const key = (k: string) => m.handleKey({ key: k, defaultPrevented: false } as any, false);
    expect(key('x')).toBe(false);
    expect(make).not.toHaveBeenCalled();
    expect(key('l')).toBe(true);
    expect(make).toHaveBeenCalledTimes(1);
    expect(key('l')).toBe(true);
    expect(make).toHaveBeenCalledTimes(1);
    lazy.stale();
    lazy.dispose();
    key('l');
    expect(make).toHaveBeenCalledTimes(1);
  });
});
