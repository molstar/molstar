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
import type { CommandArgs, CommandDefinition, CommandFormats } from '../types.js';

// The invocation's API for the menu (commands §3.2–3.4): every param editable, the active field's choices, the sentence
// as it reads, the options and what differs from the defaults; and the manager's registration as one change (§4.1).

interface Ctx {
  log: string[];
}

const formats: CommandFormats<Ctx> = {
  get: (name) =>
    name === 'color'
      ? {
          hint: () => 'a name or #rrggbb',
          check: (v) =>
            /^(red|blue|#[0-9a-f]{6})$/i.test(v)
              ? { ok: true, value: v.toLowerCase() }
              : { ok: false, message: `"${v}" is not a color.` },
          complete: (text) =>
            [
              { value: 'red', label: 'red' },
              { value: 'blue', label: 'blue' },
            ].filter((c) => c.label.startsWith(text.toLowerCase())),
        }
      : undefined,
};

const REPRS: [string, string][] = [
  ['cartoon', 'Cartoon'],
  ['sticks', 'Sticks'],
  ['spacefill', 'Spacefill'],
];

function setup(extra: Partial<CommandDefinition<Ctx>> = {}) {
  const context = new State<Ctx>({ log: [] });
  const preview = { set: vi.fn(), clear: vi.fn() };
  const ran: CommandArgs[] = [];
  const checked: CommandArgs[] = [];
  const show: CommandDefinition<Ctx> = {
    id: 'look.show',
    title: 'Show…',
    category: 'Look',
    params: () => ({
      repr: PD.Select('cartoon', REPRS),
      color: PD.Optional(PD.Text('red')) as PD.Any,
      whole: PD.Boolean(true),
      opacity: PD.Numeric(1, { min: 0, max: 1, step: 0.1 }),
      size: PD.Numeric(0.2),
      quality: PD.Select<string>('auto', [
        ['auto', 'Auto'],
        ['high', 'High'],
      ]),
      tint: PD.Color(0xffffff as never),
      hidden: PD.Text('x', { isHidden: true }),
    }),
    fields: {
      repr: { before: 'as' },
      color: { before: 'in', format: 'color' },
      opacity: { before: 'at', presets: [0.25, [0.5, 'half'], 1, 2] },
    },
    preview: (args) => ({ ...args }),
    check: (args) => (checked.push(args), args.size === 99 ? 'Too big.' : undefined),
    run: (args) => (ran.push(args), {}),
    ...extra,
  };
  const manager = new CommandManager<Ctx>({ context, formats, preview });
  manager.register([show]);
  const r = manager.open('look.show');
  if (r.kind !== 'form') throw new Error('expected an invocation');
  return { manager, inv: r.invocation as CommandInvocation<Ctx>, preview, ran, checked };
}

describe('the sentence and the options', () => {
  it('splits the params: inline fields in the sentence, the rest as options (hidden ones nowhere)', () => {
    const { inv } = setup();
    expect(inv.sentence.value.map((f) => [f.id, f.kind, f.before])).toEqual([
      ['repr', 'select', 'as'],
      ['color', 'text', 'in'],
      ['whole', 'boolean', undefined],
      ['opacity', 'slider', 'at'],
      ['size', 'number', undefined],
      ['quality', 'select', undefined],
    ]);
    expect(inv.state.value.formOnly).toEqual(['tint']);
    expect(Object.keys(inv.options.value.params)).toEqual(['tint']);
    expect(inv.options.value.values).toEqual({ tint: 0xffffff });
    expect(inv.options.value.changed).toEqual([]);
  });

  it('tells each field empty, required, active, its default and whether it differs from it', () => {
    const { inv } = setup();
    const field = (id: string) => inv.sentence.value.find((f) => f.id === id);
    expect(field('repr')).toMatchObject({ empty: true, required: true, active: true, defaultDisplay: 'Cartoon' });
    expect(field('color')).toMatchObject({ empty: true, required: false, defaultValue: 'red', defaultDisplay: 'red' });
    inv.set('repr', 'cartoon');
    expect(field('repr')).toMatchObject({ value: 'cartoon', display: 'Cartoon', differs: false, active: false });
    inv.set('repr', 'sticks');
    expect(field('repr')).toMatchObject({ display: 'Sticks', differs: true });
    inv.set('whole', false);
    expect(field('whole')).toMatchObject({ display: 'off', differs: true });
    inv.set('opacity', 3);
    expect(field('opacity')).toMatchObject({ empty: true, error: 'Not a valid Opacity.' });
  });

  it('a form-only param is set from the invocation: in the run, the preview and the check; cleared back to its default', async () => {
    const { inv, preview, ran, checked } = setup();
    expect(inv.set('tint', 0xff0000)).toBe(true);
    expect(inv.state.value.optionValues).toEqual({ tint: 0xff0000 });
    expect(inv.options.value.values.tint).toBe(0xff0000);
    expect(inv.options.value.changed).toEqual([{ name: 'tint', label: 'Tint', value: 0xff0000 }]);
    // The active field is not moved by an option.
    expect(inv.state.value.field).toBe('repr');
    // Drawn once the required fields are set.
    for (const [k, v] of [
      ['repr', 'sticks'],
      ['whole', true],
      ['opacity', 1],
      ['size', 1],
    ] as const)
      inv.set(k, v);
    expect(preview.set).not.toHaveBeenCalled();
    inv.set('quality', 'high');
    expect(preview.set).toHaveBeenLastCalledWith(expect.objectContaining({ repr: 'sticks', tint: 0xff0000 }));
    inv.set('size', 99);
    await inv.confirm();
    expect(checked.at(-1)).toMatchObject({ tint: 0xff0000, size: 99 });
    expect(inv.state.value.error?.message).toBe('Too big.');
    expect(ran).toEqual([]);
    inv.set('size', 1);
    inv.clear('tint');
    expect(inv.options.value.changed).toEqual([]);
    inv.set('tint', 0x00ff00);
    await inv.confirm();
    expect(ran).toEqual([expect.objectContaining({ repr: 'sticks', size: 1, tint: 0x00ff00 })]);
  });

  it('a hidden param is set too; an unknown one is not', () => {
    const { inv } = setup();
    expect(inv.set('hidden', 'y')).toBe(true);
    expect(inv.args().hidden).toBe('y');
    expect(inv.set('nope', 1)).toBe(false);
    expect(inv.clear('nope')).toBe(false);
  });

  it('a form-only param that reshapes builds the params again, keeping what fits', () => {
    const context = new State<Ctx>({ log: [] });
    const manager = new CommandManager<Ctx>({ context, formats });
    manager.register([
      {
        id: 'export',
        title: 'Export…',
        category: 'Data',
        params: (_c, v) => ({
          name: PD.Text('model'),
          format: PD.Group({}, { label: 'format' }),
          kind: PD.Value<string>('cif', { isHidden: true }),
          ...(v.kind === 'bcif' ? { compress: PD.Boolean(true) } : {}),
        }),
        reshapes: ['kind'],
        run: () => ({}),
      },
    ]);
    const r = manager.open('export');
    if (r.kind !== 'form') throw new Error('expected an invocation');
    const inv = r.invocation;
    inv.set('name', 'out');
    expect(inv.state.value.fields.map((f) => f.id)).toEqual(['name']);
    inv.set('kind', 'bcif');
    expect(inv.state.value.fields.map((f) => f.id)).toEqual(['name', 'compress']);
    expect(inv.values()).toEqual({ name: 'out' });
    expect(inv.args()).toMatchObject({ name: 'out', kind: 'bcif', compress: true });
    inv.clear('kind');
    expect(inv.state.value.fields.map((f) => f.id)).toEqual(['name']);
  });

  it('views are told only when what they read changes: typing tells the choices, not the sentence or the options', () => {
    const { inv } = setup();
    const sentence = vi.fn();
    const options = vi.fn();
    const choices = vi.fn();
    const off = [inv.sentence.listen(sentence), inv.options.listen(options), inv.choices.listen(choices)];
    inv.type('st');
    inv.type('sti');
    expect(sentence).not.toHaveBeenCalled();
    expect(options).not.toHaveBeenCalled();
    expect(choices).toHaveBeenCalledTimes(2);
    for (const f of off) f();
  });
});

describe('choices of the active field', () => {
  it('a select: every option with the current one marked, filtered by what is typed, an exact name first', () => {
    const { inv } = setup();
    expect(inv.choices.value).toMatchObject({ field: 'repr', kind: 'select' });
    expect(inv.choices.value?.items.map((c) => [c.value, c.current])).toEqual([
      ['cartoon', false],
      ['sticks', false],
      ['spacefill', false],
    ]);
    inv.set('repr', 'sticks', { stay: true });
    expect(inv.choices.value?.items.filter((c) => c.current).map((c) => c.value)).toEqual(['sticks']);
    inv.type('s');
    expect(inv.choices.value?.items.map((c) => c.value)).toEqual(['sticks', 'spacefill']);
    expect(inv.state.value.pick).toBe(0);
    inv.type('spacefill');
    expect(inv.choices.value?.items[0]?.value).toBe('spacefill');
  });

  it('a boolean: its two values; a format: its completions, what is typed as read first, and its hint', () => {
    const { inv } = setup();
    inv.focus('whole');
    expect(inv.choices.value?.items.map((c) => c.label)).toEqual(['on', 'off']);
    inv.accept(1);
    expect(inv.values().whole).toBe(false);
    inv.focus('whole');
    expect(inv.choices.value?.items.map((c) => c.current)).toEqual([false, true]);
    inv.focus('color');
    expect(inv.choices.value).toMatchObject({ kind: 'text', hint: 'a name or #rrggbb' });
    expect(inv.choices.value?.items.map((c) => c.value)).toEqual(['red', 'blue']);
    inv.type('#FF0000');
    expect(inv.choices.value?.items).toEqual([
      { value: '#ff0000', label: '#ff0000', detail: 'as typed', current: false },
    ]);
    inv.accept();
    expect(inv.values().color).toBe('#ff0000');
  });

  it('a number: its presets that fit, its range for a slider; typing takes the number, not a preset', () => {
    const { inv } = setup();
    inv.focus('opacity');
    const choices = inv.choices.value;
    expect(choices?.range).toEqual({ min: 0, max: 1, step: 0.1 });
    expect(choices?.items.map((c) => [c.value, c.label])).toEqual([
      ['0.25', '0.25'],
      ['0.5', 'half'],
      ['1', '1'],
    ]);
    inv.type('60%');
    expect(inv.state.value.pick).toBe(-1);
    inv.accept();
    expect(inv.values().opacity).toBe(0.6);
    inv.focus('opacity');
    inv.move(1);
    inv.move(1);
    inv.accept();
    expect(inv.values().opacity).toBe(0.5);
    inv.focus('size');
    expect(inv.choices.value).toMatchObject({ kind: 'number', items: [] });
    expect(inv.choices.value?.range).toBeUndefined();
  });

  it('a text without a format lists nothing; the highlighted choice is previewed', () => {
    const { inv, preview } = setup({ params: () => ({ name: PD.Text(''), repr: PD.Select('cartoon', REPRS) }) });
    expect(inv.choices.value).toMatchObject({ field: 'name', items: [] });
    inv.set('name', 'x');
    inv.move(1);
    expect(preview.set).toHaveBeenLastCalledWith({ name: 'x', repr: 'cartoon' });
    inv.move(1);
    expect(preview.set).toHaveBeenLastCalledWith({ name: 'x', repr: 'sticks' });
    inv.point(0);
    expect(preview.set).toHaveBeenLastCalledWith({ name: 'x', repr: 'cartoon' });
  });
});

describe('moving between fields', () => {
  it('focus, Tab and Shift+Tab (wrapping), Backspace clearing then stepping back', async () => {
    const { inv } = setup();
    const active = () => inv.sentence.value.filter((f) => f.active).map((f) => f.id);
    expect(active()).toEqual(['repr']);
    await inv.tab(-1);
    expect(active()).toEqual(['quality']);
    await inv.tab(1);
    expect(active()).toEqual(['repr']);
    inv.focus('color');
    inv.type('blu');
    await inv.tab();
    expect(inv.values().color).toBe('blue');
    // On to the next empty required field after it.
    expect(inv.state.value.field).toBe('whole');
    inv.focus('color');
    expect(inv.back()).toBe(true);
    expect(inv.values().color).toBeUndefined();
    expect(inv.back()).toBe(true);
    expect(inv.state.value.field).toBe('repr');
    expect(inv.back()).toBe(false);
  });

  it('taking a value moves on to the next empty required field, then the next empty one', () => {
    const { inv } = setup();
    inv.accept(1);
    expect(inv.values().repr).toBe('sticks');
    expect(inv.state.value.field).toBe('whole');
    expect(inv.state.value.missing).toEqual(['whole', 'opacity', 'size', 'quality']);
  });
});

describe('registration (commands §4.1)', () => {
  const cmd = (id: string): CommandDefinition<Ctx> => ({ id, title: id, category: 'Test', run: () => ({}) });

  it('registering many commands, and taking them away, is one change each', () => {
    const manager = new CommandManager<Ctx>({ context: new State<Ctx>({ log: [] }), formats });
    const told = vi.fn();
    const off = manager.commands.listen(told);
    const drop = manager.register(Array.from({ length: 50 }, (_, i) => cmd(`c${i}`)));
    expect(told).toHaveBeenCalledTimes(1);
    drop();
    expect(told).toHaveBeenCalledTimes(2);
    drop();
    expect(told).toHaveBeenCalledTimes(2);
    off();
  });

  it('a lazy contribution made again is one change of the commands', () => {
    const manager = new CommandManager<Ctx>({ context: new State<Ctx>({ log: [] }), formats });
    let ids = ['a', 'b'];
    const lazy = manager.registerLazy({
      make: () => ({ commands: ids.map(cmd), groups: [{ id: `g.${ids[0]}`, title: 'G' }] }),
    });
    manager.get('a');
    const told = vi.fn();
    const off = manager.commands.listen(told);
    ids = ['c', 'd'];
    lazy.stale();
    manager.get('c');
    expect(told).toHaveBeenCalledTimes(1);
    expect(manager.commands.value.map((c) => c.id)).toEqual(['c', 'd']);
    expect(manager.groups.value.map((g) => g.id)).toEqual(['g.c']);
    off();
  });
});

describe('retiring the manager', () => {
  it('closes the open invocation and completes the outcomes', () => {
    const { manager, inv, preview } = setup();
    inv.set('tint', 1);
    const complete = vi.fn();
    manager.outcomes.subscribe({ complete });
    manager.retire();
    expect(inv.state.value.status).toBe('closed');
    expect(manager.active.value).toBeUndefined();
    expect(complete).toHaveBeenCalledTimes(1);
    expect(inv.set('tint', 2)).toBe(false);
    expect(preview.set).not.toHaveBeenCalled();
  });
});

describe('keys (commands §6.1)', () => {
  it('a key no command is on reads nothing of the event but its key', () => {
    const manager = new CommandManager<Ctx>({ context: new State<Ctx>({ log: [] }), formats });
    manager.register([{ id: 'a', title: 'A', category: 'T', shortcut: { key: 'a', alt: true }, run: () => ({}) }]);
    const composedPath = vi.fn(() => []);
    expect(manager.handleKey({ key: 'x', altKey: true, composedPath }, false)).toBe(false);
    expect(composedPath).not.toHaveBeenCalled();
    expect(manager.handleKey({ key: 'a', altKey: true, composedPath }, false)).toBe(true);
    expect(composedPath).toHaveBeenCalledTimes(1);
  });
});
