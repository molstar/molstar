/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { observableState } from '../../reactive/observable-state.js';
import { CommandManager } from '../manager.js';
import type { CommandDefinition, CommandFormats, CommandOutcome, CommandResult } from '../types.js';
import { inline, type TestField } from './inline.js';

interface Ctx {
  current: string;
  selection: string[];
  working: boolean;
  alignments: number;
  log: string[];
}

const none: CommandFormats<Ctx> = { get: () => undefined };
const tick = () => new Promise((r) => setTimeout(r, 0));

function make(
  commands: CommandDefinition<Ctx>[],
  options: { formats?: CommandFormats<Ctx>; previewDeps?: (c: Ctx) => readonly unknown[] } = {},
) {
  const context = observableState<Ctx>({ current: 's1', selection: ['A:1'], working: false, alignments: 0, log: [] });
  const preview = { set: vi.fn(), clear: vi.fn() };
  const manager = new CommandManager<Ctx>({
    context: context.asReadonly,
    formats: options.formats ?? none,
    preview,
    ...(options.previewDeps ? { previewDeps: options.previewDeps } : {}),
  });
  const unregister = manager.register(commands);
  const outcomes: CommandOutcome<Ctx>[] = [];
  manager.outcomes.subscribe((o) => outcomes.push(o));
  return { manager, context, preview, unregister, outcomes };
}

function form(manager: CommandManager<Ctx>, id: string, tail?: string) {
  const r = manager.open(id, tail);
  if (r.kind !== 'form')
    throw new Error(`expected a form, got ${r.kind}${r.kind === 'blocked' ? `: ${r.reason}` : ''}`);
  return r.invocation;
}

/** A command whose run waits for `finish()`, honoring its signal or not. */
function slow(id: string, honorSignal: boolean, options: { form?: boolean } = { form: true }) {
  let finish: () => void = () => {};
  const command: CommandDefinition<Ctx> = {
    id,
    title: 'Load…',
    category: 'Data',
    ...(options.form ? { ...inline([{ id: 'entry', label: 'entry' }], 'Load') } : {}),
    defaults: () => ({ entry: '1abc' }),
    run: (v, c, signal) =>
      new Promise<CommandResult>((resolve, reject) => {
        finish = () => {
          c.log.push(`loaded ${String(v.entry)}`);
          resolve({ message: `loaded ${String(v.entry)}`, undoable: true });
        };
        if (honorSignal) signal.addEventListener('abort', () => reject(new Error('Load cancelled.')));
      }),
  };
  return { command, finish: () => finish() };
}

const focus: CommandDefinition<Ctx> = {
  id: 'view.focus',
  title: 'Focus on selection',
  category: 'View',
  shortcut: { key: 'f', alt: true },
  run: (_v, c) => {
    c.log.push('focus');
    return { message: 'focused' };
  },
};

describe('values', () => {
  const WHERE = [
    { value: 'here', label: 'in this scene' },
    { value: 'everywhere', label: 'everywhere' },
  ];
  const color: CommandDefinition<Ctx> = {
    id: 'look.color',
    title: 'Color…',
    category: 'Look',
    alias: 'color',
    ...inline([{ id: 'where', label: 'where', options: WHERE }], 'Color'),
    run: () => ({ message: 'colored' }),
  };

  it("typing an option's exact value and accepting takes that option, as the tail does", () => {
    const { manager } = make([color]);
    const inv = form(manager, 'look.color');
    inv.type('here');
    expect(inv.state.value.completions[0]?.value).toBe('here');
    inv.accept();
    expect(inv.state.value.values.where).toBe('here');
    inv.close();
    expect(form(manager, 'look.color', 'here').state.value.values.where).toBe('here');
  });

  it('a highlight the person moved is what they chose', () => {
    const { manager } = make([color]);
    const inv = form(manager, 'look.color');
    inv.type('here');
    inv.move(1);
    inv.accept();
    expect(inv.state.value.values.where).toBe('everywhere');
  });

  it('defaults that do not fit their params are dropped, not thrown from open', () => {
    const long: CommandDefinition<Ctx> = {
      id: 't.name',
      title: 'Name',
      category: 'Test',
      params: () => ({
        name: PD.Select('a', [
          ['a', 'A'],
          ['b', 'B'],
        ]),
      }),
      defaults: () => ({ name: 'not-an-option' }),
      run: () => ({}),
    };
    const fraction: CommandDefinition<Ctx> = {
      id: 't.count',
      title: 'Count',
      category: 'Test',
      params: () => ({ name: PD.Numeric(0, { min: 0, max: 10 }) }),
      defaults: () => ({ name: 20 }),
      run: () => ({}),
    };
    const { manager } = make([long, fraction]);
    expect(form(manager, 't.name').values().name).toBeUndefined();
    manager.close();
    expect(form(manager, 't.count').values().name).toBeUndefined();
  });

  it("a tail word a format reads as '' for a required field is not stored", () => {
    const formats: CommandFormats<Ctx> = {
      get: (name) =>
        name === 't.opt'
          ? {
              check: (v) => (v === 'none' ? { ok: true, value: '' } : { ok: false, message: 'no' }),
              complete: () => [],
            }
          : undefined,
    };
    const pick: CommandDefinition<Ctx> = {
      id: 't.pick',
      title: 'Pick',
      category: 'Test',
      alias: 'pick',
      ...inline([{ id: 'what', before: 'by', label: 'what', format: 't.opt' }], 'Go'),
      run: () => ({}),
    };
    const { manager } = make([pick], { formats });
    expect(form(manager, 't.pick', 'by none').values().what).toBeUndefined();
  });

  it('while the assistant resolves, edits do nothing', async () => {
    let answer: (v: unknown) => void = () => {};
    const show: CommandDefinition<Ctx> = {
      id: 'look.show',
      title: 'Show…',
      category: 'Look',
      ...inline([{ id: 'what', label: 'what' }], 'Show'),
      run: () => ({}),
      assist: { instruction: 'What to show.' },
    };
    const context = observableState<Ctx>({ current: 's1', selection: [], working: false, alignments: 0, log: [] });
    const manager = new CommandManager<Ctx>({
      context: context.asReadonly,
      formats: none,
      assistant: { propose: () => new Promise((r) => (answer = r as typeof answer)) },
    });
    manager.register([show]);
    const inv = form(manager, 'look.show');
    inv.type('the pocket');
    const resolving = inv.assist();
    expect(inv.state.value.status).toBe('resolving');
    expect(inv.set('what', 'typed meanwhile')).toBe(false);
    answer({ kind: 'values', values: { what: 'A:1-5' } });
    await resolving;
    expect(inv.values().what).toBe('A:1-5');
  });
});

describe('the open form', () => {
  const REPRS = ['cartoon', 'sticks'].map((r) => ({ value: r, label: r }));
  const show: CommandDefinition<Ctx> = {
    id: 'look.show',
    title: 'Show as…',
    category: 'Look',
    alias: 'show as',
    when: (c) => (c.working ? 'Waits for the agent' : true),
    ...inline([{ id: 'repr', before: 'as', label: 'representation', options: REPRS }], 'Apply'),
    preview: (v) => `${v.repr}`,
    run: (v, c) => {
      c.log.push(`show ${v.repr}`);
      return { message: 'shown' };
    },
  };

  it('blocked appears and clears with `when`; the preview and the run come back', async () => {
    const { manager, context, preview } = make([show]);
    const inv = form(manager, 'look.show', 'sticks');
    expect(preview.set).toHaveBeenLastCalledWith('sticks');
    context.set({ ...context.value, working: true });
    expect(inv.state.value.blocked).toBe('Waits for the agent');
    expect(preview.clear).toHaveBeenCalled();
    context.set({ ...context.value, working: false });
    expect(inv.state.value.blocked).toBeUndefined();
    expect(preview.set).toHaveBeenLastCalledWith('sticks');
    expect((await inv.confirm())?.message).toBe('shown');
  });

  const SELECTION = 'selection';
  const describeOp = (v: Record<string, unknown>, c: Ctx) =>
    `${c.current} ${v.target === SELECTION ? c.selection.join(',') : String(v.target)} ${String(v.color)}`;
  const color: CommandDefinition<Ctx> = {
    id: 'look.color',
    title: 'Color…',
    category: 'Look',
    alias: 'color',
    ...inline(
      [
        { id: 'target', label: 'what', options: [{ value: SELECTION, label: 'the selection' }] },
        { id: 'color', before: 'in', label: 'color', options: ['red', 'blue'].map((c) => ({ value: c, label: c })) },
      ],
      'Apply',
    ),
    defaults: (c) => (c.selection.length ? { target: SELECTION } : {}),
    preview: (v, c) => describeOp(v, c),
    run: (v, c) => {
      c.log.push(describeOp(v, c));
      return { message: 'colored' };
    },
  };

  it('the preview follows a selection and scene change, so it shows what Enter does', async () => {
    const { manager, context, preview } = make([color]);
    const inv = form(manager, 'look.color', 'in red');
    expect(preview.set).toHaveBeenLastCalledWith('s1 A:1 red');
    context.set({ ...context.value, current: 's2', selection: ['B:99'] });
    const shown = preview.set.mock.calls.at(-1)?.[0];
    await inv.enter();
    expect(context.value.log).toEqual(['s2 B:99 red']);
    expect(shown).toBe('s2 B:99 red');
  });

  it('with previewDeps, a context change outside them leaves the preview alone', () => {
    const { manager, context, preview } = make([color], { previewDeps: (c) => [c.current, c.selection.join()] });
    form(manager, 'look.color', 'in red');
    const drawn = preview.set.mock.calls.length;
    context.set({ ...context.value, alignments: 3 });
    context.set({ ...context.value, selection: ['A:1'] });
    expect(preview.set.mock.calls.length).toBe(drawn);
    context.set({ ...context.value, selection: ['A:2'] });
    expect(preview.set).toHaveBeenLastCalledWith('s1 A:2 red');
  });

  it('a throwing form or default is a blocked result, not an exception', () => {
    const broken: CommandDefinition<Ctx> = {
      id: 't.broken',
      title: 'Broken',
      category: 'Test',
      params: () => {
        throw new Error('No form here.');
      },
      run: () => ({}),
    };
    const { manager } = make([show, broken]);
    form(manager, 'look.show');
    expect(manager.open('t.broken')).toEqual({ kind: 'blocked', reason: 'No form here.' });
    expect(manager.active.value).toBeUndefined();
  });

  it('any param name is a field id', () => {
    const field = (fields: TestField[]): CommandDefinition<Ctx> => ({
      id: 't.cmd',
      title: 'Cmd',
      category: 'Test',
      ...inline(fields, 'Go'),
      run: () => ({}),
    });
    for (const fields of [
      [{ id: 'run', label: 'x' }],
      [{ id: 'form', label: 'x' }],
      [
        { id: 'a', before: 'by', label: 'a' },
        { id: 'a.before', label: 'b' },
        { id: 'before.a', label: 'c' },
      ],
    ]) {
      const { manager } = make([field(fields)]);
      expect(form(manager, 't.cmd').state.value.fields.map((f) => f.id)).toEqual(fields.map((f) => f.id));
    }
  });
});

describe('runs and outcomes', () => {
  async function running(manager: CommandManager<Ctx>, id: string) {
    const inv = form(manager, id);
    const done = inv.confirm();
    await tick();
    expect(inv.state.value.status).toBe('running');
    return { inv, done };
  }

  it('while a form runs, opening, keys and close are refused; the form stays', async () => {
    const load = slow('data.load', false);
    const { manager, context, outcomes, preview } = make([load.command, focus]);
    const { inv, done } = await running(manager, 'data.load');
    expect(manager.running.value?.command.id).toBe('data.load');
    expect(manager.open('view.focus')).toEqual({ kind: 'blocked', reason: 'Load is still running.' });
    expect(manager.handleKey({ key: 'ƒ', code: 'KeyF', altKey: true })).toBe(true);
    expect(outcomes.at(-1)).toMatchObject({ kind: 'blocked', message: 'Load is still running.' });
    expect(manager.close()).toBe(false);
    expect(manager.active.value).toBe(inv);
    load.finish();
    expect((await done)?.message).toBe('loaded 1abc');
    expect(context.value.log).toEqual(['loaded 1abc']);
    expect(outcomes.at(-1)).toMatchObject({ kind: 'done', result: { undoable: true } });
    expect(manager.running.value).toBeUndefined();
    expect(manager.active.value).toBeUndefined();
    expect(preview.clear).toHaveBeenCalled();
  });

  it('cancel aborts the signal and reports cancelled when the run gives up', async () => {
    const load = slow('data.load', true);
    const { manager, outcomes } = make([load.command]);
    const { inv, done } = await running(manager, 'data.load');
    expect(manager.cancel()).toBe(true);
    expect(inv.state.value.status).toBe('cancelling');
    expect(manager.running.value?.cancelling).toBe(true);
    expect(await done).toBeUndefined();
    expect(outcomes.map((o) => o.kind)).toEqual(['cancelled']);
    expect(manager.active.value).toBeUndefined();
    expect(manager.running.value).toBeUndefined();
  });

  it('a cancelled run that completes anyway reports done, with its Undo', async () => {
    const load = slow('data.load', false);
    const { manager, outcomes } = make([load.command]);
    const { done } = await running(manager, 'data.load');
    manager.cancel();
    load.finish();
    await done;
    expect(outcomes).toMatchObject([{ kind: 'done', result: { message: 'loaded 1abc', undoable: true } }]);
  });

  it('a run without a form can be cancelled too', async () => {
    const load = slow('data.reload', true, { form: false });
    const { manager, outcomes } = make([load.command]);
    const r = manager.open('data.reload');
    if (r.kind !== 'ran') throw new Error('expected a run');
    expect(manager.running.value?.command.id).toBe('data.reload');
    expect(manager.open('data.reload').kind).toBe('blocked');
    manager.cancel();
    expect(await r.result).toBeUndefined();
    expect(outcomes.map((o) => o.kind)).toEqual(['cancelled']);
    expect(manager.running.value).toBeUndefined();
  });

  it('a form whose run fails stays open with the error, and says so', async () => {
    const failing: CommandDefinition<Ctx> = {
      id: 't.fail',
      title: 'Fail…',
      category: 'Test',
      ...inline([{ id: 'x', label: 'x' }], 'Go'),
      defaults: () => ({ x: '1' }),
      run: () => {
        throw new Error('Nope.');
      },
    };
    const { manager, outcomes } = make([failing]);
    const inv = form(manager, 't.fail');
    await inv.confirm();
    expect(inv.state.value).toMatchObject({ status: 'editing', error: { message: 'Nope.' } });
    expect(outcomes).toMatchObject([{ kind: 'failed', message: 'Nope.', open: true }]);
    expect(manager.active.value).toBe(inv);
  });

  it('outcomes are events: a later subscriber gets none of the earlier ones', async () => {
    const { manager } = make([focus]);
    const r = manager.open('view.focus');
    if (r.kind !== 'ran') throw new Error('expected a run');
    await r.result;
    const later: unknown[] = [];
    manager.outcomes.subscribe((o) => later.push(o)).unsubscribe();
    expect(later).toEqual([]);
  });

  it('leaving a screen cancels its run; an outcome that settles later reaches no one, not the next screen', async () => {
    const load = slow('data.load', false);
    const { manager, unregister, outcomes } = make([load.command]);
    const { done } = await running(manager, 'data.load');
    unregister();
    expect(manager.active.value).toBeUndefined();
    expect(manager.running.value).toBeUndefined();
    // The next screen registers the same definitions and can run at once.
    const again = manager.register([load.command]);
    expect(manager.open('data.load').kind).toBe('form');
    load.finish();
    await done;
    expect(outcomes).toEqual([]);
    again();
  });

  it('a run without a form that settles after its screen left reports to no one', async () => {
    const load = slow('data.reload', false, { form: false });
    const { manager, unregister, outcomes } = make([load.command]);
    const r = manager.open('data.reload');
    if (r.kind !== 'ran') throw new Error('expected a run');
    unregister();
    manager.register([load.command]);
    load.finish();
    await r.result;
    expect(outcomes).toEqual([]);
  });
});

describe('matching', () => {
  it('a shared alias resolves to the command that can run, not the first registered', () => {
    const move: CommandDefinition<Ctx> = {
      id: 'look.move',
      title: 'Move by alignment…',
      category: 'Look',
      alias: ['move by alignment', 'superpose'],
      when: (c) => (c.alignments ? true : 'No alignment yet (Align… makes one).'),
      ...inline([{ id: 'by', label: 'alignment' }], 'Move'),
      run: () => ({}),
    };
    const align: CommandDefinition<Ctx> = {
      id: 'context.align',
      title: 'Align…',
      category: 'Context',
      alias: ['align', 'superpose'],
      ...inline(
        [
          { id: 'moving', label: 'structure' },
          { id: 'fixed', before: 'onto', label: 'structure' },
        ],
        'Align',
      ),
      run: () => ({}),
    };
    const { manager, context } = make([move, align]);
    expect(manager.query('superpose 1abc onto 2xyz').direct?.command.id).toBe('context.align');
    context.set({ ...context.value, alignments: 1 });
    expect(manager.query('superpose 1abc onto 2xyz').direct?.command.id).toBe('look.move');
  });

  it('a blocked direct match is listed first, then what the whole text matches', () => {
    const addRepr: CommandDefinition<Ctx> = {
      id: 'look.add',
      title: 'Add representation…',
      category: 'Look',
      alias: 'add',
      when: () => 'No scene in view.',
      ...inline([{ id: 'what', label: 'what' }], 'Add'),
      run: () => ({}),
    };
    const other: CommandDefinition<Ctx> = {
      id: 't.other',
      title: 'Load and add sticks',
      category: 'Test',
      run: () => ({}),
    };
    const newScene: CommandDefinition<Ctx> = {
      id: 'scene.new',
      title: 'New scene…',
      category: 'Scene',
      alias: ['add scene', 'new'],
      ...inline([{ id: 'title', label: 'title' }], 'Add'),
      run: () => ({}),
    };
    const { manager } = make([addRepr, other, newScene]);
    // The whole of another alias is that command, not "add" with an argument.
    expect(manager.query('add scene').direct).toBeUndefined();
    expect(manager.query('add scene').matches[0]?.command.id).toBe('scene.new');
    expect(manager.query('add scene Heme').direct?.command.id).toBe('scene.new');
    const q = manager.query('add sticks');
    expect(q.direct).toMatchObject({ command: { id: 'look.add' }, blocked: 'No scene in view.', tail: 'sticks' });
    expect(q.matches.map((m) => m.command.id)).toEqual(['look.add', 't.other']);
  });

  it("among commands with one chord, one that can run wins; else the first one's reason", () => {
    const blocked: CommandDefinition<Ctx> = { ...focus, id: 'view.blocked', when: () => 'Not now.' };
    const { manager, context, outcomes } = make([blocked, focus]);
    expect(manager.handleKey({ key: 'f', altKey: true }, false)).toBe(true);
    expect(context.value.log).toEqual(['focus']);
    const only = make([blocked]);
    expect(only.manager.handleKey({ key: 'f', altKey: true }, false)).toBe(true);
    expect(only.outcomes).toMatchObject([{ kind: 'blocked', message: 'Not now.' }]);
    expect(outcomes.filter((o) => o.kind === 'blocked')).toEqual([]);
  });
});
