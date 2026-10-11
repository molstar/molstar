/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '../../reactive/observable-state.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { formFields } from '../../forms/flat-fields.js';
import { CommandManager } from '../manager.js';
import type { CommandAssistant, CommandDefinition, CommandFormats } from '../types.js';
import { hintsOf, inline, paramsOf, type TestField } from './inline.js';

interface Ctx {
  selection: string[];
  working: boolean;
  log: string[];
}

const RANGE = /^[A-Z]:\d+(-\d+)?$/;
const TEXT = { check: (v: string) => ({ ok: true as const, value: v }), complete: () => [] };
const formats: CommandFormats<Ctx> = {
  get: (name) =>
    name === 't.text'
      ? TEXT
      : name === 't.residues'
        ? {
            check: (v) =>
              v
                .split(',')
                .map((x) => x.trim())
                .every((x) => RANGE.test(x))
                ? {
                    ok: true,
                    value: v
                      .split(',')
                      .map((x) => x.trim())
                      .join(', '),
                  }
                : { ok: false, message: `"${v}" is not residues.`, unreadable: !v.includes(':') },
            complete: (text) =>
              [{ value: 'A:248-256', label: 'P-loop' }].filter((c) =>
                c.label.toLowerCase().includes(text.toLowerCase()),
              ),
          }
        : undefined,
};

const REPRS = ['cartoon', 'sticks', 'spacefill'].map((r) => ({ value: r, label: r }));

function setup(options: { assistant?: CommandAssistant<Ctx> } = {}) {
  const context = observableState<Ctx>({ selection: ['A:315'], working: false, log: [] });
  const preview = { set: jest.fn(), clear: jest.fn() };
  const show: CommandDefinition<Ctx> = {
    id: 'look.show',
    title: 'Show as…',
    category: 'Look',
    alias: 'show as',
    when: (c) => (c.working ? 'Waits for the agent' : true),
    ...inline(
      [
        { id: 'target', label: 'what', format: 't.residues' },
        { id: 'repr', before: 'as', label: 'representation', options: REPRS },
        { id: 'color', before: 'in', label: 'colour', optional: true },
      ],
      'Apply',
    ),
    defaults: (c) => (c.selection.length ? { target: c.selection.join(', ') } : {}),
    preview: (v) => `${v.target} ${v.repr}${v.color ? ` ${v.color}` : ''}`,
    run: (v, c) => {
      c.log.push(`show ${v.target} ${v.repr} ${v.color ?? '-'}`);
      return { message: 'shown', undoable: true };
    },
    assist: { instruction: 'What to show and how.' },
  };
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
  const select: CommandDefinition<Ctx> = {
    id: 'view.select',
    title: 'Select residues…',
    category: 'View',
    alias: 'select',
    ...inline([{ id: 'residues', label: 'residues', format: 't.residues' }], 'Select'),
    run: (v, c) => {
      c.log.push(`select ${v.residues}`);
      return { message: 'selected' };
    },
    assist: { instruction: 'Residues to select.' },
  };
  const manager = new CommandManager<Ctx>({
    context: context.asReadonly,
    formats,
    preview,
    ...(options.assistant ? { assistant: options.assistant } : {}),
  });
  manager.register([show, focus, select]);
  return { manager, context, preview };
}

describe('inline fields', () => {
  it('reads params back as fields in order, with before-words, options, formats and what is required', () => {
    const spec: TestField[] = [
      { id: 'target', label: 'what', format: 't.residues' },
      { id: 'repr', before: 'as', label: 'representation', options: REPRS },
      { id: 'opacity', label: 'opacity', slider: { min: 0, max: 1, step: 0.05 }, optional: true },
    ];
    const fields = formFields(paramsOf(spec), hintsOf(spec));
    expect(fields.map((f) => [f.id, f.kind, f.required, f.before])).toEqual([
      ['target', 'text', true, undefined],
      ['repr', 'select', true, 'as'],
      ['opacity', 'slider', false, undefined],
    ]);
    expect(fields[0]?.format).toBe('t.residues');
    expect(fields[1]?.options?.map((o) => o.value)).toEqual(['cartoon', 'sticks', 'spacefill']);
  });
});

describe('CommandManager.query', () => {
  it('reports every successful run even when a command returns the same result object', async () => {
    const { manager } = setup();
    const result = { message: 'focused' };
    manager.register([{ id: 'view.focus', title: 'Focus', category: 'View', run: () => result }]);
    const outcomes = jest.fn();
    const subscription = manager.outcomes.subscribe((outcome) => {
      if (outcome.kind === 'done') outcomes(outcome);
    });
    for (let i = 0; i < 2; i++) {
      const opened = manager.open('view.focus');
      if (opened.kind !== 'ran') throw new Error('expected immediate run');
      await opened.result;
    }
    expect(outcomes).toHaveBeenCalledTimes(2);
    subscription.unsubscribe();
  });

  it('reports every blocked shortcut even when the reason is unchanged', () => {
    const { manager } = setup();
    manager.register([
      {
        id: 'view.focus',
        title: 'Focus',
        category: 'View',
        shortcut: { key: 'f', alt: true },
        when: () => 'Select something first',
        run: () => ({}),
      },
    ]);
    const errors = jest.fn();
    const subscription = manager.outcomes.subscribe((outcome) => {
      if (outcome.kind === 'blocked') errors(outcome);
    });
    expect(manager.handleKey({ key: 'f', altKey: true }, false)).toBe(true);
    expect(manager.handleKey({ key: 'f', altKey: true }, false)).toBe(true);
    expect(errors).toHaveBeenCalledTimes(2);
    subscription.unsubscribe();
  });

  it('lists recently used first, then matches best first', async () => {
    const { manager } = setup();
    await (manager.open('view.focus') as { result: Promise<unknown> }).result;
    const q = manager.query('');
    expect(q.recent.map((m) => m.command.id)).toEqual(['view.focus']);
    expect(q.matches.map((m) => m.command.id)).not.toContain('view.focus');
    expect(manager.query('sel').matches[0]?.command.id).toBe('view.select');
  });

  it("reads '<alias> <arguments>' as that command", () => {
    const { manager } = setup();
    const q = manager.query('show as sticks red');
    expect(q.direct?.command.id).toBe('look.show');
    expect(q.direct?.tail).toBe('sticks red');
  });

  it('keeps the start of a title a query; the longest alias that fits takes the arguments', () => {
    const { manager } = setup();
    manager.register([
      {
        id: 'view.focus',
        title: 'Focus on…',
        category: 'View',
        alias: ['focus on', 'focus'],
        ...inline([{ id: 'target', label: 'what' }], 'Focus'),
        run: () => ({}),
      },
    ]);
    expect(manager.query('focus on').direct).toBeUndefined();
    expect(manager.query('focus on').matches[0]?.command.id).toBe('view.focus');
    expect(manager.query('focus on A 316').direct).toMatchObject({ tail: 'A 316' });
    expect(manager.query('focus A 316').direct).toMatchObject({ tail: 'A 316' });
    expect(manager.query('show as ').direct).toBeUndefined();
  });

  it('keeps commands that cannot run, with the reason', () => {
    const { manager, context } = setup();
    context.set({ ...context.asReadonly.value, working: true });
    const m = manager.query('show').matches.find((x) => x.command.id === 'look.show');
    expect(m?.blocked).toBe('Waits for the agent');
    expect(manager.open('look.show')).toEqual({ kind: 'blocked', reason: 'Waits for the agent' });
  });
});

describe('CommandInvocation', () => {
  it('starts from defaults and the typed tail, previews once complete, runs on Enter', async () => {
    const { manager, context, preview } = setup();
    const r = manager.open('look.show', 'sticks');
    if (r.kind !== 'form') throw new Error('expected a form');
    const inv = r.invocation;
    expect(inv.state.value.values).toMatchObject({ target: 'A:315', repr: 'sticks' });
    expect(inv.state.value.missing).toEqual([]);
    expect(preview.set).toHaveBeenLastCalledWith('A:315 sticks');
    const outcomes: string[] = [];
    const sub = manager.outcomes.subscribe((o) => outcomes.push(o.kind === 'done' ? `${o.result.message}` : o.kind));
    await inv.enter();
    sub.unsubscribe();
    expect(context.asReadonly.value.log).toEqual(['show A:315 sticks -']);
    expect(outcomes).toEqual(['shown']);
    expect(manager.active.value).toBeUndefined();
    expect(preview.clear).toHaveBeenCalled();
    expect(manager.recents.value[0]).toBe('look.show');
  });

  it("reads a field's value after its before-word, up to the next one", () => {
    const { manager } = setup();
    const r = manager.open('look.show', 'A:248-256, A:300 as spacefill in dark red');
    if (r.kind !== 'form') throw new Error('expected a form');
    expect(r.invocation.state.value.values).toMatchObject({
      target: 'A:248-256, A:300',
      repr: 'spacefill',
      color: 'dark red',
    });
    // A phrase that does not read is left for the rest: "as" before no representation.
    const loose = manager.open('look.show', 'A:300 as nothing');
    if (loose.kind !== 'form') throw new Error('expected a form');
    expect(loose.invocation.state.value.values.repr).toBeUndefined();
  });

  it('completes the active field, previews the highlighted completion, and moves on when one is taken', () => {
    const { manager, preview } = setup();
    const r = manager.open('look.show');
    if (r.kind !== 'form') throw new Error('expected a form');
    const inv = r.invocation;
    expect(inv.state.value.field).toBe('repr');
    inv.type('sp');
    expect(inv.state.value.completions.map((c) => c.value)).toEqual(['spacefill']);
    expect(preview.set).toHaveBeenLastCalledWith('A:315 spacefill');
    inv.accept();
    expect(inv.state.value.values.repr).toBe('spacefill');
    expect(inv.state.value.field).toBe('color');
    inv.close();
    expect(preview.clear).toHaveBeenCalled();
  });

  it("checks typed text with the field's format, and offers the assistant only for words", () => {
    const assistant: CommandAssistant<Ctx> = { propose: jest.fn() };
    const { manager } = setup({ assistant });
    const r = manager.open('view.select');
    if (r.kind !== 'form') throw new Error('expected a form');
    const inv = r.invocation;
    inv.type('A:900x');
    expect(inv.state.value.canAssist).toBe(false);
    inv.accept();
    expect(inv.state.value.error?.message).toBe('"A:900x" is not residues.');
    inv.type('the gatekeeper');
    expect(inv.state.value.canAssist).toBe(true);
    inv.type('A:100-120');
    inv.accept();
    expect(inv.state.value.values.residues).toBe('A:100-120');
  });

  it("takes an assistant's proposal as proposed values; nothing runs until confirmed", async () => {
    const assistant: CommandAssistant<Ctx> = {
      propose: jest.fn(async () => ({
        kind: 'values' as const,
        values: { residues: 'A:248-256,A:315' },
        note: 'P-loop and T315',
      })),
    };
    const { manager, context } = setup({ assistant });
    const r = manager.open('view.select', 'the P-loop and the gatekeeper');
    if (r.kind !== 'form') throw new Error('expected a form');
    const inv = r.invocation;
    expect(inv.state.value.text).toBe('the P-loop and the gatekeeper');
    expect(inv.state.value.canAssist).toBe(true);
    await inv.tab();
    expect(assistant.propose).toHaveBeenCalledWith(
      expect.objectContaining({ text: 'the P-loop and the gatekeeper', field: 'residues' }),
      expect.any(AbortSignal),
    );
    expect(inv.state.value.values.residues).toBe('A:248-256, A:315');
    expect(inv.state.value.proposed).toEqual(['residues']);
    expect(inv.state.value.note).toBe('P-loop and T315');
    expect(context.asReadonly.value.log).toEqual([]);
    await inv.enter();
    expect(context.asReadonly.value.log).toEqual(['select A:248-256, A:315']);
  });

  it('hands words it cannot express to the agent', async () => {
    const assistant: CommandAssistant<Ctx> = {
      propose: async () => ({ kind: 'escalate', reason: 'That needs two commands.' }),
    };
    const { manager } = setup({ assistant });
    const r = manager.open('view.select', 'the pocket, and colour it red');
    if (r.kind !== 'form') throw new Error('expected a form');
    await r.invocation.assist();
    expect(r.invocation.state.value.escalated).toBe('the pocket, and colour it red');
    expect(r.invocation.state.value.error?.message).toBe('That needs two commands.');
  });

  it('blocks the run when the command cannot run any more', async () => {
    const { manager, context } = setup();
    const r = manager.open('look.show', 'sticks');
    if (r.kind !== 'form') throw new Error('expected a form');
    context.set({ ...context.asReadonly.value, working: true });
    expect(r.invocation.state.value.blocked).toBe('Waits for the agent');
    await r.invocation.enter();
    expect(context.asReadonly.value.log).toEqual([]);
    expect(r.invocation.state.value.error?.message).toBe('Waits for the agent');
  });

  it('goes back a field on backspace, clearing values first', () => {
    const { manager } = setup();
    const r = manager.open('look.show', 'sticks');
    if (r.kind !== 'form') throw new Error('expected a form');
    const inv = r.invocation;
    inv.focus('repr');
    expect(inv.back()).toBe(true);
    expect(inv.state.value.values.repr).toBeUndefined();
    expect(inv.back()).toBe(true);
    expect(inv.state.value.field).toBe('target');
  });
});

describe('keybindings', () => {
  it('runs a command from its chord (by key code, so Alt on macOS works)', async () => {
    const { manager, context } = setup();
    expect(manager.handleKey({ key: 'ƒ', code: 'KeyF', altKey: true })).toBe(true);
    await Promise.resolve();
    expect(context.asReadonly.value.log).toEqual(['focus']);
    expect(manager.handleKey({ key: 'f' })).toBe(false);
  });
});

describe('forms that follow their values', () => {
  const SOURCES = [
    { value: 'pdb', label: 'PDB entry' },
    { value: 'ihm', label: 'Integrative (PDB-IHM)' },
  ];
  const loadFields = (ihm: boolean): TestField[] => [
    { id: 'entry', label: ihm ? 'PDB-IHM id' : 'id or file', format: 't.text' },
    ...(ihm ? [{ id: 'model', before: 'model', label: 'model', optional: true }] : []),
    { id: 'source', before: 'from', label: 'source', options: SOURCES, optional: true },
  ];
  function setupLoad() {
    const context = observableState<Ctx>({ selection: [], working: false, log: [] });
    const load: CommandDefinition<Ctx> = {
      id: 'load',
      title: 'Load structure…',
      category: 'Context',
      alias: 'load',
      params: (_c, v) => paramsOf(loadFields(v.source === 'ihm')),
      fields: hintsOf(loadFields(true)),
      action: 'Load',
      reshapes: ['source'],
      presets: () =>
        SOURCES.map((s) => ({
          id: s.value,
          title: `from ${s.label}`,
          values: { source: s.value },
          keywords: [s.value],
        })),
      run: (v, c) => {
        c.log.push(`load ${v.entry} ${v.source ?? '-'} ${v.model ?? '-'}`);
        return {};
      },
    };
    const manager = new CommandManager<Ctx>({ context: context.asReadonly, formats });
    manager.register([load]);
    return { manager, context };
  }
  const ids = (fields: readonly { id: string }[]) => fields.map((f) => f.id);

  it('builds the form again when a field that reshapes it is set or cleared, keeping the values it still has', async () => {
    const { manager, context } = setupLoad();
    const r = manager.open('load', '8ZZE');
    if (r.kind !== 'form') throw new Error('expected a form');
    const inv = r.invocation;
    expect(ids(inv.state.value.fields)).toEqual(['entry', 'source']);
    inv.set('source', 'ihm');
    expect(ids(inv.state.value.fields)).toEqual(['entry', 'model', 'source']);
    expect(inv.state.value.fields[0]?.label).toBe('PDB-IHM id');
    expect(inv.state.value.values).toEqual({ entry: '8ZZE', source: 'ihm' });
    expect(inv.state.value.field).toBe('model');
    inv.set('model', '2');
    inv.clear('source');
    expect(ids(inv.state.value.fields)).toEqual(['entry', 'source']);
    // Cleared, it stays on the field (Backspace), as every other field does.
    expect(inv.state.value.field).toBe('source');
    expect(inv.state.value.values).toEqual({ entry: '8ZZE' });
    inv.set('source', 'pdb');
    await inv.confirm();
    expect(context.asReadonly.value.log).toEqual(['load 8ZZE pdb -']);
  });

  it('a form that fails to build for a value keeps the one it had, and says why', () => {
    const { manager } = setupLoad();
    const load = manager.get('load');
    if (!load) throw new Error('no load');
    const params = load.params;
    load.params = (c, v) => {
      if (v.source === 'ihm') throw new Error('PDB-IHM is not available here.');
      return params?.(c, v) as ReturnType<NonNullable<typeof params>>;
    };
    const r = manager.open('load', '8ZZE');
    if (r.kind !== 'form') throw new Error('expected a form');
    r.invocation.set('source', 'ihm');
    expect(ids(r.invocation.state.value.fields)).toEqual(['entry', 'source']);
    expect(r.invocation.state.value.error?.message).toBe('PDB-IHM is not available here.');
  });

  it('reshapes from the typed tail', () => {
    const { manager } = setupLoad();
    const r = manager.open('load', '8ZZE from ihm');
    if (r.kind !== 'form') throw new Error('expected a form');
    expect(ids(r.invocation.state.value.fields)).toEqual(['entry', 'model', 'source']);
    expect(r.invocation.state.value.values).toEqual({ entry: '8ZZE', source: 'ihm' });
  });

  it('finds presets by their words, never lists them unasked, and opens with their values', () => {
    const { manager } = setupLoad();
    expect(manager.query('').matches.map((m) => m.preset?.id)).toEqual([undefined]);
    // The command's own words find it alone; words of a preset find that preset.
    expect(manager.query('load').matches.map((m) => m.preset?.id ?? m.command.id)).toEqual(['load']);
    const ihm = manager.query('pdb-ihm').matches;
    expect(ihm[0]?.preset?.id).toBe('ihm');
    expect(manager.query('structure integrative').matches.map((m) => m.preset?.id)).toEqual(['ihm']);
    const r = manager.open('load', undefined, 'ihm');
    if (r.kind !== 'form') throw new Error('expected a form');
    expect(ids(r.invocation.state.value.fields)).toEqual(['entry', 'model', 'source']);
    expect(r.invocation.state.value.values).toEqual({ source: 'ihm' });
    expect(manager.open('load', undefined, 'nope')).toEqual({
      kind: 'blocked',
      reason: 'Load structure… has no "nope".',
    });
  });
});

describe('booleans and numbers', () => {
  it('reads on/off into a boolean field and passes form-only params with their defaults', async () => {
    const context = observableState<Ctx>({ selection: [], working: false, log: [] });
    const around: CommandDefinition<Ctx> = {
      id: 'select.around',
      title: 'Select around…',
      category: 'Select',
      alias: 'around',
      params: () => ({
        radius: PD.Numeric(5, { min: 0, max: 20 }),
        whole: PD.Boolean(true, { label: 'whole residues' }),
        quality: PD.Group({ level: PD.Text('auto') }),
      }),
      fields: { whole: { before: 'whole' } },
      defaults: () => ({ radius: 5 }),
      run: (v, c) => {
        c.log.push(`${v.radius} ${v.whole} ${(v.quality as { level: string }).level}`);
        return {};
      },
    };
    const manager = new CommandManager<Ctx>({ context: context.asReadonly, formats });
    manager.register([around]);
    const r = manager.open('select.around', 'whole off');
    if (r.kind !== 'form') throw new Error('expected a form');
    expect(r.invocation.state.value.values).toEqual({ radius: 5, whole: false });
    expect(r.invocation.state.value.formOnly).toEqual(['quality']);
    await r.invocation.confirm();
    expect(context.value.log).toEqual(['5 false auto']);
  });
});
