/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { KeyPlatform } from '../../commands/keys.js';
import { CommandManager } from '../../commands/manager.js';
import { flush, render } from '../../dom/index.js';
import { shownInOrder } from '../../testing/virtual-list.js';
import { axeViolations } from '../../testing/axe.js';
import { historyCommands } from '../history-commands.js';
import { HistoryModel, type HistoryPort, type HistoryStep, type UndoResult } from '../history-model.js';
import { createHistoryRows, HistoryView, type HistoryViewProps } from '../history-view.js';

// History's view (history §4): the steps oldest first, the position, Undo and Redo, to a step, unsafe, blocked.

const mac = KeyPlatform.mac;
let dispose: (() => void) | undefined;
afterEach(() => {
  dispose?.();
  dispose = undefined;
  KeyPlatform.mac = mac;
  document.body.innerHTML = '';
});

function store() {
  let a = 0;
  const port: HistoryPort<number> = {
    checkpoint: () => a,
    restore: async (c) => {
      a = c;
    },
  };
  return { port, get: () => a, set: (v: number) => (a = v) };
}

const person = { kind: 'person' };
const settle = async () => {
  await new Promise((r) => setTimeout(r, 0));
  flush();
};

async function mount(n: number, onResult = vi.fn(), props: Partial<HistoryViewProps> = {}) {
  KeyPlatform.mac = false;
  const s = store();
  const history = new HistoryModel(s.port, { now: () => Date.UTC(2026, 0, 1, 10, 0, 0) });
  for (let i = 1; i <= n; i++) await history.run(`Set a to ${i}`, person, async () => void s.set(i));
  const el = document.createElement('div');
  document.body.append(el);
  dispose = render(() => <HistoryView history={history} onResult={onResult} {...props} />, el);
  flush();
  return { s, history, el, onResult };
}

/** A step as the model makes them (its runs read through `lastRuns`). */
function step(id: number, runs: string[], dropped = 0): HistoryStep {
  return {
    id,
    label: `s${id}`,
    origin: person,
    time: 0,
    runs,
    runCount: runs.length + dropped,
    droppedRuns: dropped,
    lastRuns: (n) => runs.slice(Math.max(0, runs.length - n)),
  };
}

const act = (p: Promise<UndoResult>, to: (r: UndoResult) => void) => void p.then(to);
const rows = (el: HTMLElement) => shownInOrder(el, '[data-state]');
const labels = (el: HTMLElement) => rows(el).map((r) => r.querySelector('[data-tip-part]')?.textContent);
const states = (el: HTMLElement) => rows(el).map((r) => r.getAttribute('data-state'));
const button = (el: HTMLElement, name: string) =>
  [...el.querySelectorAll<HTMLButtonElement>('button')].find(
    (b) => b.getAttribute('aria-label') === name || b.getAttribute('data-tip') === name,
  )!;

describe('createHistoryRows', () => {
  it('keeps the rows of an unchanged step, and a changed step’s unchanged rows (history §4)', () => {
    const list = createHistoryRows({ runs: 3 });
    const older = step(1, ['x', 'y']);
    const turn = step(2, ['a', 'b']);
    const before = list([turn, older], []);
    // The turn takes one more run: a new step object, as the model publishes it.
    const after = list([step(2, ['a', 'b', 'c']), older], []);
    const byKey = (rows: readonly { key: string }[]) => new Map(rows.map((r) => [r.key, r]));
    const was = byKey(before);
    const now = byKey(after);
    for (const key of ['start', 's1', 'r1:0', 'r1:1', 'r2:0', 'r2:1']) expect(now.get(key)).toBe(was.get(key));
    expect(now.get('s2')).not.toBe(was.get('s2'));
    expect(now.get('r2:2')).toMatchObject({ label: 'c' });
  });

  it('marks the moves past an unsafe step unreachable', () => {
    const unsafe = { ...step(1, ['x']), unsafe: 'Too old', unsafeCode: 'too-old' as const };
    const r = createHistoryRows()([step(3, ['c']), step(2, ['b']), unsafe], []);
    expect(r.map((x) => ('reachable' in x ? x.reachable : undefined))).toEqual([false, true, true, true]);
  });

  it('lists the start, the steps oldest first with a group’s runs, then the undone ones', () => {
    const r = createHistoryRows()([step(2, ['a', 'b']), step(1, ['s1'])], [step(3, ['s3'])]);
    expect(r.map((x) => x.key)).toEqual(['start', 's1', 's2', 'r2:0', 'r2:1', 's3']);
    expect(r.map((x) => ('current' in x ? x.current : undefined))).toEqual([
      false,
      false,
      true,
      undefined,
      undefined,
      false,
    ]);
  });

  it('lists a group’s latest runs, a row for the earlier ones, all of them kept when expanded', () => {
    const runs = Array.from({ length: 30 }, (_, i) => `run ${i}`);
    const folded = createHistoryRows({ runs: 3 })([step(1, runs, 5)], []);
    expect(folded.map((x) => x.key)).toEqual(['start', 's1', 'm1', 'r1:32', 'r1:33', 'r1:34']);
    expect(folded[2]).toMatchObject({ kind: 'more', earlier: 32, expandable: true, expanded: false });
    const open = createHistoryRows({ runs: 3 })([step(1, runs, 5)], [], new Set([1]));
    expect(open).toHaveLength(3 + 30);
    expect(open[2]).toMatchObject({ kind: 'more', earlier: 5, expanded: true });
    expect(open[3]).toMatchObject({ key: 'r1:5', label: 'run 0' });
    // Few runs: no row of more.
    expect(createHistoryRows({ runs: 3 })([step(2, ['a', 'b'])], []).map((x) => x.kind)).toEqual([
      'start',
      'step',
      'run',
      'run',
    ]);
  });
});

describe('HistoryView', () => {
  it('shows the steps oldest first, the latest marked, the origin and time as the description', async () => {
    const { el } = await mount(3);
    expect(labels(el)).toEqual(['Start', 'Set a to 1', 'Set a to 2', 'Set a to 3']);
    expect(states(el)).toEqual(['done', 'done', 'done', 'current']);
    expect(rows(el)[3].getAttribute('aria-current')).toBe('step');
    expect(rows(el)[1].textContent).toContain('You · ');
    expect(el.textContent).toContain('3 steps');
  });

  it('undoes and redoes from its toolbar; the undone steps muted after the position', async () => {
    const { el, s, history } = await mount(2);
    const undo = button(el, 'Undo');
    expect(undo.getAttribute('data-tip')).toBe('Undo Set a to 2 (Ctrl+Z)');
    expect(undo.getAttribute('aria-keyshortcuts')).toBe('Control+Z');
    undo.click();
    await settle();
    expect(s.get()).toBe(1);
    expect(states(el)).toEqual(['done', 'current', 'undone']);
    const redo = button(el, 'Redo');
    expect(redo.disabled).toBe(false);
    expect(redo.getAttribute('data-tip')).toBe('Redo Set a to 2 (Ctrl+Shift+Z)');
    redo.click();
    await settle();
    expect(s.get()).toBe(2);
    expect(history.redoable.value).toEqual([]);
    expect(states(el)).toEqual(['done', 'done', 'current']);
  });

  it('undoes to a step and redoes to one undone, from the rows', async () => {
    const { el, s } = await mount(4);
    // "Undo to <step>" on the first step: it stays, the three after it are undone.
    (rows(el)[1].querySelector('button[data-tip="Undo to Set a to 1"]') as HTMLButtonElement).click();
    await settle();
    expect(s.get()).toBe(1);
    expect(states(el)).toEqual(['done', 'current', 'undone', 'undone', 'undone']);
    expect(rows(el)[3].querySelector('button')?.getAttribute('aria-label')).toBe('Redo to Set a to 3');
    (rows(el)[3].querySelector('button[data-tip="Redo to Set a to 3"]') as HTMLButtonElement).click();
    await settle();
    expect(s.get()).toBe(3);
    expect(states(el)).toEqual(['done', 'done', 'done', 'current', 'undone']);
    // The start undoes everything.
    (rows(el)[0].querySelector('button[data-tip="Undo every step"]') as HTMLButtonElement).click();
    await settle();
    expect(s.get()).toBe(0);
    expect(states(el)[0]).toBe('current');
  });

  it('warns on an unsafe step and says why Undo is blocked', async () => {
    const { el, history, onResult } = await mount(1);
    history.invalidate('Changed outside the history');
    flush();
    expect(rows(el)[1].hasAttribute('data-unsafe')).toBe(true);
    expect(rows(el)[1].querySelector('[data-tip="Changed outside the history"]')).not.toBeNull();
    expect(el.querySelector('[role=note]')?.textContent).toBe('Changed outside the history');
    const undo = button(el, 'Undo');
    expect(undo.disabled).toBe(true);
    expect(undo.getAttribute('data-tip')).toBe('Changed outside the history');
    // "Undo every step" cannot reach the start: greyed, and the model refuses it up front.
    expect((rows(el)[0].querySelector('button') as HTMLButtonElement).disabled).toBe(true);
    act(history.moveTo(undefined), onResult);
    await settle();
    expect(onResult).toHaveBeenCalledWith({ ok: false, code: 'unsafe', reason: 'Changed outside the history' });
  });

  it('greys the moves that cannot be reached: past an unsafe step (history §2)', async () => {
    const s = store();
    const history = new HistoryModel(s.port, { checkpoints: 2 });
    for (let i = 1; i <= 4; i++) await history.run(`Set a to ${i}`, person, async () => void s.set(i));
    const el = document.createElement('div');
    document.body.append(el);
    dispose = render(() => <HistoryView history={history} />, el);
    flush();
    // Steps 1 and 2 are too old: the start and "Undo to Set a to 1" (which undoes step 2) cannot run; to step 2 can.
    const action = (i: number) => rows(el)[i].querySelector('button') as HTMLButtonElement;
    expect(action(0).disabled).toBe(true);
    expect(action(1).disabled).toBe(true);
    expect(action(2).disabled).toBe(false);
    expect(action(3).disabled).toBe(false);
    action(2).click();
    await settle();
    expect(s.get()).toBe(2);
  });

  it('says nothing while a change runs: the buttons are greyed', async () => {
    const { el, history } = await mount(1);
    let open!: () => void;
    const run = history.run('Slow', person, () => new Promise<void>((r) => (open = r)));
    flush();
    expect(history.undoBlocked.value?.code).toBe('running');
    expect(el.querySelector('[role=note]')).toBeNull();
    expect(button(el, 'Undo').disabled).toBe(true);
    await settle();
    open();
    await run;
    flush();
    expect(el.querySelector('[role=note]')).toBeNull();
  });

  it('says the keys of the app’s commands, a remap too; by platform without them', async () => {
    const plain = await mount(1);
    expect(button(plain.el, 'Undo').getAttribute('aria-keyshortcuts')).toBe('Control+Z');
    dispose?.();
    document.body.innerHTML = '';
    KeyPlatform.mac = true;
    const s = store();
    const history = new HistoryModel(s.port);
    await history.run('Set a to 1', person, async () => void s.set(1));
    const commands = new CommandManager<unknown>({
      context: new State<unknown>(undefined),
      formats: { get: () => undefined },
    });
    commands.register(historyCommands(history));
    const el = document.createElement('div');
    document.body.append(el);
    dispose = render(() => <HistoryView history={history} commands={commands} />, el);
    flush();
    const undo = button(el, 'Undo');
    expect(undo.getAttribute('aria-keyshortcuts')).toBe('Meta+Z');
    expect(undo.getAttribute('data-tip')).toBe('Undo Set a to 1 (⌘Z)');
    commands.keymap.remap('edit.undo', { key: 'u', alt: true });
    flush();
    expect(undo.getAttribute('aria-keyshortcuts')).toBe('Alt+U');
    expect(undo.getAttribute('data-tip')).toBe('Undo Set a to 1 (⌥U)');
    commands.keymap.remap('edit.undo', null);
    flush();
    expect(undo.hasAttribute('aria-keyshortcuts')).toBe(false);
    expect(undo.getAttribute('data-tip')).toBe('Undo Set a to 1');
  });

  it('folds a long turn’s runs and opens them on demand', async () => {
    const s = store();
    const history = new HistoryModel(s.port);
    const agent = { kind: 'agent', turn: 't' };
    const end = history.beginGroup('Turn', agent);
    for (let i = 0; i < 15; i++) await history.run(`call ${i}`, agent, async () => void s.set(i));
    end();
    const el = document.createElement('div');
    document.body.append(el);
    dispose = render(() => <HistoryView history={history} runs={4} />, el);
    flush();
    expect(labels(el)).toEqual(['Start', 'Turn', '11 earlier runs', 'call 11', 'call 12', 'call 13', 'call 14']);
    const more = [...el.querySelectorAll<HTMLButtonElement>('button')].find((b) => b.textContent === 'Show all')!;
    expect(more.getAttribute('aria-pressed')).toBe('false');
    expect(more.getAttribute('aria-label')).toBe('Show all runs of Turn');
    more.click();
    flush();
    expect(labels(el)).toHaveLength(2 + 1 + 15);
    expect(labels(el)[2]).toBe('All 15 runs');
    const pressed = [...el.querySelectorAll('button')].find((b) => b.textContent === 'Show all');
    expect(pressed?.getAttribute('aria-pressed')).toBe('true');
  });

  it('speaks the app’s words', async () => {
    const { el } = await mount(2, vi.fn(), {
      labels: { undo: 'Annuler', undoTo: (s) => `Annuler jusqu’à ${s}`, count: (n) => `${n} étapes` },
    });
    expect(button(el, 'Annuler')).toBeDefined();
    expect(el.textContent).toContain('2 étapes');
    expect(rows(el)[1].querySelector('button')?.getAttribute('aria-label')).toBe('Annuler jusqu’à Set a to 1');
  });

  it('says an origin and a time in the app’s words, an icon by kind', async () => {
    const { el } = await mount(1, vi.fn(), {
      labels: {
        origin: (kind) => (kind === 'person' ? 'Vous' : kind),
        description: (origin, at) => `${origin}, ${new Date(at).toISOString().slice(11, 16)}`,
      },
      originIcons: { person: 'bot' },
    });
    expect(rows(el)[1].textContent).toContain('Vous, 10:00');
  });

  it('says what goes here when nothing has happened', async () => {
    const { el } = await mount(0);
    expect(el.textContent).toContain('Nothing has happened yet.');
    expect(button(el, 'Undo').disabled).toBe(true);
    expect(button(el, 'Redo').disabled).toBe(true);
  });

  it('is accessible (axe): steps done and undone', async () => {
    const { el } = await mount(3);
    button(el, 'Undo').click();
    await settle();
    expect(await axeViolations(el)).toBe('');
  });
});
