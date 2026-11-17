/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush, render } from '../../dom/index.js';
import { shownInOrder } from '../../testing/virtual-list.js';
import { axeViolations } from '../../testing/axe.js';
import { HistoryModel, type HistoryPort } from '../history-model.js';
import { historyRows, HistoryView } from '../history-view.js';

// History's view (history §4): the steps oldest first, the position, Undo and Redo, to a step, unsafe, blocked.

afterEach(() => {
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

async function mount(n: number, onResult = vi.fn()) {
  const s = store();
  const history = new HistoryModel(s.port, { now: () => Date.UTC(2026, 0, 1, 10, 0, 0) });
  for (let i = 1; i <= n; i++) await history.run(`Set a to ${i}`, person, async () => void s.set(i));
  const el = document.createElement('div');
  document.body.append(el);
  render(() => <HistoryView history={history} onResult={onResult} />, el);
  flush();
  return { s, history, el, onResult };
}

const rows = (el: HTMLElement) => shownInOrder(el, '[data-state]');
const labels = (el: HTMLElement) => rows(el).map((r) => r.querySelector('[data-tip-part]')?.textContent);
const states = (el: HTMLElement) => rows(el).map((r) => r.getAttribute('data-state'));
const button = (el: HTMLElement, name: string) =>
  [...el.querySelectorAll<HTMLButtonElement>('button')].find(
    (b) => b.getAttribute('aria-label') === name || b.getAttribute('data-tip') === name,
  )!;

describe('historyRows', () => {
  it('lists the start, the steps oldest first with a group’s runs, then the undone ones', () => {
    const step = (id: number, runs: string[]) => ({ id, label: `s${id}`, origin: person, time: 0, runs });
    const r = historyRows([step(2, ['a', 'b']), step(1, ['s1'])], [step(3, ['s3'])]);
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
    // "Undo to here" on the first step: it stays, the three after it are undone.
    (rows(el)[1].querySelector('button[data-tip="Undo to here"]') as HTMLButtonElement).click();
    await settle();
    expect(s.get()).toBe(1);
    expect(states(el)).toEqual(['done', 'current', 'undone', 'undone', 'undone']);
    (rows(el)[3].querySelector('button[data-tip="Redo to here"]') as HTMLButtonElement).click();
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
    (rows(el)[0].querySelector('button') as HTMLButtonElement).click();
    await settle();
    expect(onResult).toHaveBeenCalledWith({ ok: false, reason: 'Changed outside the history' });
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
