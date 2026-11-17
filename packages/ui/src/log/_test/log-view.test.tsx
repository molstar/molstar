/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { DomHost, flush, render } from '../../dom/index.js';
import { shownInOrder } from '../../testing/virtual-list.js';
import { axeViolations } from '../../testing/axe.js';
import { LogModel } from '../log-model.js';
import { LogView } from '../log-view.js';

// The log's view (log §3): recycled one-line rows, newest last and kept at the bottom, the filters and the text filter.

const host = DomHost.requestAnimationFrame;
let queue: ((t: number) => void)[] = [];
/** Runs the frames asked for. */
const tick = () => {
  for (let i = 0; i < 5 && queue.length; i++) {
    const q = queue;
    queue = [];
    for (const f of q) f(DomHost.now());
  }
  flush();
};

beforeEach(() => {
  DomHost.requestAnimationFrame = (f) => void queue.push(f);
});
afterEach(() => {
  document.body.innerHTML = '';
  tick();
  DomHost.requestAnimationFrame = host;
  queue = [];
});

function mount(log: LogModel) {
  const el = document.createElement('div');
  document.body.append(el);
  const dispose = render(() => <LogView log={log} />, el);
  return { el, dispose };
}

const texts = (el: HTMLElement) => {
  flush();
  return shownInOrder(el, '[data-kind]').map((e) => e.querySelector('[data-tip-part]')?.textContent ?? '');
};

const pick = (el: HTMLElement, label: string) =>
  ([...el.querySelectorAll('button')].find((b) => b.getAttribute('aria-label') === label) as HTMLButtonElement).click();

describe('LogView', () => {
  it('lists the entries oldest first, newest last, each a one-line row with the side tooltip', () => {
    const log = new LogModel();
    const { el } = mount(log);
    log.write({ text: 'first' });
    log.write({ text: 'second', kind: 'error', source: 'Plugin', detail: '3 ms' });
    log.flush();
    expect(texts(el)).toEqual(['first', 'second']);
    const row = shownInOrder(el, '[data-kind]')[1];
    expect(row.getAttribute('data-kind')).toBe('error');
    expect(row.hasAttribute('data-tip-row')).toBe(true);
    expect(row.textContent).toContain('3 ms · Plugin');
    const list = el.querySelector('[role=log]')!;
    expect(list.getAttribute('aria-live')).toBe('off');
  });

  it('renders a window of many entries in recycled rows, kept at the bottom as entries come', async () => {
    const log = new LogModel();
    const { el } = mount(log);
    for (let i = 0; i < 300; i++) log.write({ text: `Created ${i}`, source: 'Transforms' });
    await Promise.resolve();
    tick();
    const rows = [...el.querySelectorAll('[data-kind]')];
    expect(rows.length).toBeGreaterThan(20);
    expect(rows.length).toBeLessThan(120);
    expect(texts(el).at(-1)).toBe('Created 299');
    expect(texts(el)).not.toContain('Created 0');
    for (let i = 0; i < 10; i++) log.write({ text: `Updated ${i}`, source: 'Transforms' });
    await Promise.resolve();
    tick();
    expect(texts(el).at(-1)).toBe('Updated 9');
    // The same rows, re-pointed.
    expect([...el.querySelectorAll('[data-kind]')].every((r) => rows.includes(r))).toBe(true);
  });

  it('stays where it is while scrolled up, and follows again back at the end', async () => {
    const log = new LogModel();
    const { el } = mount(log);
    for (let i = 0; i < 300; i++) log.write({ text: `Created ${i}` });
    await Promise.resolve();
    tick();
    const viewport = el.querySelector('[role=log]')!.parentElement as HTMLElement;
    const at = (top: number) => {
      Object.defineProperty(viewport, 'scrollTop', { value: top, configurable: true, writable: true });
      Object.defineProperty(viewport, 'clientHeight', { value: 240, configurable: true });
      Object.defineProperty(viewport, 'scrollHeight', { value: 300 * 24, configurable: true });
      viewport.dispatchEvent(new Event('scroll'));
      tick();
    };
    at(0);
    log.write({ text: 'Newest' });
    await Promise.resolve();
    tick();
    // The bound (300) dropped the oldest; the window stayed at the top.
    expect(texts(el)[0]).toBe('Created 1');
    expect(texts(el)).not.toContain('Newest');
    at(300 * 24 - 240);
    log.write({ text: 'Newer still' });
    await Promise.resolve();
    tick();
    expect(texts(el).at(-1)).toBe('Newer still');
  });

  it('filters by kind (toggles), by source (the messages, a source) and by text; Escape empties the text', () => {
    const log = new LogModel();
    const { el } = mount(log);
    log.write({ text: 'Colored the polymer', detail: '12 ms' });
    log.write({ text: 'No such chain', kind: 'warning' });
    log.write({ text: 'Created Model', source: 'Transforms', detail: '3 ms' });
    log.write({ text: 'Created Structure', source: 'Transforms' });
    log.flush();
    flush();
    expect(texts(el)).toHaveLength(4);
    const pressed = () => [...el.querySelectorAll('[aria-pressed=true]')].map((b) => b.getAttribute('aria-label'));
    expect(pressed()).toEqual(['All']);
    pick(el, 'Transforms');
    expect(log.filter.value).toEqual({ source: 'Transforms' });
    expect(texts(el)).toEqual(['Created Model', 'Created Structure']);
    pick(el, 'Messages');
    expect(texts(el)).toEqual(['Colored the polymer', 'No such chain']);
    pick(el, 'Warnings');
    expect(texts(el)).toEqual(['No such chain']);
    expect(pressed()).toEqual(['Warnings', 'Messages']);
    pick(el, 'Info');
    expect(texts(el)).toEqual(['Colored the polymer', 'No such chain']);
    // A pressed source again: every source.
    pick(el, 'Messages');
    expect(pressed()).toEqual(['Info', 'Warnings']);
    expect(texts(el)).toHaveLength(4);

    pick(el, 'All');
    const input = el.querySelector('input') as HTMLInputElement;
    input.value = 'model';
    input.dispatchEvent(new Event('input', { bubbles: true }));
    expect(texts(el)).toEqual(['Created Model']);
    input.value = 'nothing like it';
    input.dispatchEvent(new Event('input', { bubbles: true }));
    flush();
    expect(el.textContent).toContain('Nothing matches the filter.');
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(log.query.value).toBe('');
    expect(texts(el)).toHaveLength(4);
  });

  it('says when nothing is logged yet', () => {
    const { el } = mount(new LogModel());
    flush();
    expect(el.textContent).toContain('Nothing is logged yet.');
  });

  it('is accessible (axe): empty, and with entries filtered', async () => {
    const log = new LogModel();
    const { el } = mount(log);
    flush();
    expect(await axeViolations(el)).toBe('');
    log.write({ text: 'loaded', source: 'Plugin' });
    log.write({ text: 'failed', kind: 'error', source: 'Plugin', detail: 'no such file' });
    tick();
    pick(el, 'Errors');
    tick();
    expect(await axeViolations(el)).toBe('');
  });
});
