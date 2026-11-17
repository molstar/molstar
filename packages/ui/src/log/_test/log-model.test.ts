/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { boundLog, type LogEntry, LogModel, LogAll, LogProblems } from '../log-model.js';

// The log in plain Node (log §1): a burst is one change, one bound with the problems dropped last, the sources as
// written, the filters, the data out.

const texts = (log: LogModel) => log.entries.value.map((e) => e.text);

describe('LogModel', () => {
  it('writes a burst once a task, oldest first, each with its id and time', async () => {
    const log = new LogModel({ now: () => 42 });
    const changes: number[] = [];
    const stop = log.entries.listen!((list) => changes.push(list.length));
    const first = log.write({ text: 'Created 0', source: 'Transforms', detail: '3 ms' });
    for (let i = 1; i < 5; i++) log.write({ text: `Created ${i}`, source: 'Transforms' });
    expect(changes).toEqual([]);
    expect(first).toMatchObject({ id: 1, kind: 'info', time: 42, detail: '3 ms' });
    await Promise.resolve();
    expect(changes).toEqual([5]);
    expect(texts(log)).toEqual(['Created 0', 'Created 1', 'Created 2', 'Created 3', 'Created 4']);
    stop();
  });

  it('writes at once on flush or save, and the scheduled write then finds nothing', async () => {
    const log = new LogModel();
    log.write({ text: 'a' });
    log.flush();
    expect(texts(log)).toEqual(['a']);
    log.write({ text: 'b', kind: 'error' });
    expect(log.save().map((e) => [e.text, e.kind])).toEqual([
      ['a', 'info'],
      ['b', 'error'],
    ]);
    await Promise.resolve();
    expect(texts(log)).toEqual(['a', 'b']);
    // The data is plain: it survives JSON.
    expect(JSON.parse(JSON.stringify(log.save()))).toEqual(log.save());
  });

  it('keeps the warnings and errors within the bound while anything else is left to drop', () => {
    const log = new LogModel({ size: 300 });
    for (let i = 0; i < 5; i++) log.write({ text: `message ${i}` });
    log.write({ text: 'No such chain', kind: 'warning' });
    log.write({ text: 'Created Model', source: 'Transforms' });
    log.flush();
    for (let i = 0; i < 400; i++) log.write({ text: `Updated ${i}`, source: 'Transforms' });
    log.flush();
    expect(log.entries.value).toHaveLength(300);
    expect(log.entries.value[0].text).toBe('No such chain');
    expect(log.entries.value.at(-1)?.text).toBe('Updated 399');
    // Problems alone fill the bound: the newest kept.
    for (let i = 0; i < 300; i++) log.write({ text: `warning ${i}`, kind: 'warning', source: 'Plugin' });
    log.flush();
    const list = log.entries.value;
    expect(list).toHaveLength(300);
    expect(list.every((e) => e.kind === 'warning')).toBe(true);
    expect(list.at(-1)?.text).toBe('warning 299');
    expect(list[0].text).toBe('warning 0');
  });

  it('is one bound for every source', () => {
    const log = new LogModel({ size: 10 });
    log.write({ text: 'Plugin says so', kind: 'warning', source: 'Plugin' });
    log.write({ text: 'Failed to load', kind: 'error' });
    for (let i = 0; i < 20; i++) log.write({ text: `message ${i}` });
    log.flush();
    expect(texts(log)).toEqual([
      'Plugin says so',
      'Failed to load',
      ...Array.from({ length: 8 }, (_, i) => `message ${12 + i}`),
    ]);
  });

  it('lists the sources as first written, told only when one is new', () => {
    const log = new LogModel();
    const told: (readonly string[])[] = [];
    const stop = log.sources.listen!((s) => told.push(s));
    log.write({ text: 'a', source: 'Transforms' });
    log.write({ text: 'b', source: 'Plugin' });
    log.write({ text: 'c' });
    log.flush();
    log.write({ text: 'd', source: 'Transforms' });
    log.flush();
    expect(told).toEqual([['Transforms', 'Plugin']]);
    stop();
  });

  it('filters by kind, source and text', () => {
    const log = new LogModel();
    log.write({ text: 'Colored the polymer', detail: '12 ms' });
    log.write({ text: 'No such chain', kind: 'warning' });
    log.write({ text: 'Created Model', source: 'Transforms', detail: '3 ms' });
    log.write({ text: 'Created Structure', source: 'Transforms' });
    log.write({ text: 'Recomputed', kind: 'debug', source: 'Plugin' });
    log.flush();
    const shown = () => log.filtered.value.map((e) => e.text);
    expect(shown()).toHaveLength(5);
    log.setSource('Transforms');
    expect(shown()).toEqual(['Created Model', 'Created Structure']);
    log.setSource(null);
    expect(shown()).toEqual(['Colored the polymer', 'No such chain']);
    log.setFilter(LogProblems);
    expect(shown()).toEqual(['No such chain']);
    // Kinds toggle: from every kind, one; then added, taken out; the last out is every kind.
    log.setFilter(LogAll);
    log.toggleKind('info');
    expect(shown()).toEqual(['Colored the polymer', 'Created Model', 'Created Structure']);
    log.toggleKind('warning');
    expect(log.filter.value).toEqual({ kinds: ['info', 'warning'] });
    log.setSource('Transforms');
    expect(shown()).toEqual(['Created Model', 'Created Structure']);
    log.setSource(undefined);
    log.toggleKind('info');
    log.toggleKind('warning');
    expect(log.filter.value).toEqual({});
    log.setFilter(LogAll);
    log.setQuery('  MODEL ');
    expect(shown()).toEqual(['Created Model']);
    // The detail and the source match too.
    log.setQuery('12 ms');
    expect(shown()).toEqual(['Colored the polymer']);
    log.setQuery('plugin');
    expect(shown()).toEqual(['Recomputed']);
  });

  it('clears the entries and the sources; a source filter goes back to all', () => {
    const log = new LogModel();
    log.write({ text: 'a', source: 'X' });
    log.flush();
    log.setSource('X');
    log.write({ text: 'pending' });
    log.clear();
    expect(log.entries.value).toEqual([]);
    expect(log.sources.value).toEqual([]);
    expect(log.filter.value).toEqual({});
    log.flush();
    expect(log.entries.value).toEqual([]);
    log.write({ text: 'b', source: 'X' });
    log.flush();
    expect(log.sources.value).toEqual(['X']);
  });

  it('takes an injected scheduler, and writes nothing once disposed', () => {
    const queued: (() => void)[] = [];
    const log = new LogModel({ schedule: (fn) => queued.push(fn) });
    log.write({ text: 'a' });
    log.write({ text: 'b' });
    expect(queued).toHaveLength(1);
    queued[0]();
    expect(texts(log)).toEqual(['a', 'b']);
    log.dispose();
    log.write({ text: 'c' });
    expect(queued).toHaveLength(1);
    expect(texts(log)).toEqual(['a', 'b']);
  });
});

describe('boundLog', () => {
  const e = (id: number, kind: LogEntry['kind']): LogEntry => ({ id, kind, text: String(id), time: 0 });

  it('returns the list as it is within the bound', () => {
    const list = [e(1, 'info'), e(2, 'error')];
    expect(boundLog(list, 2)).toBe(list);
  });

  it('drops the oldest non-problems first, then the oldest problems', () => {
    const list = [e(1, 'error'), e(2, 'info'), e(3, 'warning'), e(4, 'debug'), e(5, 'success')];
    expect(boundLog(list, 3).map((x) => x.id)).toEqual([1, 3, 5]);
    expect(boundLog(list, 1).map((x) => x.id)).toEqual([3]);
  });
});
