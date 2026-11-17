/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { LogModel } from '../../log/log-model.js';
import { StatusModel, type StatusTimers } from '../status-model.js';

// The status in plain Node (status §1): its fade on an injected timer, the live region's text, tasks from two writers
// kept by id, its messages in the log.

/** Timers advanced by hand. */
function manualTimers() {
  let time = 0;
  let serial = 0;
  const waiting = new Map<number, { at: number; fn: () => void }>();
  const timers: StatusTimers = {
    setTimeout: (fn, ms) => {
      waiting.set(++serial, { at: time + ms, fn });
      return serial;
    },
    clearTimeout: (id) => void waiting.delete(id as number),
  };
  const advance = (ms: number) => {
    time += ms;
    for (const [id, t] of [...waiting]) {
      if (t.at > time) continue;
      waiting.delete(id);
      t.fn();
    }
  };
  return { timers, advance, pending: () => waiting.size };
}

describe('StatusModel messages', () => {
  it('shows a message until it fades; an error stays until dismissed', () => {
    const t = manualTimers();
    const status = new StatusModel({ timers: t.timers, fadeMs: 1000 });
    status.post('Loaded 1tqn', 'success');
    expect(status.current.value).toMatchObject({ text: 'Loaded 1tqn', kind: 'success' });
    t.advance(999);
    expect(status.current.value).toBeDefined();
    t.advance(1);
    expect(status.current.value).toBeUndefined();
    const error = status.post('Could not load', 'error');
    t.advance(60_000);
    expect(status.current.value?.id).toBe(error.id);
    status.dismiss(error.id + 1);
    expect(status.current.value?.id).toBe(error.id);
    status.dismiss(error.id);
    expect(status.current.value).toBeUndefined();
  });

  it('restarts the fade for a newer message; the older one fading does not hide it', () => {
    const t = manualTimers();
    const status = new StatusModel({ timers: t.timers, fadeMs: 1000 });
    status.post('a');
    t.advance(800);
    status.post('b');
    t.advance(800);
    expect(status.current.value?.text).toBe('b');
    t.advance(200);
    expect(status.current.value).toBeUndefined();
    expect(t.pending()).toBe(0);
  });

  it('says each message, a repeat again (a no-break space), and announces without showing or logging', () => {
    const status = new StatusModel({ timers: manualTimers().timers });
    status.post('Colored the polymer', 'info', { detail: '12 ms' });
    expect(status.spoken.value).toBe('Colored the polymer');
    status.post('Colored the polymer');
    expect(status.spoken.value).toBe('Colored the polymer ');
    status.announce('3 residues selected');
    expect(status.spoken.value).toBe('3 residues selected');
    expect(status.current.value?.text).toBe('Colored the polymer');
    status.log.flush();
    expect(status.log.entries.value.map((e) => e.text)).toEqual(['Colored the polymer', 'Colored the polymer']);
    expect(status.log.entries.value[0]).toMatchObject({ detail: '12 ms' });
    expect(status.log.entries.value[0].source).toBeUndefined();
  });

  it('writes into a log it is given, after what was written before', () => {
    const log = new LogModel();
    const status = new StatusModel({ log, timers: manualTimers().timers });
    log.write({ text: 'Created 5', source: 'Transforms' });
    status.post('b');
    status.post('c', 'warning');
    log.flush();
    expect(log.entries.value.map((e) => e.text)).toEqual(['Created 5', 'b', 'c']);
    log.setSource(null);
    expect(log.filtered.value.map((e) => e.text)).toEqual(['b', 'c']);
  });

  it('clears its fade when retired', () => {
    const t = manualTimers();
    const status = new StatusModel({ timers: t.timers });
    status.post('a');
    expect(t.pending()).toBe(1);
    status.retire();
    expect(t.pending()).toBe(0);
  });
});

describe('StatusModel tasks', () => {
  it('keeps a task put one by one when another source sets its own, until it is removed', () => {
    const status = new StatusModel();
    const ids = () => status.tasks.value.map((t) => t.id);
    status.putTask({ id: 'extension:video', label: 'Loading video export' });
    status.setTasks([{ id: '1', label: 'Downloading' }]);
    expect(ids()).toEqual(['1', 'extension:video']);
    status.setTasks([]);
    expect(ids()).toEqual(['extension:video']);
    status.removeTask('extension:video');
    status.setTasks([{ id: '2', label: 'Parsing' }]);
    expect(ids()).toEqual(['2']);
  });

  it('tells a task by id on its progress, and the ids only when they change; overlay tasks are not the line’s', () => {
    const status = new StatusModel();
    const one = status.task('1');
    const seen: (number | undefined)[] = [];
    const stopOne = one.listen!((t) => seen.push(t?.progress));
    let idsTold = 0;
    const stopIds = status.taskIds.listen!(() => idsTold++);
    status.setTasks([
      { id: '1', label: 'Loading', progress: 0 },
      { id: 'x', label: 'Export', overlay: true },
    ]);
    expect(status.taskIds.value).toEqual(['1']);
    expect(idsTold).toBe(1);
    for (let i = 1; i <= 3; i++) status.setTasks([{ id: '1', label: 'Loading', progress: i / 10 }]);
    expect(idsTold).toBe(1);
    expect(seen).toEqual([0, 0.1, 0.2, 0.3]);
    status.setTasks([]);
    expect(one.value).toBeUndefined();
    expect(idsTold).toBe(2);
    stopOne();
    stopIds();
  });

  it('replaces a put task in place', () => {
    const status = new StatusModel();
    status.putTask({ id: 'a', label: 'A' });
    status.putTask({ id: 'b', label: 'B' });
    status.putTask({ id: 'a', label: 'A', progress: 0.5 });
    expect(status.tasks.value.map((t) => [t.id, t.progress])).toEqual([
      ['a', 0.5],
      ['b', undefined],
    ]);
  });

  it('sets the chips', () => {
    const status = new StatusModel();
    status.setChips([{ id: 'sel', label: '3 residues' }]);
    expect(status.chips.value.map((c) => c.label)).toEqual(['3 residues']);
  });
});
