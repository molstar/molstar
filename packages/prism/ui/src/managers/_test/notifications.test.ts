/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { NotificationsManager } from '../notifications.js';

// Tasks from two writers: one source mirrored whole (`setTasks`: the plugin's tasks) and tasks put one by one
// (`putTask`: an extension loading); the mirror does not wipe the others.

describe('the notifications tasks', () => {
  it('keeps a task put one by one when another source sets its own, until it is removed', () => {
    const n = new NotificationsManager();
    const ids = () => n.tasks.value.map((t) => t.id);
    n.putTask({ id: 'extension:video', label: 'Loading video export' });
    n.setTasks([{ id: '1', label: 'Downloading' }]);
    expect(ids()).toEqual(['1', 'extension:video']);
    n.setTasks([]);
    expect(ids()).toEqual(['extension:video']);
    n.removeTask('extension:video');
    n.setTasks([{ id: '2', label: 'Parsing' }]);
    expect(ids()).toEqual(['2']);
  });
});

describe('the log', () => {
  it("keeps a message's detail for the log, not the status line's text", () => {
    const n = new NotificationsManager();
    n.post('Colored the polymer', 'info', { detail: '12 ms' });
    expect(n.current.value).toMatchObject({ text: 'Colored the polymer', detail: '12 ms' });
    expect(n.spoken.value).toBe('Colored the polymer');
  });

  it('writes log-only entries once a task, never on show or said, all within one bound', async () => {
    const n = new NotificationsManager({ historySize: 3, logSize: 4 });
    n.post('a message');
    const changes: number[] = [];
    const sub = n.logged.subscribe((h) => changes.push(h.length));
    changes.length = 0;
    for (let i = 0; i < 5; i++) n.log(`Created ${i}`, 'info', 'Transforms', { detail: `${i} ms` });
    expect(changes).toEqual([]);
    await Promise.resolve();
    // One change for the burst; the newest four entries, whatever wrote them.
    expect(changes).toHaveLength(1);
    expect(n.logged.value.map((m) => m.text)).toEqual(['Created 4', 'Created 3', 'Created 2', 'Created 1']);
    expect(n.history.value.map((m) => m.text)).toEqual(['a message']);
    expect(n.logged.value[0]).toMatchObject({ source: 'Transforms', detail: '4 ms' });
    expect(n.current.value?.text).toBe('a message');
    expect(n.spoken.value).toBe('a message');
    // A message posted meanwhile comes after what was logged before it.
    n.log('Created 5', 'info', 'Transforms');
    for (const t of ['b', 'c']) n.post(t);
    expect(n.logged.value.map((m) => m.text)).toEqual(['c', 'b', 'Created 5', 'Created 4']);
    expect(n.history.value.map((m) => m.text)).toEqual(['c', 'b', 'a message']);
    sub.unsubscribe();
  });

  it('keeps the warnings and errors within the bound while anything else is left to drop (perf audit 11/B5)', async () => {
    const n = new NotificationsManager({ logSize: 300 });
    for (let i = 0; i < 5; i++) n.post(`message ${i}`);
    n.post('No such chain', 'warning');
    n.log('Created Model', 'info', 'Transforms');
    await Promise.resolve();
    // A load's burst of transforms: the problem stays, the oldest of the rest go.
    for (let i = 0; i < 400; i++) n.log(`Updated ${i}`, 'info', 'Transforms');
    await Promise.resolve();
    expect(n.logged.value).toHaveLength(300);
    expect(n.logged.value.map((m) => m.text)).toContain('No such chain');
    expect(n.logged.value.at(-1)?.text).toBe('No such chain');
    // 300 plugin warnings (the audit's case): the bound holds, the newest problems kept, the rest gone first.
    for (let i = 0; i < 300; i++) n.log(`plugin warning ${i}`, 'warning', 'Plugin');
    await Promise.resolve();
    const logged = n.logged.value;
    expect(logged).toHaveLength(300);
    expect(logged.every((m) => m.kind === 'warning')).toBe(true);
    expect(logged[0].text).toBe('plugin warning 299');
    // The sources written, a filter each: told once a source is new.
    expect(n.logSources.value).toEqual(['Plugin', 'Transforms']);
  });

  it('is one bound for every source: a burst of one fills it with the problems of all kept', async () => {
    const n = new NotificationsManager({ logSize: 10 });
    n.log('Plugin says so', 'warning', 'Plugin');
    n.post('Failed to load', 'error');
    for (let i = 0; i < 20; i++) n.post(`message ${i}`);
    await Promise.resolve();
    const texts = n.logged.value.map((m) => m.text);
    expect(texts).toHaveLength(10);
    expect(texts.slice(0, 8)).toEqual(Array.from({ length: 8 }, (_, i) => `message ${19 - i}`));
    expect(texts.slice(8)).toEqual(['Failed to load', 'Plugin says so']);
  });
});
