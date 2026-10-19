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

  it('writes log-only entries once a task, never on show or said, each source within its bound', async () => {
    const n = new NotificationsManager({ historySize: 3 });
    n.post('a message');
    const changes: number[] = [];
    const sub = n.history.subscribe((h) => changes.push(h.length));
    changes.length = 0;
    for (let i = 0; i < 5; i++) n.log(`Created ${i}`, 'info', 'Transforms', { detail: `${i} ms`, keep: 2 });
    expect(changes).toEqual([]);
    await Promise.resolve();
    // One change for the burst; the two newest transforms, and the message is still there.
    expect(changes).toHaveLength(1);
    expect(n.history.value.map((m) => m.text)).toEqual(['Created 4', 'Created 3', 'a message']);
    expect(n.history.value[0]).toMatchObject({ source: 'Transforms', detail: '4 ms' });
    expect(n.current.value?.text).toBe('a message');
    expect(n.spoken.value).toBe('a message');
    // A message posted meanwhile comes after what was logged before it; messages keep their own bound.
    n.log('Created 5', 'info', 'Transforms', { keep: 2 });
    for (const t of ['b', 'c', 'd']) n.post(t);
    expect(n.history.value.map((m) => m.text)).toEqual(['d', 'c', 'b', 'Created 5', 'Created 4']);
    sub.unsubscribe();
  });
});
