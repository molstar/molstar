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
