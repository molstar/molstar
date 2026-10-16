/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '@molstar/plugin/commands';
import { describe, expect, it } from 'vitest';
import { headlessViewer } from '../../_test/harness.js';
import { toastText, watchToasts } from '../toasts.js';

// The plugin's toasts as status messages (parity.md's "toasts" row).

describe('toasts', () => {
  it('posts each toast once, again when shown again under its key, and stops when released', async () => {
    const v = await headlessViewer();
    const posted: string[] = [];
    const stop = watchToasts(v.plugin, { post: (text) => posted.push(text) });
    await PluginCommands.Toast.Show(v.plugin, { title: 'Assembly', message: 'Showing <b>1</b> of 2', key: 'a' });
    await PluginCommands.Toast.Show(v.plugin, { title: 'Other', message: 'Hello' });
    expect(posted).toEqual(['Assembly: Showing 1 of 2', 'Other: Hello']);
    await PluginCommands.Toast.Show(v.plugin, { title: 'Assembly', message: 'Showing 2 of 2', key: 'a' });
    expect(posted).toHaveLength(3);
    await PluginCommands.Toast.Hide(v.plugin, { key: 'a' });
    expect(posted).toHaveLength(3);
    stop();
    await PluginCommands.Toast.Show(v.plugin, { title: 'Late', message: 'x' });
    expect(posted).toHaveLength(3);
  });

  it('says a toast in words', () => {
    expect(toastText('Title', 'a&amp;b<br/>c')).toBe('Title: a&b c');
    expect(toastText('Title', () => null)).toBe('Title');
    expect(toastText('', 'only')).toBe('only');
  });
});
