/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { createSettingsSources } from '../settings.js';

// Sources that follow the plugin only while watched (prism AGENTS.md: hidden UI holds no subscriptions), and only on
// the events that change them.

const frame = () => new Promise<void>((r) => setTimeout(r, 20));

describe('followed sources', () => {
  it("reads the settings' groups only while they are watched", async () => {
    const v = await headlessViewer();
    const behaviors = v.plugin.state.behaviors.events.changed;
    const listeners = () => behaviors.observers.length;
    const before = listeners();
    const { groups } = createSettingsSources(v.plugin);
    expect(listeners()).toBe(before);
    expect(groups.value.length).toBeGreaterThan(0);
    const seen: unknown[] = [];
    const sub = groups.subscribe((g) => seen.push(g));
    expect(listeners()).toBe(before + 1);
    expect(seen).toHaveLength(1);
    sub.unsubscribe();
    expect(listeners()).toBe(before);
  });

  it('does not count the selection again on a state change that leaves it as it is', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.select.set({ target: 'A:1-10' }));
    const { selection } = v.sources.interaction;
    const seen: (string | undefined)[] = [];
    const sub = selection.subscribe((s) => seen.push(s));
    expect(seen.at(-1)).toMatch(/10 residues/);
    const reads = vi.spyOn(v.plugin.managers.structure.selection.entries, 'values');
    v.plugin.state.data.events.changed.next({ state: v.plugin.state.data, inTransaction: false });
    await frame();
    expect(reads).not.toHaveBeenCalled();

    expectOk(await v.api.select.set({ target: 'A:1-20' }));
    await frame();
    expect(seen.at(-1)).toMatch(/20 residues/);
    sub.unsubscribe();
  });
});
