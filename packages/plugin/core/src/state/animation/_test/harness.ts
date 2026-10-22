/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Playback without WebGL: a plugin whose canvas is a fake that records what reaches it, and a clock that ticks the
// animation manager synchronously (`tick(t, true)` awaits each frame). What a tick did is read from the counters.

import { vi } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import type { PluginSpec } from '@molstar/plugin/spec';
import type { PluginState } from '@molstar/plugin/state';
import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { RawData } from '@molstar/plugin/state/transforms/data/fetch';

export const TickMs = 1000 / 60;

export async function createPlaybackPlugin(spec: PluginSpec = { behaviors: [] }) {
  const plugin = new PluginContext(spec);
  await plugin.init();
  const canvas3d = {
    props: {},
    setProps: vi.fn(),
    requestCameraReset: vi.fn(),
    requestDraw: vi.fn(),
    add: vi.fn(),
    remove: vi.fn(),
    update: vi.fn(),
    mark: vi.fn(),
    commit: vi.fn(),
    syncVisibility: vi.fn(),
  };
  (plugin as any).canvas3d = canvas3d;
  const counts = { setSnapshot: 0, changed: 0, updated: 0, recreated: 0, created: 0, removed: 0 };
  const setSnapshot = plugin.state.data.setSnapshot.bind(plugin.state.data);
  plugin.state.data.setSnapshot = ((...a: Parameters<typeof setSnapshot>) => {
    counts.setSnapshot++;
    return setSnapshot(...a);
  }) as typeof setSnapshot;
  plugin.managers.snapshot.events.changed.subscribe(() => counts.changed++);
  plugin.state.data.events.object.created.subscribe(() => counts.created++);
  plugin.state.data.events.object.removed.subscribe(() => counts.removed++);
  plugin.state.data.events.object.updated.subscribe((e) =>
    e.action === 'recreate' ? counts.recreated++ : counts.updated++,
  );
  return { plugin, canvas3d, counts };
}

/**
 * The data tree with `n` raw-data nodes labelled `label0`, `label1`, ...: made on the first call, updated on the next
 * ones, so trees with other labels differ in params only.
 */
export async function rawDataTree(plugin: PluginContext, n: number, label: string) {
  const b = plugin.build();
  for (let i = 0; i < n; i++) {
    const params = { data: `${i}`, label: `${label}${i}` };
    if (plugin.state.data.cells.has(`node${i}`)) b.to(`node${i}`).update(params);
    else b.toRoot().apply(RawData, params, { ref: `node${i}` });
  }
  await b.commit();
  return plugin.state.data.getSnapshot();
}

/** Adds a snapshot with these transition frames and makes it current. */
export function addTransition(
  plugin: PluginContext,
  frames: PluginState.StateTransition['frames'],
  extra: Partial<PluginState.Snapshot> = {},
) {
  const snapshot: PluginState.Snapshot = {
    id: 'transition' as any,
    data: frames[0].data,
    transition: { frames },
    ...extra,
  };
  plugin.managers.snapshot.add(PluginStateSnapshotManager.Entry(snapshot, { name: 'transition' }));
  plugin.managers.snapshot.setCurrent('transition');
  return snapshot;
}

/** Ticks the animation manager at 60 Hz from 0 until `untilMs` or until the animation stops; returns the ticks run. */
export async function play(plugin: PluginContext, untilMs: number) {
  let ticks = 0;
  for (let t = 0; t <= untilMs; t += TickMs) {
    await plugin.managers.animation.tick(t, true);
    ticks++;
    if (!plugin.managers.animation.isAnimating) break;
  }
  return ticks;
}
