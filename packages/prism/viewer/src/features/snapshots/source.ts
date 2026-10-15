/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { fromChanged, SourceScope } from '../../sources/state.js';

// The snapshots as the views read them (state.md §4, row 2: a mutable field with a void `changed`). Sources by role
// (contexts.md §7): this one is allowed to read the plugin, as `sources/` is, so the feature's views never do. The
// snapshot area owns it and hands it to the views (`api.snapshot.state`).

/** A saved snapshot as listed: what the list and the playback card show. */
export interface SnapshotInfo {
  readonly id: string;
  /** The entry's name; empty when it has none. */
  readonly name: string;
  readonly key?: string;
  readonly description?: string;
  /** When it was taken (ms since the epoch). */
  readonly timestamp: number;
  /** How long playback shows it (ms), when it says so. */
  readonly showFor?: number;
  /** A CSS `url(…)` of its thumbnail, when it was saved with an image. */
  readonly image?: string;
}

export interface SnapshotsState {
  readonly entries: readonly SnapshotInfo[];
  /** The one applied last. */
  readonly current?: string;
  readonly playing: boolean;
}

export interface SnapshotsSource {
  readonly state: ReadonlyState<SnapshotsState>;
  dispose(): void;
}

const sameInfo = (a: SnapshotInfo, b: SnapshotInfo) =>
  a.id === b.id &&
  a.name === b.name &&
  a.key === b.key &&
  a.description === b.description &&
  a.timestamp === b.timestamp &&
  a.showFor === b.showFor &&
  a.image === b.image;

function sameState(a: SnapshotsState, b: SnapshotsState): boolean {
  return (
    a.current === b.current &&
    a.playing === b.playing &&
    a.entries.length === b.entries.length &&
    a.entries.every((e, i) => sameInfo(e, b.entries[i]))
  );
}

export function createSnapshotsSource(plugin: PluginContext): SnapshotsSource {
  const manager = plugin.managers.snapshot;
  const scope = new SourceScope();
  // Thumbnails are blobs of the asset manager: one object URL per asset, released when its entry goes.
  const urls = new Map<string, string>();
  const imageOf = (e: PluginStateSnapshotManager.Entry): string | undefined => {
    if (!e.image || typeof URL === 'undefined' || typeof URL.createObjectURL !== 'function') return undefined;
    let url = urls.get(String(e.image.id));
    if (!url) {
      const file = plugin.managers.asset.get(e.image)?.file;
      if (!file) return undefined;
      urls.set(String(e.image.id), (url = URL.createObjectURL(file)));
    }
    return `url("${url}")`;
  };
  const read = (): SnapshotsState => {
    const { entries, current, isPlaying } = manager.state;
    const infos = entries.toArray().map((e): SnapshotInfo => {
      const image = imageOf(e);
      return {
        id: e.snapshot.id,
        name: e.name ?? '',
        ...(e.key ? { key: e.key } : {}),
        ...(e.description ? { description: e.description } : {}),
        timestamp: e.timestamp,
        ...(e.snapshot.durationInMs !== undefined ? { showFor: e.snapshot.durationInMs } : {}),
        ...(image ? { image } : {}),
      };
    });
    if (urls.size) {
      const alive = new Set(entries.toArray().flatMap((e) => (e.image ? [String(e.image.id)] : [])));
      for (const [id, url] of urls) {
        if (alive.has(id)) continue;
        URL.revokeObjectURL(url);
        urls.delete(id);
      }
    }
    return { entries: infos, ...(current ? { current } : {}), playing: isPlaying };
  };
  const state = fromChanged(scope, read, manager.events.changed, sameState);
  scope.add(() => {
    for (const url of urls.values()) URL.revokeObjectURL(url);
    urls.clear();
  });
  return { state, dispose: () => scope.dispose() };
}
