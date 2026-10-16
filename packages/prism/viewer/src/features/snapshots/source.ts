/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { PluginState } from '@molstar/plugin/state';
import type { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { merge } from 'rxjs';
import { fromBehavior, fromChanged, SourceScope } from '../../sources/state.js';

// The snapshots as the views read them (state.md §4, row 2: a mutable field with a void `changed`). Sources by role
// (README, Folders): this one is allowed to read the plugin, as `sources/` is, so the feature's views never do. The
// snapshot area owns it and hands it to the views (`api.snapshot.state`).

/** A saved snapshot as listed: what the list and the playback card show. */
export interface SnapshotInfo {
  readonly id: string;
  /** The entry's name; empty when it has none. */
  readonly name: string;
  readonly key?: string;
  readonly description?: string;
  /** The description is plain text, not markdown (the entry's `descriptionFormat`). */
  readonly plainText?: true;
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

/**
 * The current snapshot's transition (Mol*'s `snapshot.transition`: frames played after it is applied), as the canvas
 * scrubber shows it. Apart from the list's state: it changes on every frame while it plays.
 */
export interface SnapshotTransitionState {
  /** Its length (ms); 0 when the current snapshot has none. */
  readonly duration: number;
  /** Its shortest frame (ms): the scrubber's step. */
  readonly step: number;
  readonly loop: boolean;
  /** Where it is (ms). */
  readonly time: number;
  /** It plays (the plugin's animation manager runs it). */
  readonly playing: boolean;
}

export interface SnapshotsSource {
  readonly state: ReadonlyState<SnapshotsState>;
  readonly transition: ReadonlyState<SnapshotTransitionState>;
  dispose(): void;
}

const sameTransition = (a: SnapshotTransitionState, b: SnapshotTransitionState) =>
  a.duration === b.duration && a.step === b.step && a.loop === b.loop && a.time === b.time && a.playing === b.playing;

/** The transition of the current snapshot, where it is and whether it plays. */
export function readTransition(plugin: PluginContext): SnapshotTransitionState {
  const manager = plugin.managers.snapshot;
  const animation = plugin.managers.animation;
  const snapshot = manager.current?.snapshot;
  const duration = snapshot?.transition ? (PluginState.getStateTransitionDuration(snapshot) ?? 0) : 0;
  return {
    duration,
    step: PluginState.getMinFrameDuration(snapshot),
    loop: !!snapshot?.transition?.loop,
    time: Math.min(manager.state.currentAnimationTimeMs ?? 0, duration),
    playing: animation.isAnimatingStateTransition && animation.isAnimating,
  };
}

const sameInfo = (a: SnapshotInfo, b: SnapshotInfo) =>
  a.id === b.id &&
  a.name === b.name &&
  a.key === b.key &&
  a.description === b.description &&
  a.plainText === b.plainText &&
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
  // The last read and what it read from: the manager emits `changed` on every frame of a transition and step of the
  // scrubber, with the entries (an immutable list), the current one and `isPlaying` as they were.
  let last: { entries: unknown; current: string | undefined; playing: boolean; state: SnapshotsState } | undefined;
  const read = (): SnapshotsState => {
    const { entries, current, isPlaying } = manager.state;
    if (last && last.entries === entries && last.current === current && last.playing === isPlaying) return last.state;
    const infos = entries.toArray().map((e): SnapshotInfo => {
      const image = imageOf(e);
      return {
        id: e.snapshot.id,
        name: e.name ?? '',
        ...(e.key ? { key: e.key } : {}),
        ...(e.description ? { description: e.description } : {}),
        ...(e.descriptionFormat === 'plaintext' ? { plainText: true as const } : {}),
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
    const state: SnapshotsState = { entries: infos, ...(current ? { current } : {}), playing: isPlaying };
    last = { entries, current, playing: isPlaying, state };
    return state;
  };
  const state = fromChanged(scope, read, manager.events.changed, sameState);
  const transition = fromChanged(
    scope,
    () => readTransition(plugin),
    merge(manager.events.changed, plugin.managers.animation.events.updated),
    sameTransition,
  );
  scope.add(() => {
    for (const url of urls.values()) URL.revokeObjectURL(url);
    urls.clear();
  });
  return { state, transition, dispose: () => scope.dispose() };
}

/**
 * The plugin's markdown extensions as the descriptions read them (parity.md §4; v5's `controls/markdown.tsx`): which links and
 * image sources are commands (`!name=value&…`, by the plugin's parsers), images resolved by the plugin's resolvers (an
 * MVSX archive's files), and the audio player a command made. Running a command is an operation (`snapshot.link`).
 */
export interface MarkdownSource {
  /** A link's or an image's commands, or undefined when it is not a command. */
  command(href: string): Readonly<Record<string, string>> | undefined;
  /** An image's URL from the resolvers (a data URI of an archive's file), or undefined to load it as it is. */
  resolve(src: string): string | Promise<string> | undefined;
  /** The plugin's audio player, once a command made one. */
  readonly audio: ReadonlyState<HTMLAudioElement | null>;
  /** Stops it at the start (Mol*'s `audio.stop`). */
  stopAudio(): void;
  /** Stops it and lets it go (Mol*'s `audio.dispose`): the player is gone until a command makes another. */
  closeAudio(): void;
  dispose(): void;
}

export function createMarkdownSource(plugin: PluginContext): MarkdownSource {
  const manager = plugin.managers.markdownExtensions;
  const scope = new SourceScope();
  return {
    command: (href) => manager.parseArgs(href),
    resolve: (src) => manager.tryResolveUri(src),
    audio: fromBehavior(scope, manager.state.audioPlayer),
    stopAudio: () => manager.audio.stop(),
    closeAudio: () => manager.audio.dispose(),
    dispose: () => scope.dispose(),
  };
}
