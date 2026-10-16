/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { mapState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { MarkdownHost } from '@molstar/prism-ui';
import type { SnapshotStepperModel, SnapshotTransitionModel } from '@molstar/prism-ui/model/playback';
import type { Subscription } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { SnapshotEntry, SnapshotsModel } from '../../model/snapshots.js';
import type { SnapshotApi } from './operations.js';
import type { SnapshotInfo, SnapshotsState } from './source.js';

// The view models of the Snapshots tab and the playback card, over the snapshot area: they read its state (the source
// in this folder) and act through the person's API, never the plugin. DOM-free.

const Transition = PD.MappedStatic(
  'animate',
  {
    animate: PD.Group(
      {
        durationInMs: PD.Numeric(250, { min: 100, max: 5000, step: 50 }, { label: 'Duration (ms)' }),
        easing: PD.Select('linear', [
          ['linear', 'Linear'],
          ['ease-in-out', 'Ease in-out'],
          ['ease-out', 'Ease out'],
        ]),
      },
      { isFlat: true },
    ),
    instant: PD.Group({}, { isFlat: true }),
  },
  {
    label: 'Camera',
    options: [
      ['animate', 'Animate'],
      ['instant', 'Instant'],
    ],
  },
);

/** What a new snapshot includes: Mol*'s snapshot options (`PluginState.SnapshotParams`), as the tab shows them. */
export const SaveOptions = {
  durationInMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Show for (ms)' }),
  data: PD.Boolean(true, { label: 'Data' }),
  camera: PD.Boolean(true, { label: 'Camera' }),
  cameraTransition: Transition,
  canvas3d: PD.Boolean(true, { label: 'Canvas settings' }),
  componentManager: PD.Boolean(true, { label: 'Component options' }),
  structureSelection: PD.Boolean(false, { label: 'Selection' }),
  animation: PD.Boolean(true, { label: 'Animation' }),
  startAnimation: PD.Boolean(false, { label: 'Start animation' }),
  behavior: PD.Boolean(false, { label: 'Behaviors' }),
  interactivity: PD.Boolean(true, { label: 'Interactivity' }),
  image: PD.Boolean(false, { label: 'Thumbnail', description: 'Save a small image with each snapshot.' }),
};

/** A snapshot's own playback settings. */
export const EntryOptions = {
  durationInMs: PD.Numeric(1500, { min: 100, max: 15000, step: 100 }, { label: 'Show for (ms)' }),
};

const DefaultShowFor = EntryOptions.durationInMs.defaultValue;

/** The name the list and the card show: an unnamed snapshot is its position. */
const titleOf = (e: SnapshotInfo, index: number) => e.name || `Snapshot ${index + 1}`;

const clock = (ms: number) => new Date(ms).toLocaleTimeString([], { hour: 'numeric', minute: '2-digit' });

const toEntry = (e: SnapshotInfo, index: number): SnapshotEntry => ({
  id: e.id,
  name: titleOf(e, index),
  ...(e.key ? { key: e.key } : {}),
  ...(e.description ? { description: e.description } : {}),
  timestamp: clock(e.timestamp),
  ...(e.image ? { image: e.image } : {}),
});

const sameEntries = (a: readonly SnapshotEntry[], b: readonly SnapshotEntry[]) =>
  a.length === b.length &&
  a.every((x, i) => {
    const y = b[i];
    return (
      x.id === y.id &&
      x.name === y.name &&
      x.key === y.key &&
      x.description === y.description &&
      x.timestamp === y.timestamp &&
      x.image === y.image
    );
  });

export interface SnapshotsPanelModel extends SnapshotsModel {
  /**
   * Whether the tab shows: the entry forms follow the list only while it does (hidden UI holds no subscription;
   * AGENTS.md). Starts active.
   */
  setActive(active: boolean): void;
  dispose(): void;
}

/** What `session.save` gives (the export-session feature): the file's bytes, saved by the page. */
type Saved = { blob: Blob; filename: string };

export function createSnapshotsModel(
  api: ViewerApi,
  area: SnapshotApi,
  /** Hands a saved file to the page. */
  save: (blob: Blob, filename: string) => void,
): SnapshotsPanelModel {
  const state = area.state;
  const entries = mapState(state, (s) => s.entries.map(toEntry), sameEntries);
  const current = mapState(state, (s) => s.current, Object.is);
  const forms = new Map<string, ParamsFormModel<typeof EntryOptions>>();
  // A pending (debounced) edit of an entry's playback time, by entry.
  const timers = new Map<string, ReturnType<typeof setTimeout>>();
  const entryOf = (id: string) => state.value.entries.find((e) => e.id === id);
  const saveOptions = new ParamsFormModel(SaveOptions);

  const forget = (id: string) => {
    forms.delete(id);
    const t = timers.get(id);
    if (t) clearTimeout(t);
    timers.delete(id);
  };
  // One subscription for all the entry forms (not one per form), while the tab shows: it loads each form from the list
  // (a session opened, another edit) without reporting an edit, and drops the forms of removed entries.
  let follow: Subscription | undefined;
  let active = true;
  const startFollowing = () => {
    if (follow || !active || !forms.size) return;
    follow = state.subscribe((s) => {
      for (const [id, form] of [...forms]) {
        const e = s.entries.find((x) => x.id === id);
        // A pending edit is newer than the list: loading would pull the slider back mid-drag.
        if (e) {
          if (!timers.has(id)) form.load({ durationInMs: e.showFor ?? DefaultShowFor });
        } else forget(id);
      }
      if (!forms.size) stopFollowing();
    });
  };
  const stopFollowing = () => {
    follow?.unsubscribe();
    follow = undefined;
  };

  return {
    entries,
    current,
    saveOptions,
    entryOptions(id) {
      if (!entryOf(id)) return undefined;
      let form = forms.get(id);
      if (!form) {
        // A removed entry's form goes when the next one is made, even if the tab was hidden at the removal.
        for (const known of [...forms.keys()]) if (!entryOf(known)) forget(known);
        form = new ParamsFormModel(
          EntryOptions,
          { durationInMs: entryOf(id)?.showFor ?? DefaultShowFor },
          {
            onChange: (values) => {
              // A slider drags through many values: the edit is made when it rests.
              const pending = timers.get(id);
              if (pending) clearTimeout(pending);
              timers.set(
                id,
                setTimeout(() => {
                  timers.delete(id);
                  void api.snapshot?.update({ entry: id, showFor: values.durationInMs });
                }, 150),
              );
            },
          },
        );
        forms.set(id, form);
        startFollowing();
      }
      return form;
    },
    add: (name) => void api.snapshot?.add({ ...(name ? { name } : {}), options: { ...saveOptions.values.value } }),
    apply: (id) => void api.snapshot?.apply({ entry: id }),
    update: (id, patch) => {
      const e = entryOf(id);
      if (!e) return;
      // The operation keeps what is empty, unless `clear` says to empty it: the whole edit goes in.
      void api.snapshot?.update({
        entry: id,
        name: patch.name ?? e.name,
        description: patch.description ?? e.description ?? '',
        key: patch.key ?? e.key ?? '',
        clear: true,
      });
    },
    replace: (id) => void api.snapshot?.replace({ entry: id, options: { ...saveOptions.values.value } }),
    move: (id, by) => {
      const index = state.value.entries.findIndex((e) => e.id === id);
      const to = index + by;
      if (index < 0 || to < 0 || to >= state.value.entries.length) return;
      void api.snapshot?.move({ entry: id, position: to + 1 });
    },
    moveTo: (id, index) => {
      const count = state.value.entries.length;
      void api.snapshot?.move({ entry: id, position: Math.max(0, Math.min(count - 1, index)) + 1 });
    },
    remove: (id) => void api.snapshot?.remove({ entry: id }),
    clear: () => void api.snapshot?.clear(),
    download(kind) {
      // The export-session feature brings `session.save`; without it there is nothing to save.
      if (!api.describe('session.save')) return;
      void api.run('session.save', { kind }).then((r) => {
        const saved = r.ok ? (r.value as Saved | undefined) : undefined;
        if (saved) save(saved.blob, saved.filename);
      });
    },
    open(files) {
      if (files[0]) void api.session.open({ file: files[0] });
    },
    setActive(on) {
      active = on;
      if (on) startFollowing();
      else stopFollowing();
    },
    dispose() {
      active = false;
      stopFollowing();
      for (const t of timers.values()) clearTimeout(t);
      timers.clear();
      forms.clear();
    },
  };
}

/**
 * The playback card's model: the current snapshot, its place in the list, and stepping and playing. `markdown` is the
 * descriptions' host (Mol*'s markdown commands, `markdown.ts`), made by the views.
 */
export function createSnapshotStepper(
  api: ViewerApi,
  area: SnapshotApi,
  markdown?: MarkdownHost,
): SnapshotStepperModel {
  const state: ReadonlyState<SnapshotsState> = area.state;
  const indexOf = (s: SnapshotsState) => s.entries.findIndex((e) => e.id === s.current);
  const sameTitles = (a: readonly string[], b: readonly string[]) =>
    a.length === b.length && a.every((t, i) => t === b[i]);
  return {
    count: mapState(state, (s) => s.entries.length, Object.is),
    index: mapState(state, indexOf, Object.is),
    title: mapState(
      state,
      (s) => {
        const i = indexOf(s);
        return i < 0 ? undefined : titleOf(s.entries[i], i);
      },
      Object.is,
    ),
    titles: mapState(state, (s) => s.entries.map(titleOf), sameTitles),
    description: mapState(state, (s) => s.entries[indexOf(s)]?.description, Object.is),
    plainText: mapState(state, (s) => !!s.entries[indexOf(s)]?.plainText, Object.is),
    ...(markdown ? { markdown } : {}),
    isPlaying: mapState(state, (s) => s.playing, Object.is),
    previous: () => void api.snapshot?.previous(),
    next: () => void api.snapshot?.next(),
    goTo(index) {
      const e = state.value.entries[index];
      if (e) void api.snapshot?.apply({ entry: e.id });
    },
    togglePlay: () => void (state.value.playing ? api.snapshot?.stop() : api.snapshot?.play()),
    openLink(href) {
      // `#key` applies the snapshot with that key; a web link opens as one.
      if (href.startsWith('#')) {
        const e = state.value.entries.find((x) => x.key === href.slice(1));
        if (e) void api.snapshot?.apply({ entry: e.id });
      } else if (/^https?:\/\//i.test(href)) {
        globalThis.open?.(href, '_blank', 'noopener,noreferrer');
      }
    },
  };
}

/**
 * The current snapshot's transition as the canvas scrubber shows it (prism-ui's `SnapshotTransitionModel`), acting
 * through `snapshot.seek` and `snapshot.transition`. A drag reports every step: one seek at a time, the last wins.
 */
export function createSnapshotTransition(api: ViewerApi, area: SnapshotApi): SnapshotTransitionModel {
  const state = area.transitionState;
  let seeking = false;
  let wanted: number | undefined;
  const seek = async () => {
    if (seeking) return;
    seeking = true;
    try {
      while (wanted !== undefined) {
        const time = wanted;
        wanted = undefined;
        await api.snapshot?.seek({ time });
      }
    } finally {
      seeking = false;
    }
  };
  return {
    duration: mapState(state, (s) => s.duration, Object.is),
    step: mapState(state, (s) => s.step, Object.is),
    time: mapState(state, (s) => s.time, Object.is),
    isPlaying: mapState(state, (s) => s.playing, Object.is),
    seek(time) {
      wanted = time;
      void seek();
    },
    togglePlay: () => void api.snapshot?.transition({ play: !state.value.playing }),
  };
}

/** The playback card's part id. */
export const PlaybackCard = 'snapshots.overlay.playback';

/** Whether the card has something to show: a story to step through, or one snapshot's description. */
export const presentable = (s: SnapshotsState): boolean =>
  s.entries.length >= 2 || !!s.entries.find((e) => e.id === s.current)?.description?.trim();

/** The operations that load data or a session: a story comes with them (MVS, a session, a state file). */
const loads = (id: string) => id.startsWith('load.') || id === 'session.open' || id === 'session.remote.open';
const listOf = (s: SnapshotsState) => s.entries.map((e) => e.id).join(' ');

/**
 * M3 (parity.md §4, the playback row): a load that brings a story (two or more snapshots, or one with a description)
 * opens the playback card, once per load, as v5 shows the description and the stepper at once; one that finished before
 * the views came (their chunk loads lazily) opens it when they do. A load that leaves the
 * list as it was (another structure) opens nothing. Returns its disposer.
 */
export function openPlaybackOnLoad(
  api: Pick<ViewerApi, 'events' | 'view'> & Partial<Pick<ViewerApi, 'running'>>,
  area: Pick<SnapshotApi, 'state'>,
): () => void {
  const before = new Map<number, string>();
  // The views come in their own chunk, maybe while a load runs (a URL's, a story set up at once): what that load
  // brings is unknown (its snapshots may already be in), so when it ends a story opens the card.
  const loading = !!api.running?.value.some((r) => loads(r.operation));
  const mounted = loading ? '' : listOf(area.state.value);
  // A load inside another (files that hold an MVS) opens the card once.
  let opened = '';
  const started = api.events.started.subscribe((r) => {
    if (loads(r.operation)) before.set(r.id, listOf(area.state.value));
  });
  const finished = api.events.finished.subscribe((r) => {
    if (!loads(r.operation)) return;
    const was = before.get(r.id) ?? mounted;
    before.delete(r.id);
    const s = area.state.value;
    const now = listOf(s);
    if (!r.result?.ok || now === was || now === opened || !presentable(s)) return;
    opened = now;
    void api.view.open({ card: PlaybackCard });
  });
  // A load that brought a story before the views came (a page's API load, a story set up at once): with no load
  // running to wait for, it opens now, once.
  const s = area.state.value;
  if (!loading && mounted && presentable(s)) {
    opened = mounted;
    void api.view.open({ card: PlaybackCard });
  }
  return () => {
    started.unsubscribe();
    finished.unsubscribe();
  };
}
