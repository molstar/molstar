/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rafScheduler, type Scheduler } from '@molstar/seq-viewer';
import { DefaultSequenceTheme } from './document.js';
import type { SequenceEdits, SequenceOverride, SequenceStyle, SequenceViewer } from './model.js';
import { DefaultTracksRail, type TracksRail } from './presets.js';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/prism-ui/reactive/frame-state';

// The sequence's runtime (sequence.md §2, viewer-spec.md §2.3): Seq* and everything over it, one chunk the spec loads
// (`sequence.runtime`) when a structure is first shown, or an operation or a saved state needs it. The viewer's core
// holds only what the plugin side needs (`sources/sequence.ts`), so a page shows its first frame without it, and a spec
// without a sequence has none of it in its bundle.

export { fromJSON, readDocument, readSvsx } from '@molstar/seq-core/format';
export { svsStructureSource } from './document.js';
export { createSequenceViewers } from './model.js';
export { sequenceShellUI } from './shell-ui.js';

/** Frames for the sequence view: the page's, or timers headless. */
export function sequenceFrames(page: boolean): Scheduler {
  return rafScheduler({
    requestAnimationFrame: page
      ? (cb) => requestFrameCallback(cb)
      : (cb) => setTimeout(() => cb(performance.now()), 16) as unknown as number,
    cancelAnimationFrame: page ? (id) => cancelFrameCallback(id) : (id) => clearTimeout(id),
    setTimeout: (cb, ms) => setTimeout(cb, ms) as unknown as number,
    clearTimeout: (id) => clearTimeout(id),
    now: () => performance.now(),
  });
}

/**
 * The sequence's state in snapshots (sequence.md §7.1): version 2 adds the loaded SeqViewSpec, version 3 the person's
 * edits over the document (hidden tracks, by ref), version 4 makes them each structure's (tracks put in, taken out or
 * hidden, and the chain filter, by structure ref).
 */
export interface SavedSequence {
  readonly version: 1 | 2 | 3 | 4;
  readonly style: SequenceStyle;
  readonly colors: string | null;
  readonly rail: TracksRail;
  readonly structure?: string;
  readonly override?: SequenceOverride;
  readonly edits?: SequenceEdits;
}

/** Version 3's edits: one set of hidden tracks over whatever document was shown. */
interface EditsV3 {
  readonly tracks?: Readonly<Record<string, { readonly hidden: boolean }>>;
}

/**
 * The viewer's state as saved (`prism-viewer.sequence`); a structure that is gone (`exists` false) leaves its edits
 * out.
 */
export function savedOf(
  viewer: SequenceViewer,
  structure: string | undefined,
  exists: (ref: string) => boolean = () => true,
): SavedSequence {
  const kept = Object.entries(viewer.edits.value).filter(([ref]) => exists(ref));
  return {
    version: 4,
    style: viewer.style.value,
    colors: viewer.theme.value.letters,
    rail: viewer.rail.value,
    structure,
    ...(viewer.override.value ? { override: viewer.override.value } : {}),
    ...(kept.length ? { edits: Object.fromEntries(kept) } : {}),
  };
}

/** A saved state's edits as version 4 has them: version 3's are the structure's it names, else none. */
function editsOf(v: SavedSequence): SequenceEdits {
  if (v.version === 4) return v.edits ?? {};
  if (v.version !== 3) return {};
  const old = v.edits as EditsV3 | undefined;
  if (!old?.tracks || !v.structure || !Object.keys(old.tracks).length) return {};
  return { [v.structure]: { tracks: old.tracks } };
}

/**
 * A saved state applied (a snapshot, a session, Undo): none, the defaults (the spec's initial style). What the sequence
 * follows is the source's (`setFollowed`).
 */
export function applySaved(
  viewer: SequenceViewer,
  follow: (ref: string | undefined) => void,
  value: unknown,
  initial: 'tracks' | 'text' | undefined,
): void {
  // Version 1 is version 2 without a loaded SeqViewSpec, version 2 is version 3 with no edits, version 3 is version 4
  // with its edits the followed structure's (the migrations).
  const v = value as SavedSequence | undefined;
  if (!v || (v.version !== 1 && v.version !== 2 && v.version !== 3 && v.version !== 4)) {
    viewer.setStyle(initial === 'text' ? 'compact' : 'tracks');
    viewer.setTheme(DefaultSequenceTheme);
    viewer.setRail(DefaultTracksRail);
    follow(undefined);
    void viewer.setEdits({});
    void viewer.setOverride(undefined);
    return;
  }
  viewer.setStyle(v.style);
  viewer.setTheme({ letters: v.colors });
  if (v.rail) viewer.setRail(v.rail);
  // Undo's checkpoints leave out what the sequence follows (view state); a session has it.
  if ('structure' in v) follow(v.structure);
  // The edits before the override: the document it rebuilds is written with them.
  void viewer.setEdits(editsOf(v));
  void viewer.setOverride(v.override);
}
