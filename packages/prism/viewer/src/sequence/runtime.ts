/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rafScheduler, type Scheduler } from '@molstar/seq-viewer';
import { DefaultSequenceTheme } from './document.js';
import type { SequenceOverride, SequenceStyle, SequenceViewer } from './model.js';
import type { TracksRail } from './presets.js';

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
      ? (cb) => requestAnimationFrame(cb)
      : (cb) => setTimeout(() => cb(performance.now()), 16) as unknown as number,
    cancelAnimationFrame: page ? (id) => cancelAnimationFrame(id) : (id) => clearTimeout(id),
    setTimeout: (cb, ms) => setTimeout(cb, ms) as unknown as number,
    clearTimeout: (id) => clearTimeout(id),
    now: () => performance.now(),
  });
}

/** The sequence's state in snapshots (sequence.md §7.1): version 2 adds the loaded SeqViewSpec. */
export interface SavedSequence {
  readonly version: 1 | 2;
  readonly style: SequenceStyle;
  readonly colors: string | null;
  readonly rail: TracksRail;
  readonly structure?: string;
  readonly override?: SequenceOverride;
}

/** The viewer's state as saved (`prism-viewer.sequence`). */
export function savedOf(viewer: SequenceViewer, structure: string | undefined): SavedSequence {
  return {
    version: 2,
    style: viewer.style.value,
    colors: viewer.theme.value.letters,
    rail: viewer.rail.value,
    structure,
    ...(viewer.override.value ? { override: viewer.override.value } : {}),
  };
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
  // Version 1 is version 2 without a loaded SeqViewSpec (the migration).
  const v = value as SavedSequence | undefined;
  if (!v || (v.version !== 1 && v.version !== 2)) {
    viewer.setStyle(initial === 'text' ? 'compact' : 'tracks');
    viewer.setTheme(DefaultSequenceTheme);
    follow(undefined);
    void viewer.setOverride(undefined);
    return;
  }
  viewer.setStyle(v.style);
  viewer.setTheme({ letters: v.colors });
  if (v.rail) viewer.setRail(v.rail);
  // Undo's checkpoints leave out what the sequence follows (view state); a session has it.
  if ('structure' in v) follow(v.structure);
  void viewer.setOverride(v.override);
}
