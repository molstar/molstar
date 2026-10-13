/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AnimateStateSnapshotTransition } from '../animation/built-in/state-snapshots.js';
import type { MarkdownExtension } from '../manager/markdown-extensions.js';

export const ApplySnapshotMarkdownExtension: MarkdownExtension = {
  name: 'apply-snapshot',
  execute: ({ event, args, manager }) => {
    if (event !== 'click') return;
    const key = args['apply-snapshot'];
    if (!key) return;
    manager.plugin.managers.snapshot.applyKey(key);
  },
};

export const NextSnapshotMarkdownExtension: MarkdownExtension = {
  name: 'next-snapshot',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('next-snapshot' in args)) return;
    let dir: -1 | 1 = (+args['next-snapshot'] || 1) as -1 | 1;
    if (!dir) return;
    if (dir < 0) dir = -1;
    else dir = 1;
    manager.plugin.managers.snapshot.applyNext(dir);
  },
};

export const PlayTransitionMarkdownExtension: MarkdownExtension = {
  name: 'play-transition',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('play-transition' in args)) return;
    manager.plugin.managers.animation.play(AnimateStateSnapshotTransition, {});
  },
};

export const PlaySnapshotsMarkdownExtension: MarkdownExtension = {
  name: 'play-snapshots',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('play-snapshots' in args)) return;
    manager.plugin.managers.snapshot.play({ restart: true });
  },
};

export const StopAnimationMarkdownExtension: MarkdownExtension = {
  name: 'stop-animation',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('stop-animation' in args)) return;
    manager.plugin.managers.snapshot.stop();
  },
};
