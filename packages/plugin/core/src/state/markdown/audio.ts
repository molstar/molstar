/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MarkdownExtension } from '../manager/markdown-extensions.js';

export const PlayAudioMarkdownExtension: MarkdownExtension = {
  name: 'play-audio',
  execute: ({ event, args, manager }) => {
    if (event !== 'click') return;

    const src = args['play-audio'];
    if (!src?.length) return;
    manager.audio.play(src);
  },
};

export const ToggleAudioMarkdownExtension: MarkdownExtension = {
  name: 'toggle-audio',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('toggle-audio' in args)) return;

    const src = args['toggle-audio'];
    manager.audio.play(src, { toggle: true });
  },
};

export const PauseAudioMarkdownExtension: MarkdownExtension = {
  name: 'pause-audio',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('pause-audio' in args)) return;
    manager.audio.pause();
  },
};

export const StopAudioMarkdownExtension: MarkdownExtension = {
  name: 'stop-audio',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('stop-audio' in args)) return;
    manager.audio.stop();
  },
};

export const DisposeAudioMarkdownExtension: MarkdownExtension = {
  name: 'dispose-audio',
  execute: ({ event, args, manager }) => {
    if (event !== 'click' || !('dispose-audio' in args)) return;
    manager.audio.dispose();
  },
};
