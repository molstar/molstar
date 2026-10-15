/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ObservableState, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { MarkdownSource } from './source.js';

// A story's audio (Mol*'s markdown audio commands: play, toggle, pause, stop; MVS stories narrate with them) as the views
// read it: whether there is a player, whether it plays, where it is and how long it is; and its controls. The player is
// the plugin's one element (`markdownExtensions.state.audioPlayer`); this model follows its media events. Playback is
// not part of the scene: no history step.

export interface AudioModel {
  /** A player exists (a command made one, and it was not closed). */
  readonly present: ReadonlyState<boolean>;
  readonly playing: ReadonlyState<boolean>;
  /** Seconds played. */
  readonly time: ReadonlyState<number>;
  /** Seconds in all; 0 while unknown. */
  readonly duration: ReadonlyState<number>;
  /** Plays or pauses. */
  toggle(): void;
  /** Back to the start, paused. */
  stop(): void;
  /** To a fraction (0 to 1) of the length. */
  seek(fraction: number): void;
  /** Stops and lets the player go. */
  close(): void;
  dispose(): void;
}

export function createAudioModel(source: Pick<MarkdownSource, 'audio' | 'stopAudio' | 'closeAudio'>): AudioModel {
  const present = observableState(false);
  const playing = observableState(false);
  const time = observableState(0);
  const duration = observableState(0);
  let element: HTMLAudioElement | null = null;
  let detach = () => {};
  const set = <T>(state: ObservableState<T>, value: T) => {
    if (state.value !== value) state.set(value);
  };

  const attach = (audio: HTMLAudioElement | null) => {
    detach();
    element = audio;
    set(present, !!audio);
    if (!audio) {
      set(playing, false);
      set(time, 0);
      set(duration, 0);
      return;
    }
    const read = () => {
      set(playing, !audio.paused && !audio.ended);
      set(time, audio.currentTime || 0);
      set(duration, Number.isFinite(audio.duration) ? audio.duration : 0);
    };
    const events = ['play', 'pause', 'ended', 'timeupdate', 'durationchange', 'loadedmetadata', 'seeked', 'emptied'];
    for (const e of events) audio.addEventListener(e, read);
    detach = () => {
      for (const e of events) audio.removeEventListener(e, read);
    };
    read();
  };
  const sub = source.audio.subscribe(attach);

  return {
    present: present.asReadonly,
    playing: playing.asReadonly,
    time: time.asReadonly,
    duration: duration.asReadonly,
    toggle() {
      if (!element) return;
      if (element.paused || element.ended) void element.play().catch(() => {});
      else element.pause();
    },
    stop: () => source.stopAudio(),
    seek(fraction) {
      if (!element || !(element.duration > 0)) return;
      element.currentTime = Math.max(0, Math.min(1, fraction)) * element.duration;
    },
    close: () => source.closeAudio(),
    dispose() {
      sub.unsubscribe();
      detach();
    },
  };
}

/** "1:05": seconds as minutes and seconds. */
export function formatTime(seconds: number): string {
  const s = Math.max(0, Math.floor(seconds));
  return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, '0')}`;
}
