/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { describe, expect, it, vi } from 'vitest';
import { createAudioModel, formatTime } from '../audio.js';

// A story's audio as the views read it (features/snapshots/audio.ts): the plugin's player followed through its media
// events, and its controls.

/** An audio element whose playback is set by hand (jsdom plays nothing). */
function player() {
  const audio = document.createElement('audio');
  let paused = true;
  Object.defineProperty(audio, 'paused', { get: () => paused });
  Object.defineProperty(audio, 'duration', { value: 120, configurable: true });
  audio.play = vi.fn(async () => {
    paused = false;
    audio.dispatchEvent(new Event('play'));
  });
  audio.pause = vi.fn(() => {
    paused = true;
    audio.dispatchEvent(new Event('pause'));
  });
  return audio;
}

describe('the audio model', () => {
  it('follows the player: present, playing, time and length; toggles, seeks, stops and closes', async () => {
    const element = observableState<HTMLAudioElement | null>(null);
    const stopAudio = vi.fn();
    const closeAudio = vi.fn(() => element.set(null));
    const model = createAudioModel({ audio: element.asReadonly, stopAudio, closeAudio });
    expect(model.present.value).toBe(false);

    const audio = player();
    element.set(audio);
    expect(model.present.value).toBe(true);
    expect(model.duration.value).toBe(120);
    model.toggle();
    await Promise.resolve();
    expect(model.playing.value).toBe(true);
    model.seek(0.5);
    audio.dispatchEvent(new Event('timeupdate'));
    expect(model.time.value).toBe(60);
    model.toggle();
    expect(model.playing.value).toBe(false);
    model.stop();
    expect(stopAudio).toHaveBeenCalled();
    model.close();
    expect(model.present.value).toBe(false);
    model.dispose();
  });

  it('writes times as minutes and seconds', () => {
    expect(formatTime(0)).toBe('0:00');
    expect(formatTime(65.4)).toBe('1:05');
  });
});
