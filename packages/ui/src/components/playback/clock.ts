/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ticks of the reference playback models (playback §2): injected, so a test drives them with fake timers and a
// headless host with its own loop. DOM-free.

/** Calls `tick` every `ms` until the returned stop. */
export interface PlaybackClock {
  every(ms: number, tick: () => void): () => void;
}

/** The host's timers (`setInterval`, read at each call: fake timers drive it). */
export const TimerClock: PlaybackClock = {
  every(ms, tick) {
    const id = setInterval(tick, ms);
    return () => clearInterval(id);
  },
};
