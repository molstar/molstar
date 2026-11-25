/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The playback views' text (playback §3): every word they show or say, English by default. DOM-free: what a row says is
// a pure function of the model and these (`timelineValueText`), so a headless host says the same. What one frame or
// step is called ("model", "scene") is the app's: it sets the words that name it.

/** The playback views' own text (`DefaultPlaybackLabels` English); parametrized ones are functions. */
export interface PlaybackLabels {
  /** The deck's name. */
  deck: string;
  /** A stepper's name when its row or card is given none (`label`). */
  steps: string;
  /** The position button's name when past `MaxDots`: the row's name and the current step's title (if any). */
  stepPosition: (name: string, title: string | undefined) => string;
  /** The position button's tooltip when no step is current. */
  noStep: (name: string) => string;
  /** The list's filter: its placeholder. */
  filter: string;
  /** The list's filter: its accessible name. */
  filterSteps: (name: string) => string;
  previous: string;
  next: string;
  play: string;
  pause: string;
  /** A timed row's name when given none. */
  timeline: string;
  /** A counted row's name when given none. */
  frames: string;
  /** A durationless row's name when given none. */
  durationless: string;
  /** A timed scrubber's name. */
  time: string;
  /** A counted scrubber's name. */
  frame: string;
  previousFrame: string;
  nextFrame: string;
  /** A timed row's position as shown, from formatted seconds: "0.4 / 2.0 s". */
  timedPosition: (time: string, duration: string) => string;
  /** A timed scrubber's value as said: "1.2 s of 4.0 s". */
  timedValue: (time: string, duration: string) => string;
  /** A counted scrubber's value as said, 1-based: "frame 3 of 10". */
  countedValue: (frame: number, count: number) => string;
  /** How long a durationless model ran (seconds, whole): "12 s", "1:05". */
  elapsed: (seconds: number) => string;
  /** How far a durationless model ran: "step 340". */
  stepCount: (steps: number) => string;
  /** A durationless status's name, after its row. */
  status: (name: string) => string;
  /** A durationless status's value while running, from what it ran ("12 s, step 340", or empty). */
  running: (ran: string) => string;
  /** A durationless status's value while paused, from what it ran. */
  paused: (ran: string) => string;
  /** One step, while paused. */
  stepOnce: string;
  /** One step's tooltip while it is disabled (playing). */
  stepOncePlaying: string;
  restart: string;
  /** The speed choice's name. */
  speed: string;
  /** A speed as a choice says it: "1×". */
  speedValue: (speed: number) => string;
  mute: string;
}

/** Time run, whole: "12 s", "1:05", "1:02:05". */
export function formatElapsed(s: number): string {
  const t = Math.floor(Math.max(0, s) || 0);
  if (t < 60) return `${t} s`;
  const h = Math.floor(t / 3600);
  const m = Math.floor((t % 3600) / 60);
  const sec = String(t % 60).padStart(2, '0');
  return h ? `${h}:${String(m).padStart(2, '0')}:${sec}` : `${m}:${sec}`;
}

export const DefaultPlaybackLabels: PlaybackLabels = {
  deck: 'Playback',
  steps: 'Steps',
  stepPosition: (name, title) => (title ? `${name}: ${title}` : name),
  noStep: (name) => `No ${name.toLowerCase()}`,
  filter: 'Filter',
  filterSteps: (name) => `Filter ${name.toLowerCase()}`,
  previous: 'Previous',
  next: 'Next',
  play: 'Play',
  pause: 'Pause',
  timeline: 'Timeline',
  frames: 'Frames',
  durationless: 'Playback',
  time: 'Time',
  frame: 'Frame',
  previousFrame: 'Previous frame',
  nextFrame: 'Next frame',
  timedPosition: (time, duration) => `${time} / ${duration} s`,
  timedValue: (time, duration) => `${time} s of ${duration} s`,
  countedValue: (frame, count) => `frame ${frame} of ${count}`,
  elapsed: formatElapsed,
  stepCount: (steps) => `step ${steps}`,
  status: (name) => `${name} status`,
  running: (ran) => (ran ? `running, ${ran}` : 'running'),
  paused: (ran) => (ran ? `paused at ${ran}` : 'paused'),
  stepOnce: 'One step',
  stepOncePlaying: 'One step (pause first)',
  restart: 'Restart',
  speed: 'Speed',
  speedValue: (speed) => `${speed}×`,
  mute: 'Mute',
};
