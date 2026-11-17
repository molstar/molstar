/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Timed and durationless playback (playback §1): one interface, three kinds by what the model knows. A `duration`
// (seconds) is timed, a `count` of frames is counted: both seek, and the row has a scrubber. Neither is durationless:
// something that runs until stopped (a spin, a live simulation, a stream): play and pause, how long it ran or how many
// steps it took, one step while paused, restart, a speed; nothing to seek, no end. Small reference models for stories,
// tests and a headless host. DOM-free.

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { type PlaybackClock, TimerClock } from './clock.js';
import { DefaultPlaybackLabels, type PlaybackLabels } from './labels.js';

/**
 * What `TimelineRow` reads and does. A model has the fields its kind needs (`timelineKind`). It holds no words: the row's
 * name is the view's `label`, the rest its `labels` (playback §3).
 */
export interface TimelineModel {
  readonly playing: ReadonlyState<boolean>;
  togglePlay(): void;

  // Timed: a length in seconds and where it is.
  readonly duration?: ReadonlyState<number>;
  /** Seconds. */
  readonly time?: ReadonlyState<number>;
  /** The scrubber's step in seconds (default: any). */
  readonly step?: ReadonlyState<number>;

  // Counted: frames and the current one.
  readonly count?: ReadonlyState<number>;
  /** 0-based. */
  readonly index?: ReadonlyState<number>;

  /** To a time (timed) or a frame (counted). */
  seek?(at: number): void;
  /** Past the end goes to the start (playing, and a counted model's next and previous). */
  readonly periodic?: boolean;

  // Durationless: how long and how far it ran.
  /** Seconds it ran. */
  readonly elapsed?: ReadonlyState<number>;
  /** Steps it took. */
  readonly steps?: ReadonlyState<number>;
  /** Back to its start (it keeps playing or paused). */
  restart?(): void;
  /** One step while paused. */
  stepOnce?(): void;

  // Any kind: a speed, from a few choices.
  readonly speed?: ReadonlyState<number>;
  readonly speeds?: readonly number[];
  setSpeed?(speed: number): void;
}

export type TimelineKind = 'timed' | 'counted' | 'durationless';

/** A model's kind: by what it knows, a duration first. */
export function timelineKind(model: TimelineModel): TimelineKind {
  if (model.duration) return 'timed';
  if (model.count) return 'counted';
  return 'durationless';
}

/** Does the row show (a timed model with a length, a counted one with more than a frame, a durationless one always)? */
export function timelineShown(model: TimelineModel): boolean {
  const kind = timelineKind(model);
  if (kind === 'timed') return (model.duration?.value ?? 0) > 0;
  if (kind === 'counted') return (model.count?.value ?? 0) > 1;
  return true;
}

/** Seconds with a tenth ("1.2"): a timed model's time. */
export function formatSeconds(s: number): string {
  return (Math.max(0, s) || 0).toFixed(1);
}

/**
 * Where it is, shown in the row: "0.4 / 0.8 s", "3 / 10", "12 s · step 340" (`labels`: the words, English by default).
 */
export function timelinePosition(model: TimelineModel, labels: PlaybackLabels = DefaultPlaybackLabels): string {
  const kind = timelineKind(model);
  if (kind === 'timed')
    return labels.timedPosition(formatSeconds(model.time?.value ?? 0), formatSeconds(model.duration!.value));
  if (kind === 'counted') return `${(model.index?.value ?? 0) + 1} / ${model.count!.value}`;
  return ranParts(model, labels, ' · ');
}

/**
 * Where it is, said (a scrubber's or a durationless status's `aria-valuetext`): "1.2 s of 4.0 s", "frame 3 of 10";
 * durationless "running, 12 s", "paused at step 340" (with both: "running, 12 s, step 340").
 */
export function timelineValueText(model: TimelineModel, labels: PlaybackLabels = DefaultPlaybackLabels): string {
  const kind = timelineKind(model);
  if (kind === 'timed')
    return labels.timedValue(formatSeconds(model.time?.value ?? 0), formatSeconds(model.duration!.value));
  if (kind === 'counted') return labels.countedValue((model.index?.value ?? 0) + 1, model.count!.value);
  const ran = ranParts(model, labels, ', ');
  return model.playing.value ? labels.running(ran) : labels.paused(ran);
}

function ranParts(model: TimelineModel, labels: PlaybackLabels, separator: string): string {
  const elapsed = model.elapsed ? labels.elapsed(model.elapsed.value) : '';
  const steps = model.steps ? labels.stepCount(model.steps.value) : '';
  return elapsed && steps ? `${elapsed}${separator}${steps}` : elapsed || steps;
}

// --- Reference models ---------------------------------------------------------------------------------------------

/** The speeds the reference models offer unless told. */
export const DefaultSpeeds: readonly number[] = [0.25, 0.5, 1, 2, 4];

interface Common {
  readonly clock?: PlaybackClock;
  /** A speed choice: these (default none: no choice); `speed` the one to start at (default 1). */
  readonly speeds?: readonly number[];
  readonly speed?: number;
}

export interface TimelineOptions extends Common {
  /** Seconds: a timed model. */
  readonly duration?: number;
  /** Frames: a counted model (when no `duration`). */
  readonly count?: number;
  /** Ticks a second while playing (default 30 timed, 10 counted: a frame each). */
  readonly fps?: number;
  /** Past the end to the start (default); otherwise playing stops at the end. */
  readonly periodic?: boolean;
  /** Told each new time (timed) or frame (counted): the app shows it. */
  readonly onSeek?: (at: number) => void;
}

/**
 * A timed (`duration`) or counted (`count`) reference model: seeks, plays at `fps` times its speed, wraps past the end
 * (`periodic`) or stops there; played from the end, it starts over. Ticks only while playing.
 */
export class Timeline implements TimelineModel {
  readonly periodic: boolean;
  readonly duration: ReadonlyState<number> | undefined;
  readonly time: ReadonlyState<number> | undefined;
  readonly count: ReadonlyState<number> | undefined;
  readonly index: ReadonlyState<number> | undefined;
  readonly speeds: readonly number[] | undefined;
  readonly #options: TimelineOptions;
  readonly #at = new State(0);
  readonly #length: State<number>;
  readonly #timed: boolean;
  readonly #playing = new State(false);
  readonly #speed: State<number>;
  #stop: (() => void) | undefined;

  constructor(options: TimelineOptions) {
    this.#options = options;
    this.periodic = options.periodic !== false;
    this.#timed = options.duration !== undefined;
    this.#length = new State(this.#timed ? options.duration! : (options.count ?? 0));
    if (this.#timed) {
      this.duration = this.#length.asReadonly;
      this.time = this.#at.asReadonly;
    } else {
      this.count = this.#length.asReadonly;
      this.index = this.#at.asReadonly;
    }
    this.speeds = options.speeds;
    this.#speed = new State(options.speed ?? 1);
  }

  get playing(): ReadonlyState<boolean> {
    return this.#playing.asReadonly;
  }
  get speed(): ReadonlyState<number> | undefined {
    return this.speeds ? this.#speed.asReadonly : undefined;
  }

  /** A time (s) or a frame, held to the length; not a number: nothing. */
  seek(at: number): void {
    if (Number.isNaN(at)) return;
    const max = this.#timed ? this.#length.value : this.#length.value - 1;
    const v = Math.min(Math.max(this.#timed ? at : Math.round(at), 0), Math.max(max, 0));
    if (v === this.#at.value) return;
    this.#at.set(v);
    this.#options.onSeek?.(v);
  }

  /**
   * Changes the length (more frames loaded, a shorter animation); where it is stays where it can. Nothing left to play
   * (no length, one frame) stops it.
   */
  setLength(length: number): void {
    this.#length.set(Math.max(0, length) || 0);
    this.seek(this.#at.value);
    if (this.#end() <= 0) this.pause();
  }

  play(): void {
    if (this.#playing.value) return;
    const end = this.#end();
    if (end <= 0) return;
    if (this.#at.value >= end && !this.periodic) this.seek(0);
    this.#playing.set(true);
    this.#start();
  }

  pause(): void {
    this.#stop?.();
    this.#stop = undefined;
    this.#playing.set(false);
  }

  togglePlay(): void {
    if (this.#playing.value) this.pause();
    else this.play();
  }

  setSpeed(speed: number): void {
    if (!(speed > 0) || speed === this.#speed.value) return;
    this.#speed.set(speed);
    if (this.#stop) this.#start();
  }

  dispose(): void {
    this.pause();
  }

  /** The last time or frame: 0 or less, nothing to play. */
  #end(): number {
    return this.#timed ? this.#length.value : this.#length.value - 1;
  }

  #fps(): number {
    return this.#options.fps ?? (this.#timed ? 30 : 10);
  }

  #start(): void {
    this.#stop?.();
    const clock = this.#options.clock ?? TimerClock;
    // Timed: ticks at `fps`, each a frame's time times the speed. Counted: a frame a tick, ticks `fps` times the speed.
    const ms = this.#timed ? 1000 / this.#fps() : 1000 / (this.#fps() * this.#speed.value);
    this.#stop = clock.every(ms, () => this.#tick());
  }

  #tick(): void {
    // Shortened under it to nothing to play (`setLength` stops it; an app's own length state may not).
    if (this.#end() <= 0) return this.pause();
    const length = this.#length.value;
    if (this.#timed) {
      const t = this.#at.value + this.#speed.value / this.#fps();
      if (t < length) return this.seek(t);
      if (this.periodic) return this.seek(t - length * Math.floor(t / length));
      this.seek(length);
      return this.pause();
    }
    const i = this.#at.value + 1;
    if (i < length) return this.seek(i);
    if (this.periodic) return this.seek(0);
    this.pause();
  }
}

export interface DurationlessOptions extends Common {
  /** Ms a step takes at speed 1 (default 100). */
  readonly stepMs?: number;
  /** Shows the time it ran (default true). */
  readonly elapsed?: boolean;
  /** Shows the steps it took (default true). */
  readonly steps?: boolean;
  /** Offers Restart (default true). */
  readonly restart?: boolean;
  /** Offers one step while paused (default true). */
  readonly stepOnce?: boolean;
  /** Told each step (the app does it: turns the camera, advances the simulation). */
  readonly onStep?: (steps: number) => void;
  /** Told a restart. */
  readonly onRestart?: () => void;
}

/**
 * A durationless reference model: playing takes a step every `stepMs` over its speed, so a step is that long of run
 * time (`elapsed`: half as long at 2×); paused, one step at a time; restart goes back to none. Ticks only while playing.
 */
export class Durationless implements TimelineModel {
  readonly elapsed: ReadonlyState<number> | undefined;
  readonly steps: ReadonlyState<number> | undefined;
  readonly speeds: readonly number[] | undefined;
  readonly restart: (() => void) | undefined;
  readonly stepOnce: (() => void) | undefined;
  readonly #options: DurationlessOptions;
  readonly #elapsed = new State(0);
  readonly #steps = new State(0);
  readonly #playing = new State(false);
  readonly #speed: State<number>;
  /** The run time and the steps at the last speed change: `elapsed` counts on from there (playback §2). */
  #baseElapsed = 0;
  #baseSteps = 0;
  #stop: (() => void) | undefined;

  constructor(options: DurationlessOptions = {}) {
    this.#options = options;
    this.elapsed = options.elapsed === false ? undefined : this.#elapsed.asReadonly;
    this.steps = options.steps === false ? undefined : this.#steps.asReadonly;
    this.speeds = options.speeds;
    this.#speed = new State(options.speed ?? 1);
    this.restart = options.restart === false ? undefined : () => this.#restart();
    this.stepOnce = options.stepOnce === false ? undefined : () => this.#stepOnce();
  }

  get playing(): ReadonlyState<boolean> {
    return this.#playing.asReadonly;
  }
  get speed(): ReadonlyState<number> | undefined {
    return this.speeds ? this.#speed.asReadonly : undefined;
  }

  play(): void {
    if (this.#playing.value) return;
    this.#playing.set(true);
    this.#start();
  }

  pause(): void {
    this.#stop?.();
    this.#stop = undefined;
    this.#playing.set(false);
  }

  togglePlay(): void {
    if (this.#playing.value) this.pause();
    else this.play();
  }

  setSpeed(speed: number): void {
    if (!(speed > 0) || speed === this.#speed.value) return;
    this.#baseElapsed = this.#elapsed.value;
    this.#baseSteps = this.#steps.value;
    this.#speed.set(speed);
    if (this.#stop) this.#start();
  }

  dispose(): void {
    this.pause();
  }

  #restart(): void {
    this.#baseElapsed = 0;
    this.#baseSteps = 0;
    this.#steps.set(0);
    this.#elapsed.set(0);
    this.#options.onRestart?.();
  }

  #stepOnce(): void {
    if (!this.#playing.value) this.#step();
  }

  #step(): void {
    const n = this.#steps.value + 1;
    this.#steps.set(n);
    // A step is `stepMs` over the speed it ran at: counted on from the last speed change, not summed (no drift).
    const stepS = (this.#options.stepMs ?? 100) / 1000 / this.#speed.value;
    this.#elapsed.set(this.#baseElapsed + (n - this.#baseSteps) * stepS);
    this.#options.onStep?.(n);
  }

  #start(): void {
    this.#stop?.();
    const clock = this.#options.clock ?? TimerClock;
    this.#stop = clock.every((this.#options.stepMs ?? 100) / this.#speed.value, () => this.#step());
  }
}
