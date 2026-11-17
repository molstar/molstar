/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Stepped playback (playback §1): a count of steps with titles, the current one, previous, next, go to, play. The
// interface an app implements over what it steps through, and a small reference model (stories, tests, a headless
// host). DOM-free.

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import type { MarkdownHost } from '../markdown/markdown.js';
import { type PlaybackClock, TimerClock } from './clock.js';

/** What the stepper views read (`StepperRow`, `StepCard`) and do. */
export interface StepperModel {
  readonly count: ReadonlyState<number>;
  /** 0-based; -1 when no step is current. */
  readonly index: ReadonlyState<number>;
  /** Every step's title, in order (the dots' names, the list). */
  readonly titles: ReadonlyState<readonly string[]>;
  /** The current step's title. */
  readonly title: ReadonlyState<string | undefined>;
  /** The current step's description: markdown, unless `plainText`. */
  readonly description?: ReadonlyState<string | undefined>;
  /** The description is plain text, shown as it is (its line breaks kept). */
  readonly plainText?: ReadonlyState<boolean>;
  /** The page's part of the markdown (links it runs, images it resolves, media). */
  readonly markdown?: MarkdownHost;
  /** Steps on by itself (a slideshow). */
  readonly playing: ReadonlyState<boolean>;
  previous(): void;
  next(): void;
  goTo(index: number): void;
  togglePlay(): void;
  /** Follows a link in the description that the markdown does not (`#key`, an app's own). */
  openLink?(href: string): void;
}

/** A step of the reference model. */
export interface Step {
  readonly title: string;
  readonly description?: string;
}

export interface StepperOptions {
  readonly steps?: readonly Step[];
  /** The step to start at (default 0; -1 for none). */
  readonly index?: number;
  /** Ms a step shows while playing (default 4000). */
  readonly intervalMs?: number;
  /** Playing past the last step goes to the first (default); otherwise it stops there. */
  readonly loop?: boolean;
  readonly clock?: PlaybackClock;
  readonly plainText?: boolean;
  readonly markdown?: MarkdownHost;
  /** A description's link: `#n` goes to step n (1-based) when not given. */
  readonly onLink?: (href: string) => void;
}

/**
 * A stepper over a list of steps: previous and next wrap around; playing steps on every `intervalMs` (from two steps
 * on), a move by hand starting its interval again. The ticks are the model's (an injected clock), only while playing.
 */
export class Stepper implements StepperModel {
  readonly #steps: State<readonly Step[]>;
  readonly #index: State<number>;
  readonly #playing = new State(false);
  readonly #options: StepperOptions;
  #stop: (() => void) | undefined;

  readonly count: ReadonlyState<number>;
  readonly titles: ReadonlyState<readonly string[]>;
  readonly title: ReadonlyState<string | undefined>;
  readonly description: ReadonlyState<string | undefined>;
  readonly plainText: ReadonlyState<boolean>;
  readonly markdown: MarkdownHost | undefined;

  constructor(options: StepperOptions = {}) {
    this.#options = options;
    const steps = options.steps ?? [];
    this.#steps = new State<readonly Step[]>(steps);
    this.#index = new State(steps.length ? Math.min(options.index ?? 0, steps.length - 1) : -1);
    this.count = derive(this.#steps, (s) => s.length);
    this.titles = derive(this.#steps, (s) => s.map((x) => x.title));
    this.title = derive([this.#steps, this.#index], (s, i) => s[i]?.title);
    this.description = derive([this.#steps, this.#index], (s, i) => s[i]?.description);
    this.plainText = new State(!!options.plainText);
    this.markdown = options.markdown;
  }

  get index(): ReadonlyState<number> {
    return this.#index.asReadonly;
  }
  get playing(): ReadonlyState<boolean> {
    return this.#playing.asReadonly;
  }

  /** Replaces the steps; the current one stays where it can. */
  setSteps(steps: readonly Step[]): void {
    this.#steps.set(steps);
    const i = this.#index.value;
    this.#index.set(steps.length ? Math.min(Math.max(i, 0), steps.length - 1) : -1);
    if (steps.length < 2) this.pause();
  }

  goTo(index: number): void {
    const n = this.#steps.value.length;
    if (!n) return;
    this.#index.set(Math.min(Math.max(Math.trunc(index), 0), n - 1));
    if (this.#stop) this.#start();
  }

  next(): void {
    const n = this.#steps.value.length;
    if (n) this.goTo((this.#index.value + 1) % n);
  }

  previous(): void {
    const n = this.#steps.value.length;
    if (n) this.goTo((this.#index.value - 1 + n) % n);
  }

  play(): void {
    if (this.#playing.value || this.#steps.value.length < 2) return;
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

  openLink(href: string): void {
    if (this.#options.onLink) return this.#options.onLink(href);
    const m = /^#(\d+)$/.exec(href);
    if (m) this.goTo(Number(m[1]) - 1);
  }

  /** Stops its ticks. */
  dispose(): void {
    this.pause();
  }

  #start(): void {
    this.#stop?.();
    const clock = this.#options.clock ?? TimerClock;
    this.#stop = clock.every(this.#options.intervalMs ?? 4000, () => this.#tick());
  }

  #tick(): void {
    const n = this.#steps.value.length;
    const i = this.#index.value;
    if (i + 1 < n) this.#index.set(i + 1);
    else if (this.#options.loop !== false && n > 1) this.#index.set(0);
    else this.pause();
  }
}

/**
 * A description without its first heading when that says nothing the step's title (shown above it) does not: the same
 * words, or fewer of them ("ATP Binding" under "ATP Binding [1/2]"), case, markdown marks, punctuation and a part
 * counter ("[1/2]", "(2/5)") aside. A heading that says more stays.
 */
export function withoutTitle(markdown: string, title: string): string {
  const m = /^\s*#{1,6}[ \t]+(.+?)[ \t]*#*[ \t]*(?:\r?\n|$)/.exec(markdown);
  if (!m) return markdown;
  const heading = words(m[1]);
  const said = words(title);
  if (
    !heading ||
    !(
      said === heading ||
      said.startsWith(`${heading} `) ||
      said.includes(` ${heading} `) ||
      said.endsWith(` ${heading}`)
    )
  )
    return markdown;
  return markdown.slice(m[0].length).replace(/^\s+/, '');
}

/** A heading's or a title's words: lower case, no markdown marks, punctuation or part counter. */
function words(text: string): string {
  return text
    .replace(/[[(]\s*\d+\s*\/\s*\d+\s*[\])]/g, ' ')
    .replace(/[*_`~]/g, '')
    .toLowerCase()
    .replace(/[^\p{L}\p{N}]+/gu, ' ')
    .trim();
}
