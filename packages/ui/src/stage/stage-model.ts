/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { type CanvasInk, inkForColor } from './canvas-ink.js';

// The stage's state (stage §1): DOM-free, on core's states. The views (`Stage`, `EmptyStage`, `HoverLabel`) read it
// and write what they measure (`setTopInset`, `setSize`); an app, a script or an agent drives the same methods with no
// document. Its one wait (the empty mark's) runs on injectable timers.

/** A run of a formatted hover line: text only (shown by `textContent`, never as markup). */
export interface HoverRun {
  readonly text: string;
  readonly bold?: boolean;
  readonly italic?: boolean;
  readonly small?: boolean;
  readonly code?: boolean;
  readonly sub?: boolean;
  readonly sup?: boolean;
  /** Link text: shown as a link, not followed (a label takes no pointer). */
  readonly link?: boolean;
}

/** A formatted line: its runs and its role (a title, text, a list item). */
export interface HoverRichLine {
  readonly runs: readonly HoverRun[];
  readonly role: 'title' | 'text' | 'item';
}

/** A hover label's line: plain text (the first is the title), or formatted runs. */
export type HoverLine = string | HoverRichLine;

/** What a hover label shows: its lines, empty when nothing is hovered. */
export interface HoverLabelModel {
  readonly lines: ReadonlyState<readonly HoverLine[]>;
}

/**
 * Where the hover label shows: the stage's bottom right corner, above the toolbar; beside the pointer; or `auto` (the
 * default): the corner on a roomy stage, beside the pointer on a narrow or short one (a small pane, a window), where a
 * corner label would sit on the canvas's middle or be squeezed.
 */
export type HoverPlacement = 'auto' | 'corner' | 'cursor';

/** Where the label is now (`auto` decided by the stage's size). */
export type HoverAt = 'corner' | 'cursor';

/**
 * The empty mark's phase: `gone`, `in` (the stage has been empty and idle a while: the wordmark is up), `out` (fading
 * out as something came, only after it showed).
 */
export type StageMarkPhase = 'gone' | 'in' | 'out';

/** The stage's size in px; 0 while unknown (headless, not yet measured). */
export interface StageSize {
  readonly width: number;
  readonly height: number;
}

/** The timers the empty mark's wait runs on (the globals by default, read at each call: fake timers drive them). */
export interface StageTimers {
  setTimeout(fn: () => void, ms: number): unknown;
  clearTimeout(id: unknown): void;
}

const GlobalTimers: StageTimers = {
  setTimeout: (fn, ms) => setTimeout(fn, ms),
  clearTimeout: (id) => clearTimeout(id as ReturnType<typeof setTimeout>),
};

/**
 * How long the stage stays empty and idle before the mark shows (ms): clearing a scene to load the next never flashes
 * it, and work restarts the wait (a long load shows nothing meanwhile).
 */
export const StageEmptyDelayMs = 700;
/** How long the mark takes to fade out when something comes (faster than it came). */
export const StageMarkFadeOutMs = 250;
/** Below this width (px) the stage is compact: its toolbar folds into one button (`collapsed: 'auto'`). */
export const StageCompactWidth = 480;
/** Below this height (px) the stage is short: an `auto` hover label follows the pointer. */
export const StageShortHeight = 360;

export interface StageModelOptions {
  /** The canvas's background, a CSS color (default white). */
  background?: string;
  /** Nothing drawn yet (default true: a stage starts empty, and quiet). */
  empty?: boolean;
  /** The mark on an empty, idle stage (default true); off, the stage never goes quiet. */
  mark?: boolean;
  hoverPlacement?: HoverPlacement;
  emptyDelayMs?: number;
  markFadeOutMs?: number;
  timers?: StageTimers;
}

/** Shallow equality of two line lists (an equal list set anew changes nothing). */
function sameLines(a: readonly HoverLine[], b: readonly HoverLine[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/**
 * The stage (stage §1): the app's canvas and the chrome over it. Its background and the ink derived from it, whether
 * anything is drawn (`empty`) or work is under way (`busy`), the hover label's lines, the band's height (`topInset`:
 * the app frames its scene below it) and the stage's size, and the empty mark's phase with the canvas's `quiet`.
 */
export class StageModel {
  readonly #background: State<string>;
  readonly #empty: State<boolean>;
  readonly #busy = new State(false);
  readonly #lines = new State<readonly HoverLine[]>([], sameLines);
  readonly #placement: State<HoverPlacement>;
  readonly #topInset = new State(0);
  readonly #band = new State(true);
  readonly #size = new State<StageSize>(
    { width: 0, height: 0 },
    (a, b) => a.width === b.width && a.height === b.height,
  );
  readonly #mark = new State<StageMarkPhase>('gone');
  /** Whether the canvas shows: false while the mark is up, and before anything was drawn. */
  readonly #shown: State<boolean>;

  /** The canvas's background (a CSS color). */
  readonly background: ReadonlyState<string>;
  /** The chrome's theme and color on the background (`inkForColor`). */
  readonly ink: ReadonlyState<CanvasInk>;
  readonly empty: ReadonlyState<boolean>;
  readonly busy: ReadonlyState<boolean>;
  /** The hover label's model: its lines. */
  readonly hover: HoverLabelModel;
  readonly hoverPlacement: ReadonlyState<HoverPlacement>;
  /**
   * Where the hover label is (stage §3.1): `hoverPlacement`, with `auto` decided by the size: beside the pointer while
   * the stage is narrower than `StageCompactWidth` or shorter than `StageShortHeight`, else the corner (the corner too
   * while the size is unknown: headless).
   */
  readonly hoverAt: ReadonlyState<HoverAt>;
  /** The band's height (px) over the canvas, its fade not counted; 0 without one or headless. */
  readonly topInset: ReadonlyState<number>;
  /** Whether the band (the `top` slot) shows. */
  readonly band: ReadonlyState<boolean>;
  readonly size: ReadonlyState<StageSize>;
  /** Narrower than `StageCompactWidth` (false while the width is unknown). */
  readonly compact: ReadonlyState<boolean>;
  readonly mark: ReadonlyState<StageMarkPhase>;
  /** The canvas is faded out: empty, under the mark (or never drawn yet). */
  readonly quiet: ReadonlyState<boolean>;

  readonly #markOn: boolean;
  readonly #delay: number;
  readonly #fadeOut: number;
  readonly #timers: StageTimers;
  #timer: unknown;
  #pending: StageMarkPhase = 'gone';
  #disposed = false;

  constructor(options: StageModelOptions = {}) {
    const empty = options.empty ?? true;
    this.#background = new State(options.background ?? '#ffffff');
    this.#empty = new State(empty);
    this.#markOn = options.mark ?? true;
    this.#shown = new State(!empty || !this.#markOn);
    this.#placement = new State(options.hoverPlacement ?? 'auto');
    this.#delay = options.emptyDelayMs ?? StageEmptyDelayMs;
    this.#fadeOut = options.markFadeOutMs ?? StageMarkFadeOutMs;
    this.#timers = options.timers ?? GlobalTimers;

    this.background = this.#background;
    this.ink = derive(this.#background, inkForColor, (a, b) => a.theme === b.theme && a.color === b.color);
    this.empty = this.#empty;
    this.busy = this.#busy;
    this.hover = { lines: this.#lines };
    this.hoverPlacement = this.#placement;
    this.hoverAt = derive([this.#placement, this.#size], (p, s): HoverAt => {
      if (p !== 'auto') return p;
      const small = (s.width > 0 && s.width < StageCompactWidth) || (s.height > 0 && s.height < StageShortHeight);
      return small ? 'cursor' : 'corner';
    });
    this.topInset = this.#topInset;
    this.band = this.#band;
    this.size = this.#size;
    this.compact = derive(this.#size, (s) => s.width > 0 && s.width < StageCompactWidth);
    this.mark = this.#mark;
    this.quiet = derive(this.#shown, (s) => !s);
    this.#decide();
  }

  /** The canvas's background, a CSS color: the chrome follows it (`ink`). */
  setBackground(color: string): void {
    this.#background.set(color);
  }

  /** Whether anything is drawn: empty and idle a while, the mark shows; something drawn, the canvas shows. */
  setEmpty(empty: boolean): void {
    if (this.#empty.value === empty) return;
    batch(() => {
      this.#empty.set(empty);
      if (!empty) this.#shown.set(true);
    });
    this.#decide();
  }

  /** Work under way (a load, a task): an empty stage's wait for the mark starts over when it ends. */
  setBusy(busy: boolean): void {
    if (this.#busy.value === busy) return;
    this.#busy.set(busy);
    this.#decide();
  }

  /** The hover label's lines (empty: nothing hovered). */
  setHover(lines: readonly HoverLine[]): void {
    this.#lines.set(lines);
  }

  clearHover(): void {
    if (this.#lines.value.length) this.#lines.set([]);
  }

  setHoverPlacement(placement: HoverPlacement): void {
    this.#placement.set(placement);
  }

  /** The band's height (px, its fade not counted), as measured by the view. */
  setTopInset(px: number): void {
    this.#topInset.set(Math.max(0, Math.round(px)));
  }

  /** Shows or hides the band (the `top` slot): hidden, it is paused and inert, and `topInset` is 0. */
  showBand(on: boolean): void {
    this.#band.set(on);
    if (!on) this.#topInset.set(0);
  }

  /** The stage's size (px), as measured by the view. */
  setSize(width: number, height: number): void {
    this.#size.set({ width: Math.round(width), height: Math.round(height) });
  }

  /** Stops the mark's wait. */
  dispose(): void {
    this.#disposed = true;
    this.#cancel();
  }

  #cancel(): void {
    if (this.#timer !== undefined) this.#timers.clearTimeout(this.#timer);
    this.#timer = undefined;
  }

  #after(ms: number, next: StageMarkPhase): void {
    this.#cancel();
    this.#pending = next;
    this.#timer = this.#timers.setTimeout(() => {
      this.#timer = undefined;
      this.#step(this.#pending);
    }, ms);
  }

  /** Fading out only what showed; the mark up quiets the canvas, going out shows it. */
  #step(next: StageMarkPhase): void {
    const prev = this.#mark.value;
    const phase = next === 'out' && prev !== 'in' ? 'gone' : next;
    if (phase === prev) return;
    batch(() => {
      this.#mark.set(phase);
      if (phase === 'in') this.#shown.set(false);
      else if (phase === 'out') this.#shown.set(true);
    });
  }

  /**
   * One wait at a time (the next phase): drawn, the mark fades out; empty and idle, it comes after the delay; busy,
   * nothing changes (a mark up stays, one not yet up waits) and the wait starts over after.
   */
  #decide(): void {
    this.#cancel();
    if (this.#disposed || !this.#markOn) return;
    if (!this.#empty.value) {
      this.#step('out');
      if (this.#mark.value === 'out') this.#after(this.#fadeOut, 'gone');
    } else if (!this.#busy.value && this.#mark.value !== 'in') this.#after(this.#delay, 'in');
  }
}
