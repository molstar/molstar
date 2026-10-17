/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The splash fields' rest and wake (David, 2026-10-17: prism AGENTS.md, idle is idle; the splash no longer the
// exception): the light's particles and the wave run at full rate while something happens, and fade out when nothing
// has for a while, then stop for good: no frame is asked for at rest (the glass wordmark stays as the last frame drew
// it). A pointer move or press, a resize wakes them: they fade back in and the steps resume where they were.
//
// Rest is all of: no input for `ActiveMs`; the glass not animating (a shatter, a gathering, a lean settling); at least
// `RestAfterMs` of running since the start (the intro needs time to play; time hidden does not count, a frame only
// comes while shown). The fade is by alpha in the field's own draw (`level`, 0..1), in `FadeMs`.

/** Full rate this long after the last input. */
export const ActiveMs = 1500;
/** The splash runs at least this long before it may rest (the intro plays out). */
export const RestAfterMs = 12_000;
/** How long the particles and the wave take to fade out or in. */
export const FadeMs = 400;
/** The longest a frame counts for the fade and the clock (a stalled tab does not skip the fade or the intro). */
const MaxFrameMs = 100;

export interface AmbientPaceOptions {
  /** Defaults to `ActiveMs`, `RestAfterMs`, `FadeMs` (tests shorten them). */
  activeMs?: number;
  restAfterMs?: number;
  fadeMs?: number;
}

// One clock for input and frames (`performance.now()`): a frame's own timestamp may be on another (jsdom's is).
const now = () => performance.now();

export class AmbientPace {
  readonly #activeMs: number;
  readonly #restAfterMs: number;
  readonly #fadeMs: number;
  #input = Number.NEGATIVE_INFINITY;
  #tick = Number.NaN;
  /** Milliseconds run (frames that came, each clamped). */
  #age = 0;
  #level = 1;
  #asleep = false;
  /** Steps taken (a probe: the guards count them). */
  steps = 0;

  constructor(options: AmbientPaceOptions = {}) {
    this.#activeMs = options.activeMs ?? ActiveMs;
    this.#restAfterMs = options.restAfterMs ?? RestAfterMs;
    this.#fadeMs = options.fadeMs ?? FadeMs;
  }

  /** How much of the particles and the wave shows (0..1): the alpha of their draw. */
  get level(): number {
    return this.#level;
  }

  /** Faded out and stopped: nothing to draw or ask a frame for until `input`. */
  get asleep(): boolean {
    return this.#asleep;
  }

  /** Input came (the pointer moved or pressed, the field was resized): awake, at full rate for a while. */
  input(): void {
    this.#input = now();
    this.#asleep = false;
  }

  /**
   * One frame: moves the fade toward where it goes (out when at rest, else in) and returns its level; at 0 the pace
   * is asleep. `animating`: the glass is shattering, gathering or settling (keeps it awake).
   */
  advance(animating: boolean): number {
    const t = now();
    const dt = Number.isNaN(this.#tick) ? 0 : Math.min(MaxFrameMs, Math.max(0, t - this.#tick));
    this.#tick = t;
    this.#age += dt;
    const rest = !animating && t - this.#input >= this.#activeMs && this.#age >= this.#restAfterMs;
    if (rest) {
      this.#level = Math.max(0, this.#level - dt / this.#fadeMs);
      if (this.#level === 0) this.#asleep = true;
    } else {
      this.#level = Math.min(1, this.#level + dt / this.#fadeMs);
    }
    return this.#level;
  }

  stepped(): void {
    this.steps++;
  }

  /**
   * Started again after a pause (hidden, then shown): the first frame measures from here, not from before the pause;
   * one that was awake gets its active window again (it comes back, then rests), one asleep stays so.
   */
  reset(): void {
    this.#tick = Number.NaN;
    if (!this.#asleep) this.#input = now();
  }
}
