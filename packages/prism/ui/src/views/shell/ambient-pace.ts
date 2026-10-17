/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The splash fields' pace (David, 2026-10-17: prism AGENTS.md, idle is idle; the splash the exception): at full rate
// while the pointer moves over them and a little after, at `IdleFps` or so otherwise. Through the page's one frame
// tick, by skipping frames (a frame callback every frame, a step and a draw only when due), not by timers: the wave and
// the light move as they did (a step's `dt` is the time since the last step), only less often.

/** The rate at rest: a step when this much has passed since the last (about 20 a second at 60 Hz, 24 at 120). */
export const IdleStepMs = 1000 / 24 - 4;
/** Full rate this long after the last input. */
export const ActiveMs = 1500;

// One clock for input and frames (`performance.now()`): a frame's own timestamp may be on another (jsdom's is).
const now = () => performance.now();

export class AmbientPace {
  #input = Number.NEGATIVE_INFINITY;
  #stepped = Number.NEGATIVE_INFINITY;
  /** Steps taken (a probe: the idle rate's guard counts them). */
  steps = 0;

  /** Input came (the pointer moved, pressed; the field was resized): full rate for a while. */
  input(): void {
    this.#input = now();
  }

  /** Whether to step and draw in this frame. */
  due(): boolean {
    const t = now();
    return t - this.#input < ActiveMs || t - this.#stepped >= IdleStepMs;
  }

  stepped(): void {
    this.#stepped = now();
    this.steps++;
  }

  /** Started again (shown, a resume): the next frame steps. */
  reset(): void {
    this.#stepped = Number.NEGATIVE_INFINITY;
  }
}
