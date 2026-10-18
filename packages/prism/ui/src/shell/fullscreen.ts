/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

/**
 * True fullscreen for the plugin (`requestFullscreen`), next to Expand (which fills the page). Embedded in an iframe,
 * expanding only fills the frame, so there fullscreen is what Expand should do; it needs the embedding page to allow
 * it (`<iframe allow="fullscreen">`), and without that it is not offered (molstar#1663). Leaving fullscreen by Escape
 * (or the browser's own controls) is followed.
 *
 * Several viewers on a page (viewer README decision 29): the state is that of this model's own element only. Another
 * element going fullscreen (or leaving it) does not flip this one, and its button enters fullscreen even while another
 * viewer is in it (the browser moves it), instead of leaving the other's. `fullscreenchange` is heard on the document
 * (it bubbles from any element), but only recomputed against the own element.
 */
export class FullscreenModel {
  readonly #active = observableState(false);
  readonly #onChange = () => this.#active.set(this.#isOwn());

  /** Whether the page may go fullscreen at all (false in an iframe that does not allow it). */
  readonly available: boolean;
  /** Running inside an iframe: Expand fills only the frame. */
  readonly embedded: boolean;

  constructor(private readonly target: () => Element | undefined) {
    this.available = typeof document !== 'undefined' && !!document.fullscreenEnabled;
    this.embedded = typeof window !== 'undefined' && window.self !== window.top;
    if (typeof document !== 'undefined') document.addEventListener('fullscreenchange', this.#onChange);
  }

  /** Whether the own element is the fullscreen one: asked of its own root, so an element inside a shadow root counts. */
  #isOwn(): boolean {
    const element = this.target();
    if (!element || typeof document === 'undefined') return false;
    const root = element.getRootNode();
    const current = (root as Partial<DocumentOrShadowRoot>).fullscreenElement;
    return (current ?? document.fullscreenElement) === element;
  }

  /**
   * What the canvas's Expand button does here: `fullscreen` embedded in an iframe that allows it (expanding would only
   * fill the frame); otherwise `expand` (fill the page or the frame), with fullscreen offered beside it.
   */
  get mode(): 'expand' | 'fullscreen' {
    return this.embedded && this.available ? 'fullscreen' : 'expand';
  }

  get active(): ReadonlyState<boolean> {
    return this.#active.asReadonly;
  }

  /** Why it cannot be used here, or undefined. */
  get unavailableReason(): string | undefined {
    return this.available ? undefined : 'Not allowed by the page Mol* is embedded in';
  }

  async toggle() {
    if (!this.available) return;
    try {
      if (this.#isOwn()) await document.exitFullscreen();
      else await this.target()?.requestFullscreen();
    } catch (e) {
      console.error(e);
    }
  }

  dispose() {
    if (typeof document !== 'undefined') document.removeEventListener('fullscreenchange', this.#onChange);
  }
}
