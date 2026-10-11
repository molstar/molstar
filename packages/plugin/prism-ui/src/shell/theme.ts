/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

/** The UI themes (theme.css): `auto` follows the system's light or dark. */
export type UITheme = 'light' | 'dark' | 'blue' | 'auto';

export const UIThemes: readonly { id: UITheme; title: string; icon: IconName }[] = [
  { id: 'light', title: 'Light', icon: 'sun' },
  { id: 'dark', title: 'Dark', icon: 'moon' },
  { id: 'blue', title: 'Blue', icon: 'palette' },
  { id: 'auto', title: 'Auto (system)', icon: 'sun-moon' },
];

/** The UI theme, switched from the rail (by Settings): set as `data-msp-theme` on the element the UI is themed by. */
export class ThemeModel {
  readonly #theme;

  constructor(initial: UITheme = 'light') {
    this.#theme = observableState<UITheme>(initial);
  }

  get theme(): ReadonlyState<UITheme> {
    return this.#theme.asReadonly;
  }

  get value(): UITheme {
    return this.#theme.value;
  }

  set(theme: UITheme) {
    this.#theme.set(theme);
  }

  /** The next theme in `UIThemes` order (the rail button cycles through them). */
  next() {
    const i = UIThemes.findIndex((t) => t.id === this.#theme.value);
    this.#theme.set(UIThemes[(i + 1) % UIThemes.length].id);
  }

  /** Keeps `data-msp-theme` on `el` in step; returns the unsubscribe. */
  applyTo(el: HTMLElement): () => void {
    const sub = this.#theme.asReadonly.subscribe((t) => {
      el.dataset.mspTheme = t;
    });
    return () => sub.unsubscribe();
  }
}
