/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ElementRoot, ViewElement } from '../dom/element.js';
import '../styles/base.js';
import { adoptStyles, releaseStyles } from '../styles/registry.js';

/**
 * A `ViewElement` that brings the UI's styles (styles §8): on connect it adopts the registered styles into its
 * document or shadow root, and releases them once it has left that root for good; a top-level one (not inside another
 * UI root) is a root (`.msui-root`), which scopes the reset.
 */
export abstract class UiElement<M = unknown> extends ViewElement<M> {
  private _adopted: ElementRoot | undefined = undefined;

  protected override onConnected(root: ElementRoot): void {
    if (root !== this._adopted) {
      if (this._adopted) releaseStyles(this._adopted);
      adoptStyles(root);
      this._adopted = root;
    }
    if (!this.parentElement?.closest('.msui-root')) this.classList.add('msui-root');
  }

  protected override onReleased(root: ElementRoot): void {
    if (root !== this._adopted) return;
    releaseStyles(root);
    this._adopted = undefined;
  }
}
