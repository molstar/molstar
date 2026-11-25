/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Session Locations carried onto the view's sections (_spec/README.md "Emphasis"): the selection and the highlight
// channels projected as the core's `projectOntoView` does it, computed once per Location (by identity) and kept while
// the view and the model hold, so a frame that redraws the emphasis projects nothing again.

import type { Location } from '@molstar/seqx-core/query';
import { type ProjectOptions, projectOntoView, type RegionProjector } from '@molstar/seqx-core/view';
import type { SpecModel } from './spec.js';
import type { ViewModel } from './view.js';

/** The default options, one object: a call per frame makes none. */
const NotOwn: Pick<ProjectOptions, 'own'> = Object.freeze({ own: false });

/** The most Locations kept: past it, the cache starts again (a session holds a few). */
const MaxKept = 64;

export class Projector {
  private _view: unknown;
  private _model: unknown;
  private _held = -1;
  private readonly _kept = [new Map<Location, Location>(), new Map<Location, Location>()];

  constructor(
    private readonly spec: SpecModel,
    private readonly view: ViewModel,
    /** What carries a Location through a map read by region (a liftover): the store's runs held now. */
    private readonly regions?: RegionProjector,
  ) {}

  /**
   * `loc` on the view's sections' systems: with `own`, the parts it has there kept (a highlight's whole look), else only
   * what it adds where it has no part (a selection's projection). Undefined when no view is shown.
   */
  of(loc: Location, options: Pick<ProjectOptions, 'own'> = NotOwn): Location | undefined {
    const view = this.view.resolved.value;
    const model = this.spec.model.value;
    if (!view) return undefined;
    const held = this.regions?.reads(model) ? this.regions.version.value : -1;
    if (view !== this._view || model !== this._model || held !== this._held) {
      this._view = view;
      this._model = model;
      this._held = held;
      for (const m of this._kept) m.clear();
    }
    const kept = this._kept[options.own ? 1 : 0];
    let out = kept.get(loc);
    if (out === undefined) {
      if (kept.size >= MaxKept) kept.clear();
      out = projectOntoView(model, view, loc, {
        own: !!options.own,
        ...(this.regions ? { regions: this.regions } : {}),
      });
      kept.set(loc, out);
    }
    return out;
  }
}
