/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { railSections } from '@molstar/seqx-core/render';
import type { ResolvedNavigator, ResolvedView } from '@molstar/seqx-core/view';
import { windowIdle } from '@molstar/seqx-viewer/behaviors/window';
import { ViewerIssues } from './issues.js';
import { narrowedVersion, rowsRailsOf } from './narrowed.js';
import type { NavigatorDef, NavigatorModel } from './navigator.js';
import type { SeqNavigatorsPart } from './navigators.js';
import type { SeqViewer } from './viewer.js';

// The viewer's rails and navigators (_spec/README.md "Composition", Navigators): the ruler rails a frame draws (every
// linear section no shown navigator binds, every rows track no navigator binds along its rows), which of them can
// move, whether the reader shows them; the navigators' models by id, the view's declared ones and a host's.

export class ViewerRails {
  /**
   * The ruler rails shown: the reader's session (the embedded layout's corner button), never the spec;
   * `navigator.form: none` draws none whatever it says.
   */
  private readonly _shown = new State(true);
  readonly shown: ReadonlyState<boolean> = this._shown.asReadonly;
  /** The sections that draw a rail and which can move: kept while they hold (the frame's key reads them by identity). */
  rails: Record<string, true> | undefined;
  active: Record<string, boolean> = {};
  private readonly _railsFor: [unknown, boolean, number, number] = [undefined, false, 0, 0];
  /** Bumped when a host's navigator comes or goes: the rails it replaces are found again. */
  private _hostedVersion = 0;
  /** The rails' keys that are rows tracks' (a rows rail). */
  private _rowsRails: ReadonlySet<string> = new Set();
  /** The navigators' models by id: the view's declared ones as asked for, and a host's (`navigator(def)`). */
  private readonly _navigators = new Map<string, NavigatorModel>();

  constructor(
    private readonly viewer: SeqViewer,
    /** The navigators the spec lists (`navigators`): none, a declared navigator is not drawn and its rail stays. */
    readonly part: SeqNavigatorsPart | undefined,
    private readonly issues: State<readonly Issue[]>,
    private readonly register: (dispose: () => void) => void,
  ) {}

  /**
   * The rails this frame draws: every linear section no shown navigator binds, and every rows track no navigator binds
   * along its rows (a rows rail), while `navigator.form` is `rail` and the reader shows them; each active when the view
   * can move it (the `navigator` behavior listed, not idle: idle is the whole axis at the zoom limit, a follower under
   * `ignore`; a rows rail's every row shown, the layout's).
   */
  update(resolved: ResolvedView): void {
    const viewer = this.viewer;
    const on = viewer.config.navigator.form === 'rail' && this._shown.value;
    const f = this._railsFor;
    if (f[0] !== resolved || f[1] !== on || f[2] !== this._hostedVersion || f[3] !== narrowedVersion(viewer)) {
      f[0] = resolved;
      f[1] = on;
      f[2] = this._hostedVersion;
      f[3] = narrowedVersion(viewer);
      // a host's navigator (`navigator(def)`) replaces the rail on what it binds, as a declared one does
      const hosted: ResolvedNavigator[] = [];
      for (const m of this._navigators.values()) {
        const n = m.isHosted ? m.navigator : undefined;
        if (n) hosted.push(n);
      }
      // a declared navigator replaces a rail only when the spec lists the navigators (else it is not drawn: said once)
      const shown = this.part ? resolved : this._withoutNavigators(resolved);
      // and every rows layer's rows rail, by its track (`navigator.rowsRail`; a navigator a narrow frame hides, hidden)
      const rows = on ? rowsRailsOf(viewer, shown, hosted, viewer.config) : [];
      this._rowsRails = new Set(rows);
      this.rails = on
        ? Object.fromEntries([...railSections(shown, hosted), ...rows].map((a) => [a, true as const]))
        : undefined;
      this.active = {};
    }
    const rails = this.rails;
    if (!rails) return;
    const behavior = viewer.behaviors.has('navigator');
    let changed = false;
    // a rows rail can move while the behavior is there (idle while every row shows: the layout knows, `whole`)
    // (no closure: a frame makes none)
    for (const a in rails) if (this.active[a] !== this._can(behavior, a)) changed = true;
    if (!changed) return;
    const next: Record<string, boolean> = {};
    for (const a in rails) next[a] = this._can(behavior, a);
    this.active = next;
  }

  private _can(behavior: boolean, address: string): boolean {
    return behavior && (this._rowsRails.has(address) || !windowIdle(this.viewer.view, address, this.viewer.config));
  }

  /**
   * The view as drawn by a viewer that lists no navigators: its declared ones left out (their rails stay), and one
   * issue (`seq.viewer.navigators.none`) when it declares any shown.
   */
  private _withoutNavigators(resolved: ResolvedView): ResolvedView {
    if (!resolved.navigators.some((n) => !n.hidden)) return resolved;
    if (!this.issues.value.some((i) => i.code === 'seq.viewer.navigators.none'))
      this.issues.set([
        ...this.issues.value,
        issueOf(
          ViewerIssues,
          'seq.viewer.navigators.none',
          'The view declares navigators and the viewer lists none: they are not drawn, the rails stay',
          { address: resolved.address, hint: ViewerIssues.codes['seq.viewer.navigators.none'].hint },
        ),
      ]);
    return { ...resolved, navigators: [] };
  }

  /** Shows or hides the ruler rails (the reader's session). */
  setShown(shown: boolean): void {
    if (this._shown.value === shown) return;
    this._shown.set(shown);
    this.viewer.schedule();
  }

  /** A navigator's model by its id (`SeqViewer.navigatorModel`): made once and kept with the viewer. */
  model(id: string): NavigatorModel | undefined {
    const kept = this._navigators.get(id);
    if (kept) return kept;
    const view = this.viewer.view.resolved.value;
    const part = this.part;
    if (!part || !view?.navigators.some((n) => n.id === id || n.address === id)) return undefined;
    const m = part.model(this.viewer, id);
    this._navigators.set(id, m);
    this.register(() => m.retire());
    return m;
  }

  /** A host's navigator (`SeqViewer.navigator`): replaces one of its id, the rails it binds found again. */
  hosted(def: NavigatorDef): NavigatorModel {
    const part = this.part;
    if (!part) throw new Error(`The viewer lists no navigators (its spec's \`navigators\`): no navigator ${def.id}`);
    this._navigators.get(def.id)?.retire();
    const m = part.model(this.viewer, def.id, def);
    this._navigators.set(def.id, m);
    this._hostedVersion++;
    this.register(() => m.retire());
    this.viewer.schedule();
    return m;
  }
}
