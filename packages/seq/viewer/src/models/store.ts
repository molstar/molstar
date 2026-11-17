/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The data store's model (_spec/README.md "StoreModel"): the viewer's one data store, made with the configuration's
// bounds and timed by `DomHost` (so a test's patched host drives it), and the sources of the static model opened and
// attached under their addresses: a source kept while its entry is the same object (the resolver keeps entries by
// hash), opened again when it changes, detached when it goes. Retired with the viewer.

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import {
  type AssetReader,
  DataStore,
  type FormatRuntime,
  isAbortError,
  type SourceEntry,
  SourceIssues,
  SourceOpenError,
  type StoreClock,
} from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { DomHost } from '@molstar/ui/dom';
import { DefaultViewerConfig, type ViewerConfig } from '../config.js';

/** The store's time through `DomHost`, read at each call. */
export const DomHostClock: StoreClock = {
  now: () => DomHost.now(),
  setTimeout: (fn, ms) => DomHost.setTimeout(fn, ms),
  clearTimeout: (handle) => DomHost.clearTimeout(handle),
};

interface Opening {
  readonly entry: SourceEntry;
  readonly abort: AbortController;
}

export class StoreModel extends ReactiveModel {
  readonly store: DataStore;
  private readonly _open = new Map<string, Opening>();
  private readonly _issues = new State<readonly Issue[]>([]);
  /** Sources that could not be opened, one issue each (the store's own reports are `store.issues`). */
  readonly issues: ReadonlyState<readonly Issue[]> = this._issues.asReadonly;

  constructor(config: ViewerConfig = DefaultViewerConfig, options: { clock?: StoreClock } = {}) {
    super();
    this.store = new DataStore({
      settings: {
        maxBytes: config.store.cacheBytes,
        maxInFlight: config.store.inFlight,
        tileBins: config.store.tileBins,
      },
      clock: options.clock ?? DomHostClock,
    });
    this.registerPermanentDispose(() => {
      for (const o of this._open.values()) o.abort.abort();
      this._open.clear();
      this.store.retire();
    });
  }

  /**
   * Opens the static model's sources (each `source` bound to a runtime, and each inline access plan) and attaches each
   * to the store under its address; a source whose entry is the object it was opened from is kept, one gone is
   * detached. Resolves when the openings started now have settled.
   */
  async sync(model: StaticModel, assets: AssetReader): Promise<void> {
    const wanted = new Map<string, { entry: SourceEntry; runtime: FormatRuntime }>();
    for (const b of model.sources.values())
      if (b.runtime) wanted.set(b.address, { entry: b.entry, runtime: b.runtime });
    const plans = [...model.systems.values()].map((s) => s.access);
    for (const a of model.annotations.values()) plans.push(a.access);
    for (const p of plans)
      if (p.kind === 'inline') wanted.set(p.source.address, { entry: p.source, runtime: p.runtime });

    for (const [address, o] of [...this._open]) {
      if (wanted.get(address)?.entry === o.entry) continue;
      o.abort.abort();
      this._open.delete(address);
      this.store.detach(address);
    }
    const started: Promise<void>[] = [];
    for (const [address, { entry, runtime }] of wanted) {
      if (this._open.has(address)) continue;
      const opening: Opening = { entry, abort: new AbortController() };
      this._open.set(address, opening);
      started.push(this._attach(address, opening, runtime, assets));
    }
    await Promise.all(started);
  }

  private async _attach(address: string, o: Opening, runtime: FormatRuntime, assets: AssetReader): Promise<void> {
    try {
      const opened = await runtime.open(o.entry, { assets, signal: o.abort.signal });
      if (o.abort.signal.aborted || this._open.get(address) !== o) {
        opened.dispose();
        return;
      }
      this._report(address, undefined);
      this.store.attach(address, opened);
    } catch (e) {
      if (isAbortError(e) || o.abort.signal.aborted) return;
      const issue =
        e instanceof SourceOpenError && e.issues[0]
          ? e.issues[0]
          : issueOf(SourceIssues, 'seq.source.unreadable', `The source could not be opened: ${String(e)}`, {
              address,
            });
      this._report(address, { ...issue, address });
    }
  }

  private _report(address: string, issue: Issue | undefined): void {
    const rest = this._issues.value.filter((i) => i.address !== address);
    this._issues.set(issue ? [...rest, issue] : rest);
  }
}
