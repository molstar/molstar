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
import { type Issue, issueOf, type Node } from '@molstar/core/repr/tree';
import {
  type AssetReader,
  DataStore,
  type FormatRuntime,
  isAbortError,
  type OpenedSource,
  type SourceEntry,
  SourceIssues,
  SourceOpenError,
  type StoreClock,
} from '@molstar/seqx-core/data';
import { OverlaySources } from '@molstar/seqx-core/render';
import {
  EmptyOverlay,
  nodeAddress,
  type Overlay,
  OverlayBuilder,
  type OverlayNode,
  type StaticModel,
} from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import { DomHost } from '@molstar/ui/dom';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { ViewerIssues } from '@molstar/seqx-viewer/issues';

/** The store's time through `DomHost`, read at each call. */
export const DomHostClock: StoreClock = {
  now: () => DomHost.now(),
  setTimeout: (fn, ms) => DomHost.setTimeout(fn, ms),
  clearTimeout: (handle) => DomHost.clearTimeout(handle),
};

/** A source's nodes the store holds now, by key. */
function nodesOf(store: DataStore, source: string): Map<string, Node> {
  const out = new Map<string, Node>();
  for (const n of store.nodes.values()) if (n.source === source) out.set(n.key, n.node);
  return out;
}

/** The addresses of the nodes gone between two readings, and of the children (units) a kept node lost. */
function goneSince(before: ReadonlyMap<string, Node>, after: ReadonlyMap<string, Node>, source: string): string[] {
  const out: string[] = [];
  for (const [key, node] of before) {
    const address = nodeAddress(source, key);
    const now = after.get(key);
    if (!now) out.push(address);
    else if (now !== node) {
      const ids = new Set((now.children ?? []).map((c) => c.id));
      for (const c of node.children ?? []) if (c.id !== undefined && !ids.has(c.id)) out.push(`${address}/${c.id}`);
    }
  }
  return out;
}

/** Why a source could not be opened, at its address. */
function unopened(e: unknown, address: string): Issue {
  const issue =
    e instanceof SourceOpenError && e.issues[0]
      ? e.issues[0]
      : issueOf(SourceIssues, 'seq.source.unreadable', `The source could not be opened: ${String(e)}`, { address });
  return { ...issue, address };
}

interface Opening {
  readonly entry: SourceEntry;
  readonly abort: AbortController;
  /** Stops listening to the opened source's `changes`. */
  stop?: () => void;
  /** A reopen is asked for the next frame (changes inside one frame reopen once). */
  asked?: boolean;
  /** The latest reopen: an earlier one that settles after it is dropped. */
  reading?: number;
}

export class StoreModel extends ReactiveModel {
  readonly store: DataStore;
  private readonly _open = new Map<string, Opening>();
  private _overlays: OverlayBuilder | undefined;
  private _overlaySources: OverlaySources | undefined;
  private _synced: { overlay: Overlay; view: ResolvedView } | undefined;
  /** The addresses of the overlay's sources reported as unopened: cleared when they open or their node goes. */
  private _overlayFailed = new Set<string>();
  private readonly _issues = new State<readonly Issue[]>([]);
  /** Sources that could not be opened, one issue each, and what an opened source said (rows a file of records left
   * out); the store's own reports are `store.issues`. */
  readonly issues: ReadonlyState<readonly Issue[]> = this._issues.asReadonly;

  /**
   * `removed`: told the addresses of the nodes a host's update removed (and of the children a kept node lost), it takes
   * them out of the session's state and says whether something was (then one issue at the source).
   */
  constructor(
    config: ViewerConfig = DefaultViewerConfig,
    private readonly options: { clock?: StoreClock; removed?: (gone: readonly string[]) => boolean } = {},
  ) {
    super();
    this.store = new DataStore({
      settings: {
        maxBytes: config.store.cacheBytes,
        maxInFlight: config.store.inFlight,
        tileBins: config.store.tileBins,
        maxColumnLength: config.bounds.maxColumnLength,
        maxBlockValues: config.bounds.maxBlockValues,
      },
      clock: this.options.clock ?? DomHostClock,
    });
    this.registerPermanentDispose(() => {
      // the overlay's openings first: one that settles later is disposed, never attached
      this._overlaySources?.dispose();
      for (const held of this._sessionHolds.values()) for (const release of held.values()) release();
      this._sessionHolds.clear();
      for (const o of this._open.values()) {
        o.stop?.();
        o.abort.abort();
      }
      this._open.clear();
      this.store.retire();
    });
  }

  /**
   * Opens the static model's sources (each `source` bound to a runtime, and each inline access plan) and attaches each
   * to the store under its address; a source whose entry is the object it was opened from is kept, one gone is
   * detached (its issue cleared). A source that failed to open, or whose coverage the store refused, is tried again
   * by the next sync. Resolves when the openings started now have settled.
   */
  async sync(model: StaticModel, assets: AssetReader): Promise<void> {
    const wanted = this._detachStale(model);
    const started: Promise<void>[] = [];
    for (const [address, { entry, runtime }] of wanted) {
      if (this._open.has(address)) continue;
      const opening: Opening = { entry, abort: new AbortController() };
      this._open.set(address, opening);
      started.push(this._attach(address, opening, runtime, assets));
    }
    await Promise.all(started);
  }

  /**
   * Detaches at once what `model` no longer opens from the entry it was opened from (an edited source's old file, a
   * source gone), so no frame before the next sync reads it: what `sync` does first.
   */
  detachStale(model: StaticModel): void {
    this._detachStale(model);
  }

  private _detachStale(model: StaticModel) {
    const wanted = new Map<string, { entry: SourceEntry; runtime: FormatRuntime }>();
    for (const b of model.sources.values())
      if (b.runtime) wanted.set(b.address, { entry: b.entry, runtime: b.runtime });
    const plans = [...model.systems.values()].map((s) => s.access);
    for (const a of model.annotations.values()) {
      if (a.access) plans.push(a.access);
      // links' ends, a ranges layer's over them (core's _spec "Render", Layer types: ranges)
      if (a.ends) plans.push(...a.ends.values());
    }
    for (const p of plans)
      if (p.kind === 'inline') wanted.set(p.source.address, { entry: p.source, runtime: p.runtime });

    for (const [address, o] of [...this._open]) {
      if (wanted.get(address)?.entry === o.entry) continue;
      o.stop?.();
      o.abort.abort();
      this._open.delete(address);
      this.store.detach(address);
      this._report(address, undefined);
    }
    // a source that failed and is no longer wanted takes its issue with it
    const stale = this._issues.value.filter((i) => i.address !== undefined && !wanted.has(i.address));
    for (const i of stale) this._report(i.address!, undefined);
    return wanted;
  }

  /**
   * The coverage of a source of `model` when what the store holds under its address was opened from the entry the
   * model binds there now; undefined while it is not (an edited source's old file, one still opening).
   */
  openedCoverage(model: StaticModel, address: string) {
    const o = this._open.get(address);
    return o && o.entry === model.sources.get(address)?.entry ? this.store.coverage(address) : undefined;
  }

  /**
   * Tier 2 for a frame: the overlay of the store's dynamic nodes now (one a node version), its sources (a dynamic
   * annotation's rows, the nested sources the view's templates reach) opened when it or the view changes, and let go
   * with their nodes. What cannot be opened is an issue here, as a source of the model's.
   */
  overlay(model: StaticModel, view: ResolvedView, ctx: { runtimes: readonly FormatRuntime[]; assets: AssetReader }) {
    // the model's configuration is read through the model; other runtimes make another builder
    if (this._overlays?.options.runtimes !== ctx.runtimes)
      this._overlays = new OverlayBuilder({ runtimes: ctx.runtimes });
    const overlay = this._overlays.of(model, this.store);
    if (this._synced?.overlay !== overlay || this._synced.view !== view) {
      this._synced = { overlay, view };
      this._overlaySources ??= new OverlaySources(this.store, ctx.assets);
      void this._overlaySources.sync(overlay, view).then((issues) => {
        if (this.isRetired) return;
        const now = new Set<string>();
        for (const i of issues) {
          now.add(i.address ?? '');
          this._report(i.address ?? '', i);
        }
        for (const a of this._overlayFailed) if (!now.has(a)) this._report(a, undefined);
        this._overlayFailed = now;
      });
    }
    return overlay;
  }

  /**
   * The overlay of the nodes the store holds now over `model`, as the last frame's builder makes it (one per node
   * version), without syncing sources; no node before the first frame. What addresses resolve against.
   */
  currentOverlay(model: StaticModel): Overlay {
    return this._overlays ? this._overlays.of(model, this.store) : EmptyOverlay;
  }

  /** Holds the nodes a frame's layout shows (`ViewLayout.shown`), and lets go of those it no longer shows. */
  holdShown(nodes: readonly OverlayNode[]): void {
    this._overlaySources?.hold(nodes);
  }

  /** What each session state holds (`holdSession`), by its name: a release by `<source>#<key>`. */
  private readonly _sessionHolds = new Map<string, Map<string, () => void>>();

  /**
   * Holds the dynamic nodes a session state names (the selection's elements' owners and its systems, a highlight's),
   * with what they name (`Overlay.holds`), under `name`, and lets go of those it held before and names no more: a
   * selected node stays while it is selected, shown or not. None: lets go of all of `name`'s.
   */
  holdSession(name: string, model: StaticModel, addresses: Iterable<string>): void {
    if (this.isRetired) return;
    const nodes = this.currentOverlay(model).holds(addresses);
    let held = this._sessionHolds.get(name);
    if (!held) {
      if (!nodes.length) return;
      this._sessionHolds.set(name, (held = new Map()));
    }
    const now = new Set<string>();
    for (const n of nodes) {
      const k = `${n.source}#${n.key}`;
      now.add(k);
      if (!held.has(k)) held.set(k, this.store.hold(n.source, n.key));
    }
    for (const [k, release] of [...held]) {
      if (now.has(k)) continue;
      release();
      held.delete(k);
    }
  }

  /** The nodes a session state holds now, as `<source>#<key>`. */
  sessionHolds(name: string): readonly string[] {
    return [...(this._sessionHolds.get(name)?.keys() ?? [])];
  }

  private async _attach(address: string, o: Opening, runtime: FormatRuntime, assets: AssetReader): Promise<void> {
    try {
      const { maxRows, maxBytes } = this.store.settings;
      const opened = await runtime.open(o.entry, { assets, signal: o.abort.signal, limits: { maxRows, maxBytes } });
      if (o.abort.signal.aborted || this._open.get(address) !== o) {
        opened.dispose();
        return;
      }
      // what opening said without refusing it (rows a file of records left out), at its address
      this._report(
        address,
        opened.issues?.map((i) => ({ ...i, address })),
      );
      // a coverage the store refuses is its report (`store.issues`); forgotten here, so a later sync tries again
      if (!this.store.attach(address, opened)) {
        opened.dispose();
        this._open.delete(address);
      } else this._watch(address, o, opened, runtime, assets);
    } catch (e) {
      if (isAbortError(e) || o.abort.signal.aborted) return;
      // forgotten, so a later sync tries again
      if (this._open.get(address) === o) this._open.delete(address);
      this._report(address, unopened(e, address));
    }
  }

  /**
   * Listens to an opened source's `changes` (a host's edit): a change asks a reopen for the next frame (`DomHost`'s),
   * so changes inside one frame reopen once (_spec/README.md "StoreModel", Host updates).
   */
  private _watch(address: string, o: Opening, opened: OpenedSource, runtime: FormatRuntime, assets: AssetReader) {
    o.stop?.();
    o.stop = opened.changes?.(() => {
      if (o.asked) return;
      o.asked = true;
      DomHost.requestAnimationFrame(() => {
        o.asked = false;
        if (this._open.get(address) === o && !o.abort.signal.aborted) void this._reopen(address, o, runtime, assets);
      });
    });
  }

  /**
   * Opens a source again from its entry and reattaches it in place (`DataStore.reattach`: nodes kept by key, tiles read
   * again); a failed reopen is an issue at its address and the reading before it stays.
   */
  private async _reopen(address: string, o: Opening, runtime: FormatRuntime, assets: AssetReader): Promise<void> {
    const reading = (o.reading = (o.reading ?? 0) + 1);
    const current = () => this._open.get(address) === o && !o.abort.signal.aborted && o.reading === reading;
    try {
      const { maxRows, maxBytes } = this.store.settings;
      const opened = await runtime.open(o.entry, { assets, signal: o.abort.signal, limits: { maxRows, maxBytes } });
      if (!current()) return opened.dispose();
      this._report(
        address,
        opened.issues?.map((i) => ({ ...i, address })),
      );
      const before = nodesOf(this.store, address);
      if (!this.store.reattach(address, opened)) return opened.dispose();
      this._watch(address, o, opened, runtime, assets);
      const gone = goneSince(before, nodesOf(this.store, address), address);
      if (gone.length && this.options.removed?.(gone))
        this._report(address, [
          ...(opened.issues?.map((i) => ({ ...i, address })) ?? []),
          issueOf(
            ViewerIssues,
            'seq.viewer.selection.dropped',
            `The source’s update removed ${gone.length === 1 ? 'a node' : `${gone.length} nodes`} the selection named: dropped from it`,
            { address },
          ),
        ]);
    } catch (e) {
      if (isAbortError(e) || !current()) return;
      this._report(address, unopened(e, address));
    }
  }

  /** The issues at `address` replaced: one, a list (an opened source's), or none. */
  private _report(address: string, issues: Issue | readonly Issue[] | undefined): void {
    const rest = this._issues.value.filter((i) => i.address !== address);
    const now = issues === undefined ? [] : Array.isArray(issues) ? issues : [issues as Issue];
    this._issues.set(now.length ? [...rest, ...now] : rest);
  }
}
