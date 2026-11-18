/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import {
  type AssetReader,
  type DataStore,
  isAbortError,
  nodeKey,
  type OpenedSource,
  SourceOpenError,
} from '@molstar/seqx-core/data';
import type { Overlay, OverlayNode } from '@molstar/seqx-core/resolve';
import { type ResolvedView, templateNeeds } from '@molstar/seqx-core/view';
import { RenderIssues } from './issues.js';

// The overlay's sources (_spec/README.md "Render", The overlay's sources): what the overlay needs read goes through the
// store like any source. A dynamic annotation's inline rows are a `seq` source under the node's address; a nested
// source is opened only when a template of the view needs its depth (the overlay lists it). Each is kept while its
// node is the same entry, and detached when its node goes, so a node's data lives as long as the node. What could not be
// opened is an issue while its node is held, opened again at the next sync, and cleared when it opens or its node goes.
// On a shared store (a viewer's, an export reading through it) the overlay sources share what one of them attached:
// each claims what it needs, and the last claim let go detaches it, so an export ending never takes away (and cancels
// the reads of) a nested source the live view still shows.

interface Held {
  /** The node entry it was opened for: a node that comes back is opened again. */
  readonly identity: unknown;
  readonly abort: AbortController;
  /** Its opening has not settled. */
  opening: boolean;
  /** Opened, and the store refused its coverage (the store's report): kept, not opened again until its node changes. */
  refused?: boolean;
  /** Why it could not be opened. */
  failed?: readonly Issue[];
}

/** A store's overlay sources: their claims by address, and the coverage one of them attached under each address. */
interface Shared {
  readonly claims: Map<string, number>;
  readonly attached: Map<string, unknown>;
}

const SharedOf = new WeakMap<DataStore, Shared>();

function sharedOf(store: DataStore): Shared {
  let s = SharedOf.get(store);
  if (!s) SharedOf.set(store, (s = { claims: new Map(), attached: new Map() }));
  return s;
}

export class OverlaySources {
  private readonly held = new Map<string, Held>();
  /** The store's holds on the nodes a layout shows, by `<source>#<key>`. */
  private readonly shown = new Map<string, () => void>();
  private readonly shared: Shared;

  constructor(
    readonly store: DataStore,
    private readonly assets: AssetReader,
  ) {
    this.shared = sharedOf(store);
  }

  /** Why the sources it holds could not be opened, now. */
  get issues(): Issue[] {
    const out: Issue[] = [];
    for (const h of this.held.values()) if (h.failed) out.push(...h.failed);
    return out;
  }

  /**
   * Opens what the overlay needs now (a failed one again) and lets go of what it no longer lists. Resolves, with the
   * issues of what cannot be opened now (`issues`), when the openings this call started have settled.
   */
  async sync(overlay: Overlay, view: ResolvedView): Promise<Issue[]> {
    const wanted = new Map<string, { identity: unknown; open: (signal: AbortSignal) => Promise<OpenedSource> }>();
    for (const n of overlay.nodes) {
      const plan = overlay.annotation(n.address)?.access;
      if (plan?.kind === 'inline')
        wanted.set(plan.source.address, {
          identity: n.entry,
          open: (signal) => plan.runtime.open(plan.source, { assets: this.assets, signal }),
        });
    }
    for (const nested of overlay.nested(templateNeeds(view))) {
      const runtime = nested.runtime;
      if (runtime)
        wanted.set(nested.address, {
          identity: overlay.node(nested.address)?.entry,
          open: (signal) => runtime.open(nested.entry, { assets: this.assets, signal }),
        });
    }
    for (const [address, h] of [...this.held]) {
      if (wanted.get(address)?.identity === h.identity && !h.failed && !this.lost(address, h)) continue;
      this.release(address, h);
    }
    const started: Promise<void>[] = [];
    for (const [address, w] of wanted) {
      if (this.held.has(address)) continue;
      const present = this.store.coverage(address);
      // attached by another consumer of a shared store (not an overlay's): left alone
      if (present !== undefined && this.shared.attached.get(address) !== present) continue;
      const h: Held = { identity: w.identity, abort: new AbortController(), opening: present === undefined };
      this.held.set(address, h);
      this.shared.claims.set(address, (this.shared.claims.get(address) ?? 0) + 1);
      // another overlay's: claimed, so it stays while this one needs it
      if (present !== undefined) continue;
      started.push(
        w.open(h.abort.signal).then(
          (opened) => {
            h.opening = false;
            if (h.abort.signal.aborted || this.held.get(address) !== h) return opened.dispose();
            // one attached meanwhile (another overlay's opening settled first) is shared, never replaced
            if (this.store.coverage(address) !== undefined) return opened.dispose();
            if (this.store.attach(address, opened)) this.shared.attached.set(address, this.store.coverage(address));
            else {
              h.refused = true;
              opened.dispose();
            }
          },
          (e) => {
            h.opening = false;
            if (isAbortError(e) || h.abort.signal.aborted) return;
            if (this.held.get(address) !== h) return;
            h.failed =
              e instanceof SourceOpenError
                ? e.issues.map((i) => ({ ...i, address }))
                : [
                    issueOf(RenderIssues, 'seq.render.source', `The source could not be opened: ${String(e)}`, {
                      address,
                      hint: RenderIssues.codes['seq.render.source'].hint,
                    }),
                  ];
          },
        ),
      );
    }
    await Promise.all(started);
    return this.issues;
  }

  /** A held source opened (not failed, not refused) and no longer in the store. */
  private lost(address: string, h: Held): boolean {
    return !h.opening && !h.failed && !h.refused && this.store.coverage(address) === undefined;
  }

  /** Lets go of its claim: the last claim detaches what an overlay attached (never another consumer's). */
  private release(address: string, h: Held): void {
    h.abort.abort();
    this.held.delete(address);
    const n = (this.shared.claims.get(address) ?? 1) - 1;
    if (n > 0) {
      this.shared.claims.set(address, n);
      return;
    }
    this.shared.claims.delete(address);
    const attached = this.shared.attached.get(address);
    this.shared.attached.delete(address);
    if (attached !== undefined && attached === this.store.coverage(address)) this.store.detach(address);
  }

  /**
   * Holds in the store the nodes a layout shows (`ViewLayout.shown`), and lets go of those it no longer shows: a node
   * shown stays when the region that brought it is evicted, and goes once nothing shows it.
   */
  hold(nodes: readonly OverlayNode[]): void {
    const now = new Set<string>();
    for (const n of nodes) {
      const k = nodeKey(n.source, n.key);
      now.add(k);
      if (!this.shown.has(k)) this.shown.set(k, this.store.hold(n.source, n.key));
    }
    for (const [k, release] of [...this.shown])
      if (!now.has(k)) {
        this.shown.delete(k);
        release();
      }
  }

  /** Lets go of its holds and its claims: what it alone still needed is detached. */
  dispose(): void {
    this.hold([]);
    for (const [address, h] of [...this.held]) this.release(address, h);
  }
}
