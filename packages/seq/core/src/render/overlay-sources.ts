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
  type OpenedSource,
  SourceOpenError,
} from '@molstar/seqx-core/data';
import type { Overlay } from '@molstar/seqx-core/resolve';
import { type ResolvedView, templateNeeds } from '@molstar/seqx-core/view';
import { RenderIssues } from './issues.js';

// The overlay's sources (_spec/README.md "Render", The overlay's sources): what the overlay needs read goes through the
// store like any source. A dynamic annotation's inline rows are a `seq` source under the node's address; a nested
// source is opened only when a template of the view needs its depth (the overlay lists it). Each is kept while its
// node is the same entry, and detached when its node goes, so a node's data lives as long as the node.

interface Held {
  /** The node entry it was opened for: a node that comes back is opened again. */
  readonly identity: unknown;
  readonly abort: AbortController;
  attached: boolean;
}

export class OverlaySources {
  private readonly held = new Map<string, Held>();

  constructor(
    readonly store: DataStore,
    private readonly assets: AssetReader,
  ) {}

  /**
   * Opens what the overlay needs now and lets go of what it no longer lists. Resolves, with the issues of what could
   * not be opened, when the openings this call started have settled.
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
      if (wanted.get(address)?.identity === h.identity && (!h.attached || this.store.coverage(address))) continue;
      this.release(address, h);
    }
    const issues: Issue[] = [];
    const started: Promise<void>[] = [];
    for (const [address, w] of wanted) {
      // one held, or attached by another consumer of a shared store
      if (this.held.has(address) || this.store.coverage(address)) continue;
      const h: Held = { identity: w.identity, abort: new AbortController(), attached: false };
      this.held.set(address, h);
      started.push(
        w.open(h.abort.signal).then(
          (opened) => {
            if (h.abort.signal.aborted || this.held.get(address) !== h) return opened.dispose();
            h.attached = this.store.attach(address, opened);
            if (!h.attached) opened.dispose();
          },
          (e) => {
            if (isAbortError(e) || h.abort.signal.aborted) return;
            if (e instanceof SourceOpenError) issues.push(...e.issues.map((i) => ({ ...i, address })));
            else
              issues.push(
                issueOf(RenderIssues, 'seq.render.source', `The source could not be opened: ${String(e)}`, {
                  address,
                  hint: RenderIssues.codes['seq.render.source'].hint,
                }),
              );
          },
        ),
      );
    }
    await Promise.all(started);
    return issues;
  }

  private release(address: string, h: Held): void {
    h.abort.abort();
    this.held.delete(address);
    if (h.attached) this.store.detach(address);
  }

  /** Detaches everything it attached. */
  dispose(): void {
    for (const [address, h] of [...this.held]) this.release(address, h);
  }
}
