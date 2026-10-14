/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Introduction's two renderings of the example, linked: the page mounts the viewer and Nightingale as two parts
// with prose between them, so they meet here. `exampleLink` holds the viewer's context and Nightingale's side, and links
// them (support/seqstar-link.ts) while both are there and draw one set. The snippet the page shows
// (minimal-example.tsx) is unchanged: Nightingale mounts here by the same steps as its `mountNightingaleExample`, on the
// viewer's set (a loci names its set) and with `mount`'s events.

import { createRegistries, KernelRegistry, register, resolveView } from '@molstar/seq-core/format';
import type { Issue, SeqSet } from '@molstar/seq-core/model';
import { materialize, NightingaleSupported } from '@molstar/seq-nightingale';
import { mount, type NightingaleEvent } from '@molstar/seq-nightingale/mount';
import type { SeqViewerContext } from '@molstar/seq-viewer';
import { linkNightingale, type NightingaleSide } from '../support/seqstar-link.js';
import { doc } from './minimal-example.js';

const registries = createRegistries();
register(registries, [KernelRegistry]);

/** The example, mounted in `host` on `set` (the viewer's), with its side for the link; `issues`: what it cannot draw. */
export function mountNightingaleSide(host: HTMLElement, set: SeqSet) {
  const resolved = resolveView(doc, set, 'tracks', registries, NightingaleSupported);
  const { root, issues } = materialize(resolved.view, set);
  const events = new Set<(e: NightingaleEvent) => void>();
  const mounted = mount(root, host, { view: resolved.view, set }, (e) => {
    for (const listener of events) listener(e);
  });
  const manager = host.querySelector('nightingale-manager');
  const side: NightingaleSide = {
    view: resolved.view,
    set,
    highlight: (loci) => mounted.highlight(loci),
    // As `mount` copies a range within a scroll group: a `change` on the manager itself (an attribute set on it would
    // be re-applied over every later zoom, NG-01). It does not bubble, so `onRange` never hears it.
    setRange: (start, end) =>
      manager?.dispatchEvent(new CustomEvent('change', { detail: { 'display-start': start, 'display-end': end } })),
    onRange(listener) {
      const heard = (e: Event) => {
        const d = (e as CustomEvent).detail as Record<string, unknown> | null;
        if (d && 'display-start' in d) listener(Number(d['display-start']), Number(d['display-end']));
      };
      host.addEventListener('change', heard);
      return () => host.removeEventListener('change', heard);
    },
    onEvent(listener) {
      events.add(listener);
      return () => events.delete(listener);
    },
  };
  return {
    side,
    dispose: () => mounted.dispose(),
    /** Settles once the elements have their data (the link waits for it: Nightingale's first events are its drawing's). */
    ready: mounted.ready.then((failed): Issue[] => [...resolved.issues, ...issues, ...failed]),
  };
}

/** The viewer's context and Nightingale's side, linked while both are there on one set. */
class ExampleLink {
  private ctx: SeqViewerContext | null = null;
  private side: NightingaleSide | null = null;
  private unlink: (() => void) | null = null;
  private linked: [SeqViewerContext, NightingaleSide] | null = null;
  private watch: { unsubscribe(): void } | null = null;
  private readonly listeners = new Set<(set: SeqSet | null) => void>();

  /** The set the viewer shows now; null: none yet. */
  get set(): SeqSet | null {
    return this.ctx?.graph.assembled.value?.set ?? null;
  }
  /** Calls `listener` with the viewer's set now and whenever it changes. Returns the stop. */
  onSet(listener: (set: SeqSet | null) => void): () => void {
    this.listeners.add(listener);
    listener(this.set);
    return () => this.listeners.delete(listener);
  }
  setViewer(ctx: SeqViewerContext | null) {
    this.watch?.unsubscribe();
    this.ctx = ctx;
    // The replay is the current set: listeners hear it once, here.
    this.watch =
      ctx?.graph.assembled.subscribe(() => {
        this.relink();
        for (const l of this.listeners) l(this.set);
      }) ?? null;
    if (!ctx) {
      this.relink();
      for (const l of this.listeners) l(null);
    }
  }
  setSide(side: NightingaleSide | null) {
    this.side = side;
    this.relink();
  }
  /** Links the pair there is, once: a re-assembly of the same set (or a repeated call) changes nothing. */
  private relink() {
    const ctx = this.ctx;
    const side = this.side;
    const pair = ctx && side && side.set === this.set && !ctx.disposed ? ([ctx, side] as const) : null;
    if (pair && this.linked?.[0] === pair[0] && this.linked[1] === pair[1]) return;
    this.unlink?.();
    this.unlink = null;
    this.linked = pair && [...pair];
    if (pair) this.unlink = linkNightingale(...pair);
  }
}
export const exampleLink = new ExampleLink();
