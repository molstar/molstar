/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeProbe } from '@molstar/prism-ui/dom/scope';
import { ReactiveProbe } from '@molstar/prism-ui/reactive/state';

// The waste audit (view-models.md §3, a view follows the narrowest source that decides its output): per state change,
// the listener calls it took, the derived values it recomputed (wasted when the value came out the same), and the
// binding writes and `watch` renders that followed from it (wasted when the DOM did not change: a MutationObserver's
// records, taken right before and after each). Put down to the state that changed and to the view that did the work,
// both by where they were made (the first stack frame outside the runtime).

/** Frames that are the runtime's, not a view's or a model's. */
const Runtime =
  /\/prism\/ui\/src\/(dom\/(scope|render|control-flow|jsx-runtime|jsx-dev-runtime|element|view|writes)|reactive\/)|node_modules|\/_test\/|<anonymous>|node:/;

/** Where something was made: the first frame outside the runtime, as `package/path:line`. */
export function siteOf(stack: string | undefined): string {
  for (const line of (stack ?? '').split('\n').slice(1)) {
    const m = /\(?((?:\/|file:)[^()\s]+):(\d+):\d+\)?\s*$/.exec(line);
    if (!m || Runtime.test(m[1])) continue;
    const path = m[1].replace(/^file:\/\//, '');
    const at = path.indexOf('/packages/');
    const ext = path.indexOf('/extensions/');
    const short = at >= 0 ? path.slice(at + 10) : ext >= 0 ? path.slice(ext + 1) : path;
    return `${short.replace('/src/', '/')}:${m[2]}`;
  }
  return '(runtime)';
}

export interface SourceRow {
  source: string;
  changes: number;
  listenerCalls: number;
  derives: number;
  derivesSame: number;
  writes: number;
  writesNoDom: number;
  renders: number;
  rendersNoDom: number;
}

export interface ViewRow {
  view: string;
  calls: number;
  wasted: number;
}

export interface Audit {
  /** Starts counting (what happened before, mounting, is not counted). */
  start(): void;
  stop(): void;
  sources(): SourceRow[];
  views(): ViewRow[];
  /** Source state -> view. */
  pairs(): ViewRow[];
  /** Wasted work in all: recomputes that changed nothing, writes and renders that wrote nothing. */
  wasted(): number;
  dispose(): void;
}

/**
 * Installs the probes: call before anything to be audited is made (states and views are named where they are made).
 * `roots`: what the views write into, observed for mutations (a shadow root is not seen from the document).
 */
export function installAudit(): Audit & { observe(root: Node): void } {
  const made = new WeakMap<object, string>();
  const bySource = new Map<string, SourceRow>();
  const byView = new Map<string, ViewRow>();
  let counting = false;
  const observer = new MutationObserver(() => {});
  const name = (o: object | undefined) => (o ? (made.get(o) ?? '(unnamed)') : '(first write)');
  const row = (source: string) => {
    let r = bySource.get(source);
    if (!r)
      bySource.set(
        source,
        (r = {
          source,
          changes: 0,
          listenerCalls: 0,
          derives: 0,
          derivesSame: 0,
          writes: 0,
          writesNoDom: 0,
          renders: 0,
          rendersNoDom: 0,
        }),
      );
    return r;
  };
  const byPair = new Map<string, ViewRow>();
  const view = (v: string, wasted: boolean, source?: string) => {
    let r = byView.get(v);
    if (!r) byView.set(v, (r = { view: v, calls: 0, wasted: 0 }));
    r.calls++;
    if (wasted) r.wasted++;
    const key = `${source ?? '?'} -> ${v}`;
    let p = byPair.get(key);
    if (!p) byPair.set(key, (p = { view: key, calls: 0, wasted: 0 }));
    p.calls++;
    if (wasted) p.wasted++;
  };
  ReactiveProbe.on = true;
  ReactiveProbe.created = (o) => made.set(o, siteOf(new Error().stack));
  ReactiveProbe.told = (state, calls) => {
    if (!counting) return;
    const r = row(name(state));
    r.changes++;
    r.listenerCalls += calls;
  };
  ReactiveProbe.derive = (d, changed) => {
    if (!counting) return;
    const r = row(name(ReactiveProbe.source));
    r.derives++;
    if (!changed) r.derivesSame++;
    view(made.get(d) ?? '(derive)', !changed, name(ReactiveProbe.source));
  };
  RuntimeProbe.run = (kind, site, because, fn) => {
    // A first write or render (a view made, a row re-pointed) writes into DOM that may not be attached yet: not waste.
    if (!counting || because.every((b) => b === undefined)) return fn();
    observer.takeRecords();
    fn();
    const wrote = observer.takeRecords().length > 0;
    const r = row(name(because.find((b) => b !== undefined)));
    if (kind === 'write') {
      r.writes++;
      if (!wrote) r.writesNoDom++;
    } else {
      r.renders++;
      if (!wrote) r.rendersNoDom++;
    }
    view(siteOf(site), !wrote, name(because.find((b) => b !== undefined)));
  };
  return {
    observe(root) {
      observer.observe(root, { subtree: true, childList: true, attributes: true, characterData: true });
    },
    start() {
      bySource.clear();
      byView.clear();
      byPair.clear();
      counting = true;
    },
    stop() {
      counting = false;
    },
    sources: () =>
      [...bySource.values()].sort(
        (a, b) => b.derivesSame + b.writesNoDom + b.rendersNoDom - (a.derivesSame + a.writesNoDom + a.rendersNoDom),
      ),
    views: () => [...byView.values()].sort((a, b) => b.wasted - a.wasted),
    pairs: () => [...byPair.values()].sort((a, b) => b.wasted - a.wasted),
    wasted: () => {
      let n = 0;
      for (const r of bySource.values()) n += r.derivesSame + r.writesNoDom + r.rendersNoDom;
      return n;
    },
    dispose() {
      counting = false;
      observer.disconnect();
      ReactiveProbe.on = false;
      ReactiveProbe.created = ReactiveProbe.derive = ReactiveProbe.told = undefined;
      RuntimeProbe.run = undefined;
    },
  };
}
