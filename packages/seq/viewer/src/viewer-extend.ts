/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { BindingTable } from '@molstar/seqx-viewer/bindings';
import type { HeaderEntry } from './header-entries.js';
import type { SeqViewer, ViewerCommand } from './viewer.js';
import type { createBehaviors } from './viewer-behaviors.js';
import type { BehaviorContext, BehaviorHandle, SeqViewerExtension, SeqViewerSpec } from './viewer-spec.js';

// An extension registered at run time (_spec/README.md "Composition", Extensions): `SeqViewer.extend`'s lazy chunk.

/** What the viewer keeps for `extend`: the extensions it registered with their counts, the bindings, the context. */
export interface ViewerExtensions {
  readonly extended: Map<SeqViewerExtension, { n: number; readonly undo: () => void; readonly dispose: () => void }>;
  /** Whether the one retire hook for late behaviors is registered. */
  hooked?: boolean;
  /** The viewer's spec as composed (its extensions' parts added). */
  readonly spec: SeqViewerSpec;
  /** The binding table of a spec: the composed one with the late extensions' parts added. */
  readonly bind: (spec: SeqViewerSpec) => BindingTable;
  readonly context: BehaviorContext;
  readonly onRetire: (dispose: () => void) => void;
  /** The viewer's own maker of behaviors (passed, so this chunk shares no module with the viewer's). */
  readonly create: typeof createBehaviors;
}

/** The context's tables moved: the bindings made again, the view resolved with the new kinds (its frame follows). */
function kindsChanged(viewer: SeqViewer, x: ViewerExtensions): void {
  let spec = x.spec;
  for (const e of x.extended.keys()) spec = e.composeViewer?.(spec) ?? spec;
  viewer.bindings = x.bind(spec);
  viewer.view.setKinds(viewer.context.kinds);
  viewer.store.setRuntimes(viewer.context.runtimes);
}

/**
 * Registers `ext` on a viewer: the context's entries (its checks throw first, nothing changed), its behaviors (made
 * after the viewer's, `composed` called), its commands and bindings; the undo takes them away in reverse. The same
 * descriptor again is counted: its parts are made once and go with its last undo (one composed is the viewer's). What
 * a late registration would leave unread is refused, as a format is: a behavior with a part over the canvas (`attach`,
 * called when it is mounted); its settings compose (`_extend` reads the composition's pairs under its namespace). Its header entries go on the viewer's
 * headers (`register`), with it.
 */
export function extendViewer(viewer: SeqViewer, ext: SeqViewerExtension, x: ViewerExtensions): () => void {
  const context = viewer.context;
  let held = x.extended.get(ext);
  if (!held && context.extensions.includes(ext)) return () => {};
  if (!held) {
    if (ext.bindings && !ext.composeViewer)
      throw new Error(`${ext.id}: a viewer extension is made by defineViewerExtension`);
    // settings compose again on `_extend` (the composition's pairs under its namespace, core's "Context", Settings)
    if (ext.registry.some((e) => /^(format|query)/.test(e.kind)))
      throw new Error(`${ext.id}: formats and query entries are composed, not registered late`);
    context._extend([...context.extensions, ext]);
    const disposers: (() => void)[] = [];
    const commands: ViewerCommand[] = [...(ext.commands ?? [])];
    // only the behaviors this call made are its own (one of a name the viewer has throws, and leaves the viewer's)
    const before = new Set(viewer.behaviors.keys());
    const made = () => [...viewer.behaviors.keys()].filter((name) => !before.has(name));
    const remove = () => {
      for (const d of disposers.reverse()) d();
      disposers.length = 0;
      for (const name of made()) viewer.behaviors.delete(name);
      context._extend(context.extensions.filter((e) => e !== ext));
    };
    try {
      x.create(ext.behaviors ?? [], x.context, viewer.behaviors, commands, (d) => disposers.push(d));
      for (const name of made()) {
        const h = viewer.behaviors.get(name) as BehaviorHandle | undefined;
        if (h?.attach)
          throw new Error(`${ext.id}: a behavior with a part over the canvas is composed, not registered late`);
        h?.composed?.(viewer);
      }
    } catch (e) {
      remove();
      throw e;
    }
    viewer.commands.push(...commands);
    // its header entries on the viewer's headers (none composed: nothing reads them), once its commands are there
    // (an entry names one); taken away with it
    const headers = viewer.behavior<{ register(...e: HeaderEntry[]): () => void }>('headers');
    if (ext.headerActions?.length && headers)
      try {
        disposers.push(headers.register(...ext.headerActions));
      } catch (e) {
        remove();
        for (const c of commands) viewer.commands.splice(viewer.commands.indexOf(c), 1);
        throw e;
      }
    held = {
      n: 0,
      dispose: () => disposers.forEach((d) => d()),
      undo: () => {
        remove();
        for (const c of commands) {
          const i = viewer.commands.indexOf(c);
          if (i >= 0) viewer.commands.splice(i, 1);
        }
      },
    };
    x.extended.set(ext, held);
    // the behaviors registered late go with the viewer: one hook for them all, never one a registration
    if (!x.hooked) {
      x.hooked = true;
      x.onRetire(() => {
        for (const h of x.extended.values()) h.dispose();
      });
    }
  }
  held.n++;
  kindsChanged(viewer, x);
  let undone = false;
  return () => {
    if (undone || viewer.retired) return;
    undone = true;
    if (--held.n) return;
    x.extended.delete(ext);
    held.undo();
    kindsChanged(viewer, x);
  };
}
