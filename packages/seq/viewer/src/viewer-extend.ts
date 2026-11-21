/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { BindingTable } from '@molstar/seqx-viewer/bindings';
import type { SeqViewer, ViewerCommand } from './viewer.js';
import type { createBehaviors } from './viewer-behaviors.js';
import type { BehaviorContext, BehaviorHandle, SeqViewerExtension } from './viewer-spec.js';

// An extension registered at run time (_spec/README.md "Composition", Extensions): `SeqViewer.extend`'s lazy chunk.

/** What the viewer keeps for `extend`: the extensions it registered with their counts, the bindings, the context. */
export interface ViewerExtensions {
  readonly extended: Map<SeqViewerExtension, { n: number; readonly undo: () => void }>;
  /** The binding table with the spec's bindings and the composed and late extensions' layers over them. */
  readonly bind: (late: readonly SeqViewerExtension[]) => BindingTable;
  readonly context: BehaviorContext;
  readonly onRetire: (dispose: () => void) => void;
  /** The viewer's own maker of behaviors (passed, so this chunk shares no module with the viewer's). */
  readonly create: typeof createBehaviors;
}

/** The context's tables moved: the bindings made again, the view resolved with the new kinds (its frame follows). */
function kindsChanged(viewer: SeqViewer, x: ViewerExtensions): void {
  viewer.bindings = x.bind([...x.extended.keys()]);
  viewer.view.setKinds(viewer.context.kinds);
}

/**
 * Registers `ext` on a viewer: the context's entries (its checks throw first, nothing changed), its behaviors (made
 * after the viewer's, `composed` called), its commands and bindings; the undo takes them away in reverse. The same
 * descriptor again is counted: its parts are made once and go with its last undo (one composed is the viewer's).
 */
export function extendViewer(viewer: SeqViewer, ext: SeqViewerExtension, x: ViewerExtensions): () => void {
  const context = viewer.context;
  let held = x.extended.get(ext);
  if (!held && context.extensions.includes(ext)) return () => {};
  if (!held) {
    // a format runtime or a query entry changes what the open spec resolves to: composed, never registered late
    if (ext.registry.some((e) => /^(format|query)/.test(e.kind)))
      throw new Error(`${ext.id}: formats and query entries are composed, not registered late`);
    context._extend([...context.extensions, ext]);
    const disposers: (() => void)[] = [];
    const commands: ViewerCommand[] = [...(ext.commands ?? [])];
    const names = (ext.behaviors ?? []).map((b) => b.provider.name);
    const remove = () => {
      for (const d of disposers.reverse()) d();
      disposers.length = 0;
      for (const name of names) viewer.behaviors.delete(name);
      context._extend(context.extensions.filter((e) => e !== ext));
    };
    try {
      x.create(ext.behaviors ?? [], x.context, viewer.behaviors, commands, (d) => disposers.push(d));
      for (const name of names) (viewer.behaviors.get(name) as BehaviorHandle | undefined)?.composed?.(viewer);
    } catch (e) {
      remove();
      throw e;
    }
    viewer.commands.push(...commands);
    // its behaviors go with the viewer when it is not undone first
    x.onRetire(() => disposers.forEach((d) => d()));
    held = {
      n: 0,
      undo: () => {
        remove();
        for (const c of commands) viewer.commands.splice(viewer.commands.indexOf(c), 1);
      },
    };
    x.extended.set(ext, held);
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
