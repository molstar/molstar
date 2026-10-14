/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// "Show me" across views: focus on extents in the view that holds them. A view holds them when a block of one of its
// sections shows one of their systems directly: its unit or a merged copy, its unit's molecule, or the alignment of a
// row or of the columns. The current view when it does, else the first that does in
// the document's order, switched to; when none does, the current view, where the projection may still reach them
// through links and paths. The other views' blocks are resolved headlessly (`resolveSection`), without drawing them.

import { type BlockRow, childrenOf, type DocIndex, indexDoc, resolveSection } from '@molstar/seq-core/format';
import { type Extents, CoordinateSystemKind as K, type SeqSet } from '@molstar/seq-core/model';
import type { SeqViewerContext } from './context.js';
import type { Origin } from './state.js';

const shows = (set: SeqSet, b: BlockRow, kind: number, id: number) =>
  kind === K.row || kind === K.columns
    ? b.alignment === id
    : b.alignment < 0 &&
      (kind === K.unit ? b.unit === id || (b.members?.includes(id) ?? false) : set.units.molecule[b.unit] === id);

/** Whether a view's blocks show one of `ext`'s systems directly. */
function holds(ctx: SeqViewerContext, index: DocIndex, set: SeqSet, view: string, ext: Extents): boolean {
  const reg = ctx.registry.svs;
  const node = index.nodes.get(view);
  for (const s of node ? childrenOf(index, node) : []) {
    if (s.kind !== 'section') continue;
    for (const b of resolveSection(index, set, s.ref, reg).section.blocks) {
      for (let i = 0; i < ext.count; i++) if (shows(set, b, ext.kind[i], ext.id[i])) return true;
    }
  }
  return false;
}

/** The view that holds `ext`: the current one when it does, else the first that does; null: none, or no document. */
export function viewShowing(ctx: SeqViewerContext, ext: Extents): string | null {
  const doc = ctx.document.shown.value;
  const set = ctx.graph.assembled.value?.set;
  if (!doc || !set) return null;
  const index = indexDoc(doc, ctx.registry.svs);
  const current = ctx.views.current.value;
  if (current && holds(ctx, index, set, current, ext)) return current;
  return ctx.views.list.value.find((v) => v.ref !== current && holds(ctx, index, set, v.ref, ext))?.ref ?? null;
}

/** "Show me": focus on `ext` in the view that holds it (switched to first); false when no section shows any of it. */
export async function showExtents(
  ctx: SeqViewerContext,
  ext: Extents,
  opts: { transitionMs?: number; origin?: Origin } = {},
): Promise<boolean> {
  const view = viewShowing(ctx, ext);
  if (view && view !== ctx.views.current.value) {
    ctx.views.setCurrent(view, opts.origin);
    await ctx.graph.update();
  }
  return ctx.focus.fromExtents(ext, opts);
}
