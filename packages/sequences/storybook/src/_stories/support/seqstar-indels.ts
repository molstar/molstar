/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Case 15 (indels, `core/_test/format/cases.ts`) in the Viewer cases story: the selections its Select choices make, and
// what each lights per section of each view, written out by hand (not derived from the case's own helpers), for the
// smoke test's steps.

import { ChannelMarks, MarkBit, type SeqViewerContext } from '@molstar/seq-viewer';

/** Case 15's views: the wild type's axis (the engineer's frame), the alignment, each design on its own axis. */
export const INDEL_VIEWS = { wt: 'wt_frame', aligned: 'aligned', own: 'own_axes' } as const;

/** What the story's selection buttons select, as `seqv:select` targets. */
export const INDEL_SELECTIONS = [
  { id: 'h102', label: 'H102 (wild type)', ref: 'wt', start: 101, stop: 102 },
  { id: 'graft', label: "D3's graft", ref: 'd3', start: 91, stop: 97 },
  { id: 'tag', label: "D4's tag", ref: 'd4', start: 108, stop: 116 },
] as const;
export type IndelSelection = (typeof INDEL_SELECTIONS)[number];

export const selectCommand = (s: IndelSelection) => ({
  command: 'seqv:select',
  params: { target: { refs: [s.ref], args: { start: s.start, stop: s.stop } } },
});

/**
 * What each selection lights per section of each view (seqspec:model#projection), by hand; a caret (a boundary lit) as
 * `[at, at]`. H102 lights in every design that keeps it, at its own index (D2 lost four residues before it, D3 two and
 * gained six), and in D4, which deleted it, the caret at the junction its deletion leaves; on the alignment its
 * column. D3's graft lights in D3, and as a caret where it sits on the wild type and the other designs (after WT 93);
 * on the alignment its insert columns. D4's tag, beyond what the mappings pair, lights only in D4 and its columns.
 */
export const INDEL_LIT: Readonly<Record<string, Readonly<Record<IndelSelection['id'], Record<string, number[][]>>>>> = {
  [INDEL_VIEWS.own]: {
    h102: { wt: [[101, 102]], d1: [[101, 102]], d2: [[97, 98]], d3: [[105, 106]], d4: [[101, 101]] },
    graft: { wt: [[93, 93]], d1: [[93, 93]], d2: [[89, 89]], d3: [[91, 97]], d4: [[94, 94]] },
    tag: { wt: [], d1: [], d2: [], d3: [], d4: [[108, 116]] },
  },
  [INDEL_VIEWS.aligned]: { h102: { main: [[108, 109]] }, graft: { main: [[94, 100]] }, tag: { main: [[117, 125]] } },
  [INDEL_VIEWS.wt]: { h102: { main: [[101, 102]] }, graft: { main: [[93, 93]] }, tag: { main: [] } },
};

/** What the current selection lights per section of the current view (runs, then carets as `[at, at]`), by section. */
export function litPerSection(ctx: SeqViewerContext): Record<string, number[][]> {
  const view = ctx.graph.assembled.value?.view?.resolved;
  const set = ctx.graph.assembled.value?.set;
  const out: Record<string, number[][]> = {};
  if (!view || !set) return out;
  const marks = new ChannelMarks(ctx);
  try {
    for (let s = 0; s < view.sections.count; s++) {
      const m = marks.marks(s)?.runs(MarkBit.selection);
      const runs: number[][] = [];
      for (let r = m ? m.from(0) : 0; m && r < m.to(0); r++) runs.push([m.start[r], m.stop[r]]);
      // Its carets (a boundary lit: an insert seen from the other side of its map), as empty runs.
      for (let k = 0; m && k < m.carets; k++) runs.push([m.caretAt[k], m.caretAt[k]]);
      out[set.strings.get(view.sections.ref[s]).replace(/^[^.]*\./, '')] = runs;
    }
  } finally {
    marks.dispose();
  }
  return out;
}
