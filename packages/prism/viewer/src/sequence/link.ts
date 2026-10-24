/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { escapeHtml } from '@molstar/core/util/string';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { Extents } from '@molstar/seq-core/model';
import { StructureLink } from './molstar/link.js';
import { coordinatesOf, linkViews, type LinkViewsOptions, type StructureView3D } from './molstar/link-views.js';
import type { Mark, MarksView } from '../sources/marks.js';
import { type HoverLabel, lociSize } from '@molstar/seq-viewer';
import type { SequenceStructure, SequenceViewer } from './model.js';

/** The 3D view the sequence links to: Seq*'s contract and, optionally, the camera's reset to the whole structure. */
export interface SequenceView3D extends StructureView3D {
  /** The default framing of the structure (below the sequence band), animated; absent: the camera is left alone. */
  resetCamera?(): void;
  /**
   * Adds `read()` to the 3D label of what a mapped sequence hover lights (Mol*'s label providers, sequence.md §10):
   * the hovered track and its lines under Mol*'s label of the residue. `refresh` labels the highlight again when the
   * lines changed and the highlight did not (another track over the same residue). Absent: the label is Mol*'s alone.
   */
  annotateHover?(read: () => string | undefined): { refresh(): void; dispose(): void };
}

/** Selection origins of the person's own moves in the sequence view (the rest: 3D, data, API, replays). */
const FromSequence: ReadonlySet<string> = new Set(['pointer', 'keyboard']);

/**
 * The two-way link between a viewer and the 3D view of the structure it shows (sequence.md §4): hover and
 * selection both ways, focus. Made when the viewer has a structure and a 3D view of it, dropped and made again when
 * the structure changes. Returns the unlink: the control calls it when it goes, so nothing stays subscribed. The replay
 * rule makes remounting safe (linking sends the current selection as a selection, never a clear). A selection the person
 * empties in the sequence view (a click or Shift+drag removing the last stretch, Escape) also resets the 3D camera, once
 * (sequence.md §4); an emptying from 3D, the selection bar, the API or a data edit leaves the camera.
 */
export function linkSequence(
  viewer: SequenceViewer,
  view3d: (structure: Structure) => SequenceView3D | undefined,
  options?: LinkViewsOptions & { readonly marks?: MarksView },
): () => void {
  const seq = viewer.seq;
  let link: StructureLink | undefined;
  let shownNow: SequenceStructure | undefined;
  let stopViews: (() => void) | undefined;
  let lit: SequenceView3D | undefined;
  /** `clear`: what the sequence hovered in 3D goes with the views; kept on a frame step, where the views start again. */
  const stop = (clear = true) => {
    stopViews?.();
    stopViews = undefined;
    // (the 3D view also drops `hoverMapped`)
    if (clear) lit?.highlight(null);
    lit = undefined;
  };
  const drop = () => {
    stop();
    link?.detach();
    link = shownNow = undefined;
  };
  /** The views over the link, for the 3D view of `structure` (the one shown, or a frame of it: no `replay`). */
  const start = (structure: Structure, view: SequenceView3D, replay = true) => {
    if (!link) return;
    lit = view;
    const stopLink = linkViews(seq, view, link, options, replay);
    const stopMarks = options?.marks ? showMarks(seq, link, structure, options.marks) : undefined;
    const stopLabel = view.annotateHover ? annotateHover(seq, view.annotateHover) : undefined;
    // Whether the selection held something: the reset follows only the move that empties it, not an empty selection again.
    let held = false;
    const sub = seq.selection.current.subscribe((sel) => {
      const parts = [sel.seq, sel.aligned, sel.columns, sel.rows, sel.items, sel.boundaries] as Parameters<
        typeof lociSize
      >[0][];
      const empty = parts.every((l) => {
        const size = lociSize(l);
        return !size || (size.items === 0 && size.length === 0);
      });
      const emptied = held && empty;
      held = !empty;
      if (emptied && FromSequence.has(sel.origin)) view.resetCamera?.();
    });
    stopViews = () => {
      sub.unsubscribe();
      stopLink();
      stopMarks?.();
      stopLabel?.();
    };
  };
  const sub = viewer.structure.subscribe((shown) => {
    drop();
    const set = seq.graph.assembled.value?.set;
    if (!shown || !set) return;
    const structure = viewer.frame.value ?? shown.structure;
    const view = view3d(structure);
    if (!view) return;
    shownNow = shown;
    link = StructureLink.attach(seq.systems, '3d', structure, set, shown.complex);
    start(structure, view);
  });
  // A frame step (sequence.md §2): the link follows the new Structure, its systems stay; only the views (which map 3D
  // loci onto the structure drawn) start again, with no clear of the 3D hover and no replay of the selection, which 3D
  // holds (a replay was a recorded `select.set` per step: perf-audit-2 10/B6).
  const frames = viewer.frame.subscribe((structure) => {
    if (!link || !shownNow || !structure || structure === link.structure) return;
    const view = view3d(structure);
    if (!view || !link.rebind(structure)) return;
    stop(false);
    start(structure, view, false);
  });
  return () => {
    sub.unsubscribe();
    frames.unsubscribe();
    drop();
  };
}

/**
 * A view's marks and hover in the sequence (marks.md): each mark a highlight channel (`marks:<owner>:<id>`) in its
 * tone's style (Seq*'s `highlight` or `accent`, the skin's colors, as the 3D view paints them), on the residues it has in
 * the structure shown; the hover as the 3D view's hover is (so it is not sent back to the 3D view, which lit it).
 * Returns the stop, which clears what it set.
 */
function showMarks(seq: SequenceViewer['seq'], link: StructureLink, structure: Structure, marks: MarksView) {
  const scratch = Extents.create(16);
  const coordsOf = (loci: StructureElement.Loci) => {
    // Another structure's loci would remap by unit id and atom index onto whatever residues share them.
    if (loci.structure.root !== structure.root) return undefined;
    const here = StructureElement.Loci.remap(loci, structure);
    if (StructureElement.Loci.isEmpty(here)) return undefined;
    scratch.count = 0;
    return coordinatesOf(link, link.fromLoci(here, scratch));
  };
  let channels = new Set<string>();
  // The marks already shown, by what they are (an owner's same marks set again are not new); the first list is what was
  // there when the link was made: shown, not focused.
  const keyOf = (m: Mark) => `${m.owner}\u0000${m.id}\u0000${m.structure}\u0000${m.label}`;
  let shown: ReadonlySet<string> | undefined;
  const marked = marks.marks.subscribe((list) => {
    // New marks that focus the sequence: their owner's marks, together, in view (a box's two ranges both shown).
    const focusing = new Set(
      shown ? list.filter((m) => m.focusSequence && !shown!.has(keyOf(m))).map((m) => m.owner) : [],
    );
    shown = new Set(list.map(keyOf));
    if (focusing.size) {
      // One target of their runs, each in its system (a box's ranges are often in two chains).
      const systems: string[] = [];
      const runs: number[] = [];
      for (const m of list) {
        if (!m.focusSequence || !focusing.has(m.owner)) continue;
        const coords = coordsOf(m.loci);
        if (!coords) continue;
        for (let i = 0; i < coords.runs.length; i += 3) {
          const key = coords.systems[coords.runs[i]];
          let k = systems.indexOf(key);
          if (k < 0) k = systems.push(key) - 1;
          runs.push(k, coords.runs[i + 1], coords.runs[i + 2]);
        }
      }
      if (runs.length) seq.focus.fromLoci({ systems, runs: Int32Array.from(runs) }, { origin: 'host:3d' });
    }
    const next = new Set<string>();
    for (const m of list) {
      const coords = coordsOf(m.loci);
      if (!coords) continue;
      const channel = `marks:${m.owner}:${m.id}`;
      seq.highlight.set(channel, coords, { style: m.tone });
      next.add(channel);
    }
    for (const c of channels) if (!next.has(c)) seq.highlight.clear(c);
    channels = next;
  });
  let hovering = false;
  const hovered = marks.hover.subscribe((h) => {
    const coords = h ? coordsOf(h.loci) : undefined;
    if (!coords && !hovering) return;
    hovering = !!coords;
    seq.hover.setLoci(coords ?? null, '3d');
  });
  return () => {
    marked.unsubscribe();
    hovered.unsubscribe();
    for (const c of channels) seq.highlight.clear(c);
    if (hovering) seq.hover.setLoci(null, '3d');
  };
}

/**
 * Seq*'s label of a mapped hover as the extra lines of Mol*'s label (sequence.md §10): the track's name, then Seq*'s
 * own lines (an annotation, its fields, a value), never its title or the residue's lines (Mol* names the residue). Escaped: the
 * text comes from loaded documents and Mol*'s labels are HTML. The lines are made once per Seq* label (edge rate) and
 * read by the provider at Mol*'s label time; a new Seq* label over the same highlight labels it again.
 */
function annotateHover(seq: SequenceViewer['seq'], annotate: NonNullable<SequenceView3D['annotateHover']>) {
  let html: string | undefined;
  const handle = annotate(() => html);
  const sub = seq.labels.current.subscribe((label) => {
    const next = label ? sequenceLabelHtml(label) : undefined;
    if (next === html) return;
    html = next;
    handle.refresh();
  });
  return () => {
    sub.unsubscribe();
    handle.dispose();
  };
}

/** Escaped for Mol*'s label HTML, on one line: a line break would make the label Markdown to Prism (`hoverRich`). */
const oneLine = (s: string) => escapeHtml(s.replace(/\s+/g, ' '));

/**
 * Seq*'s providers that label the residue itself (its name and numbers, its copies, its stored letter): Mol*'s label
 * says that, so their lines are left out (sequence.md §10). What stays is the track's own: an annotation, a value; and
 * `copies: varies` (the copies differ here), which Mol*'s label of one residue cannot say.
 */
const ResidueProviders: ReadonlySet<string> = new Set(['position', 'copies']);
const residueLine = (from: string | undefined, l: { readonly key: string; readonly value: string }) =>
  !!from && ResidueProviders.has(from) && !(from === 'copies' && l.value === 'varies');

/**
 * A Seq* label's track and its own lines as Mol*'s label HTML (`<b>track</b>`, then `key: value` a line): none of the
 * residue's lines, and no track line when a line already names it (a value keyed by the track's annotation). Nothing
 * of the track's own: nothing added.
 */
export function sequenceLabelHtml(label: HoverLabel): string | undefined {
  const lines: string[] = [];
  let named = false;
  for (let i = 0; i < label.lines.length; i++) {
    const l = label.lines[i];
    const from = label.from?.[i];
    if (residueLine(from, l)) continue;
    if (label.track && l.key === label.track) named = true;
    lines.push(oneLine(l.value ? `${l.key}: ${l.value}` : l.key));
  }
  if (!lines.length) return undefined;
  if (label.track && !named) lines.unshift(`<b>${oneLine(label.track)}</b>`);
  return lines.join('<br/>');
}
