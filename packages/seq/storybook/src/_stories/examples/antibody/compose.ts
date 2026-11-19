/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Patch } from '@molstar/core/repr/patch';
import { cdrAddress, wildType } from '@molstar/seqx-examples/antibody/data';
import type { Chain } from '@molstar/seqx-examples/antibody/numbering';
import { RulerLayer, rulerLabel } from '@molstar/seqx-examples/antibody/schemes';
import { findNumbering } from '@molstar/seqx-core/query';
import type { EditingBehavior } from '@molstar/seqx-viewer/behaviors/editing';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { EditingSession } from '../../features/editing/compose.js';

// The Antibody's host panel and editing views, DOM-free (the stories and their headless twins call these): the hovered
// residue's label in both numberings, read from the viewer's query model (the source's labels once they are in);
// which scheme a ruler shows (one patch); a CDR selected by its address in Kabat's labels, bookmarked from the
// selection, renamed, recolored and removed, each one step of the host's stack.

const editing = (s: EditingSession) => s.viewer.behaviors.get('editing') as EditingBehavior;

/** One patch as one step: a preview ended at once. */
async function step(s: EditingSession, label: string, patch: Patch): Promise<string> {
  if (!editing(s).preview(patch)) return 'Refused';
  await editing(s).endPreview(label);
  return label;
}

/** A position's label in each numbering of its system (`''` where it has none, or the labels are not in yet). */
export function labelsAt(viewer: SeqViewer, system: string, position: number): { kabat: string; imgt: string } {
  const s = viewer.addresses.queryModel().system(system);
  const at = (id: string) => {
    const n = s ? findNumbering(s, id) : undefined;
    return n?.numbering.labelAt(position - n.shift) ?? '';
  };
  return { kabat: at('kabat'), imgt: at('imgt') };
}

/** The panel's lines: the hovered residue in both numberings, and the selection's address. */
export function readoutLines(viewer: SeqViewer): string[] {
  const h = viewer.hover.current;
  const lines: string[] = [];
  if (h.system !== undefined && Number.isFinite(h.position)) {
    const { kabat, imgt } = labelsAt(viewer, h.system, h.position);
    const chain = h.system.split('/').pop();
    lines.push(`Hover: ${chain} ${h.position + 1}, Kabat ${kabat || '–'}, IMGT ${imgt || '–'}`);
  } else lines.push('Hover: (a residue)');
  lines.push(`Selection: ${viewer.addresses.selectionAddress() ?? '(none)'}`);
  return lines;
}

/** The panel's lines with the bookmarks listed: each its id, label, color and what it names. */
export function bookmarkLines(viewer: SeqViewer): string[] {
  const refs = viewer.addresses.refs.value;
  return [
    ...readoutLines(viewer),
    `Bookmarks (${refs.length}):`,
    ...(refs.length
      ? refs.map((r) => `  ${r.label} [${r.id}]${r.color ? ` ${r.color}` : ''}: ${r.to.join(', ')}`)
      : ['  (none)']),
  ];
}

const SchemeNames: Readonly<Record<string, string>> = { kabat: 'Kabat', imgt: 'IMGT', ordinal: 'Ordinal' };

/** Which numbering the chains' rulers show: one patch of each ruler layer's `numbering`, one step. */
export function setRuler(
  s: EditingSession,
  numbering: string,
  chains: readonly Chain[] = ['light', 'heavy'],
): Promise<string> {
  // the track's label says the scheme its ruler shows
  const scheme = SchemeNames[numbering] ?? numbering;
  return step(s, `Ruler in ${numbering}`, {
    ops: chains.flatMap((c) => [
      { op: 'set' as const, at: RulerLayer(c), key: 'numbering', value: numbering },
      { op: 'set' as const, at: RulerLayer(c).slice(0, -1), key: 'label', value: rulerLabel(c, scheme) },
    ]),
  });
}

/** Selects a CDR of the wild type by its address in Kabat's labels (session state). */
export function selectCdr(s: EditingSession, chain: Chain, cdr: 1 | 2 | 3): string {
  return s.viewer.addresses.select(cdrAddress(chain, cdr, 'kabat', wildType()[chain].kabat)).message;
}

/** Bookmarks the selection (the selection's bookmark command): one step. */
export const bookmark = (s: EditingSession, id: string, label: string) =>
  Promise.resolve(s.viewer.execute('selection.bookmark', { id, label }));

/** Renames a bookmark (`refs.label`): one step. */
export const rename = (s: EditingSession, id: string, label: string) =>
  Promise.resolve(s.viewer.execute('refs.label', { address: `/refs/${id}`, label }));

/** Recolors a bookmark (its `color`): one patch, one step. */
export const recolor = (s: EditingSession, id: string, color: string) =>
  step(s, 'Recolor a bookmark', { ops: [{ op: 'set', at: ['refs', id], key: 'color', value: color }] });

/** Removes a bookmark (`refs.remove`): one step. */
export const remove = (s: EditingSession, id: string) =>
  Promise.resolve(s.viewer.execute('refs.remove', { address: `/refs/${id}` }));

/** The bookmarks as the spec holds them: id, label, color. */
export const bookmarks = (s: EditingSession) =>
  s.viewer.addresses.refs.value.map((r) => ({ id: r.id, label: r.label, color: r.color }));
