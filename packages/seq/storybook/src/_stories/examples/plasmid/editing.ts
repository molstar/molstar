/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Patch } from '@molstar/core/repr/patch';
import { type Node, printAddress, resolveAddress } from '@molstar/core/repr/tree';
import { PlasmidLength } from '@molstar/seqx-examples/plasmid/derive';
import type { EditingBehavior } from '@molstar/seqx-viewer/behaviors';
import type { EditingSession } from '../../features/editing/compose.js';

// The Plasmid's editing views, DOM-free (the story and its headless twin call these): the ring edited through the
// viewer's commands and its editing behavior, each edit one patch and one step of the host's stack, so Undo and Redo
// walk them. Hide a strand's track; show the sites of some enzymes only (a filter in the sites layer's data); bookmark
// bla from its start to its stop, across the origin; and, on the rotated linear view, set the start of its range.

const ring = ['views', 'main', 'ring'];
export const strandTrack = (strand: '+' | '-') => `/views/main/ring/${strand === '+' ? 'plus' : 'minus'}`;
const sitesLayer = [...ring, 'sites', 'sites'];
const editing = (s: EditingSession) => s.viewer.behaviors.get('editing') as EditingBehavior;

/** Whether a strand's track is hidden (its own `hidden`). */
export function strandHidden(s: EditingSession, strand: '+' | '-'): boolean {
  const id = strand === '+' ? 'plus' : 'minus';
  return resolveAddress(s.viewer.spec.root.value, [...ring, id])?.params?.hidden === true;
}

/** Hides or shows a strand's whole track: one patch, one step. */
export function toggleStrand(s: EditingSession, strand: '+' | '-'): Promise<string> {
  const command = strandHidden(s, strand) ? 'track.show' : 'track.hide';
  return Promise.resolve(s.viewer.execute(command, { address: strandTrack(strand) }));
}

/** The data of the sites layer for the enzymes named (all of them when none): the layer's `data` address. */
export function sitesData(enzymes?: readonly string[]): string {
  return enzymes ? `/data/sites${printAddress([{ kind: 'item', params: { id: [...enzymes] } }])}` : '/data/sites';
}

/** The sites layer's data as the spec says. */
export function sitesDataOf(s: EditingSession): unknown {
  return resolveAddress(s.viewer.spec.root.value, sitesLayer)?.params?.data;
}

/** Shows the sites of some enzymes only (a filter on the layer's data), or all: a preview that ends in one patch. */
export async function showEnzymes(s: EditingSession, enzymes?: readonly string[]): Promise<boolean> {
  const patch: Patch = { ops: [{ op: 'set', at: sitesLayer, key: 'data', value: sitesData(enzymes) }] };
  if (!editing(s).preview(patch)) return false;
  await editing(s).endPreview(enzymes ? `Show ${enzymes.join(', ')}` : 'Show every enzyme');
  return true;
}

/**
 * Bookmarks bla from its start to its stop: a range over the origin of the rotated variant (1-based, both ends
 * included), so its focus takes the shortest arc. `start` and `end` are the feature's, 0-based and end-exclusive.
 */
export function bookmarkRange(s: EditingSession, start: number, end: number, id: string, label: string) {
  const to = `/data/plasmid/(range :from ${start + 1} :to ${((end - 1) % PlasmidLength) + 1})`;
  return Promise.resolve(s.viewer.execute('refs.add', { ref: { id, to, label } }));
}

/** Sets the rotated linear view's range: its start is the one param the user sets (`length` positions from it). */
export async function setStart(s: EditingSession, start: number, length: number): Promise<boolean> {
  const end = (start + length) % PlasmidLength || PlasmidLength;
  const patch: Patch = { ops: [{ op: 'set', at: ['views', 'main', 'linear'], key: 'range', value: [start, end] }] };
  if (!editing(s).preview(patch)) return false;
  await editing(s).endPreview('Set the start');
  return true;
}

/** The linear section's range as the spec says. */
export function rangeOf(s: EditingSession): unknown {
  return resolveAddress(s.viewer.spec.root.value, ['views', 'main', 'linear'])?.params?.range;
}

/** The spec shown without its data. */
export function printedWithoutData(root: Node): Node {
  return { ...root, children: (root.children ?? []).filter((c) => c.id !== 'data') };
}
