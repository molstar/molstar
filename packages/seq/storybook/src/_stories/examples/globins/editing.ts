/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { resolveAddress } from '@molstar/core/repr/tree';
import { theme } from '@molstar/seqx-core/lang';
import { copiesPatch, HemoglobinSection } from '@molstar/seqx-examples/globins/hemoglobin';
import type { EditingBehavior } from '@molstar/seqx-viewer/behaviors/editing';
import type { EditingSession } from '../../features/editing/compose.js';

// The Globins' editing views, DOM-free (examples.md, Globins: Views): the story and its headless twin call these on an
// editing session (`composeEditing`: a host's stack). Each edit is one patch, one step: the family's residue scheme and
// its logo hidden; the folds kept (`clades.keep`) and a row pinned; hemoglobin's copies split or merged.

/** The family's rows track and layer, and its logo track. */
export const RowsTrack = '/views/main/family/rows';
export const RowsLayer = `${RowsTrack}/rows`;
export const LogoTrack = '/views/main/family/logo';

const node = (s: EditingSession, address: string) =>
  resolveAddress(s.viewer.spec.root.value, address.split('/').filter(Boolean));

/** The rows' residue scheme: one theme patch. */
export function setScheme(s: EditingSession, scheme: string): Promise<string> {
  return Promise.resolve(s.viewer.execute('layer.color', { address: RowsLayer, color: theme.scheme(scheme) }));
}

/** Whether the logo is hidden (its track's own `hidden`). */
export const logoHidden = (s: EditingSession) => node(s, LogoTrack)?.params?.hidden === true;

/** Hides or shows the logo track: one patch. */
export function toggleLogo(s: EditingSession): Promise<string> {
  return Promise.resolve(s.viewer.execute(logoHidden(s) ? 'track.show' : 'track.hide', { address: LogoTrack }));
}

/** The folds the spec opens with (the layer's `collapsed`). */
export const keptFolds = (s: EditingSession) => (node(s, RowsLayer)?.params?.collapsed as string[] | undefined) ?? [];

/** Folds or unfolds a clade: session state, no step. */
export function toggleClade(s: EditingSession, clade: string): Promise<string> {
  return Promise.resolve(s.viewer.execute('clades.toggle', { track: RowsTrack, clade }));
}

/** Writes the folds shown into the spec: one patch. */
export function keepFolds(s: EditingSession): Promise<string> {
  return Promise.resolve(s.viewer.execute('clades.keep', { track: RowsTrack }));
}

/** Pins (or unpins) data rows above the tree's order: one patch. */
export function pinRows(s: EditingSession, rows: readonly number[]): Promise<string> {
  return Promise.resolve(s.viewer.execute('rows.pin', { track: RowsTrack, rows }));
}

/** The rows pinned (the layer's `pinned`). */
export const pinnedRows = (s: EditingSession) => (node(s, RowsLayer)?.params?.pinned as string[] | undefined) ?? [];

/** Hemoglobin's copies as the spec says them (`merge` when absent). */
export const copiesOf = (s: EditingSession) =>
  (node(s, HemoglobinSection)?.params?.copies as 'merge' | 'each' | undefined) ?? 'merge';

/** Splits or merges hemoglobin's copies: one patch (the section's `copies`). */
export async function toggleCopies(s: EditingSession): Promise<string> {
  const next = copiesOf(s) === 'merge' ? 'each' : 'merge';
  const editing = s.viewer.behaviors.get('editing') as EditingBehavior;
  if (!editing.preview(copiesPatch(next))) return 'Refused';
  await editing.endPreview(next === 'each' ? 'Split the copies' : 'Merge the copies');
  return next === 'each' ? 'Split' : 'Merged';
}
