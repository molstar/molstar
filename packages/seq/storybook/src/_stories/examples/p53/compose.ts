/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { printLink } from '@molstar/seqx-core/query';
import { type EditingSession, proteinSection, Section, track } from '../../features/editing/compose.js';

// The p53 example stories' shared composition, DOM-free: each story and its headless twin call these. The pieces are the
// Features stories' (editing, addressing, selection); what is here is what the examples add on p53. Every action goes
// through the viewer's public API (commands, `spec.patch`, its input), never around it.

/** The id of the track the example adds. */
export const AddedTrack = 'domains-2';

/**
 * Adds a second domains track after the first, a copy of it with its own id and label, through `track.add`: one patch,
 * so one step of the host's stack. A second call finds the id taken and says so.
 */
export async function addTrack(s: EditingSession): Promise<string> {
  const section = proteinSection(s);
  const domains = section?.children?.find((c) => c.id === 'domains');
  if (!section || !domains) return 'No domains track to copy';
  if (section.children?.some((c) => c.id === AddedTrack)) return 'Already added';
  const node: Node = { ...domains, id: AddedTrack, params: { ...domains.params, label: 'Domains (added)' } };
  return Promise.resolve(s.viewer.execute('track.add', { address: Section, node, after: 'domains' }));
}

/** The track ids of the protein section, in order. */
export function trackOrder(s: EditingSession): string[] {
  return (proteinSection(s)?.children ?? []).map((c) => c.id as string);
}

/** A track's height as the spec shown says (a preview included), else as laid out; undefined before a layout. */
export function trackHeight(s: EditingSession, id: string): number | undefined {
  const set = proteinSection(s, true)?.children?.find((c) => c.id === id)?.params?.height;
  if (typeof set === 'number') return set;
  const laid = s.viewer.frame?.layout.sections.flatMap((x) => x.tracks).find((t) => t.track.address === track(id));
  return laid?.box.height;
}

/** A run of the short text: plain, or a link (its text and its `!seq-…` href). */
export type Run = string | { readonly text: string; readonly href: string };

const link = (text: string, address: string): Run => ({
  text,
  href: printLink({ actions: [{ action: 'focus', address }] }),
});

/** The text beside the gene's region: each link moves the view; a transcript's arrives later, and its focus waits. */
export const RegionText: readonly Run[] = [
  'TP53 sits at ',
  link('the gene', '/refs/gene'),
  ' of this region. Its second transcript, ',
  link('T2', '/refs/t2'),
  ', is read with its part of the region, so a focus on it waits until it arrives; so does one on a neighbor, ',
  link('WRAP53', '/refs/wrap53'),
  '.',
];
