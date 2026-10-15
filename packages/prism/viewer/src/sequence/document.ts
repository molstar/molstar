/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { PrismSecondaryColors } from './palette.js';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { StructureFormat } from '@molstar/seq-molstar-poc';

// The document the sequence model builds (prism-sequence.md §3): one view with one section over a `molstar:structure`
// complex whose file the host already parsed; with no blocks the section takes every polymer unit in order, then the
// ligands. The theme (§6) is written into it: Seq*'s themes draw the colors.

/** The complex's ref in the document, and the name the link maps the 3D structure onto. */
export const ComplexRef = 'structure';
/** The section's ref. */
export const SectionRef = 'main';
/** The data node's ref: behaviors add their data under it (sequence.md §3). */
export const DataRef = 'data';
/** The letters track's ref. */
export const LettersRef = 'letters';
/** The tracks' refs by what they show (the track menus are by these: `trackActions`). */
export const SecondaryRef = 'ss';
export const SecondaryTrack = `${SecondaryRef}-track`;

/** The host file the document names its structure by: the parsed one, which the loaders' `parsed` hook returns. */
export const hostFileOf = (key: string) => `structure:${key}`;

/** The color scheme the letters use, by its registered name (`KernelColorSchemes`); null: not colored. */
export interface SequenceTheme {
  readonly letters: string | null;
}
export const DefaultSequenceTheme: SequenceTheme = { letters: 'prism:residues' };

export function sequenceDocument(o: {
  title: string;
  key: string;
  format: StructureFormat;
  theme: SequenceTheme;
  /** Only the ligands of these chemical components become units (the ones the 3D view draws); absent: all. */
  ligands?: readonly string[];
}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: o.title });
  b.data.complex({
    ref: ComplexRef,
    provider: 'molstar:structure',
    source: { host_file: hostFileOf(o.key), format: o.format },
    ...(o.ligands ? { options: { ligand_comps: [...o.ligands] } } : {}),
  });
  const section = b.view({ ref: 'sequence', label: o.title }).section({ ref: SectionRef, style: 'compact' });
  const letters = o.theme.letters === null ? {} : { color: { theme: 'scheme', params: { name: o.theme.letters } } };
  // One document for both styles: Text draws the ruler nowhere, secondary structure as underlines and a behavior's values
  // (a pLDDT track) as strips (Seq*'s lanes); Tracks draws each as a track.
  section.track({ ref: 'ruler', label: '' }).ruler();
  section.track({ ref: LettersRef, label: '' }).letters(letters);
  // Secondary structure: the model's own helices, strands and turns (blocks and arrows over a coil line).
  const ss = b.data.computed({
    ref: SecondaryRef,
    provider: 'molstar:secondary_structure',
    options: { of: ComplexRef },
  });
  section.track({ ref: SecondaryTrack, label: '' }).secondary({ data: ss, color: PrismSecondaryColors });
  // Behaviors add under the data node by its ref (the builder gives it none).
  const doc = b.getDocument();
  const data = doc.root.children?.find((n) => n.kind === 'data');
  if (data) data.ref = DataRef;
  return doc;
}
