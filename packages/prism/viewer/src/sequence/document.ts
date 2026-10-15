/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { PrismPlddtBands, PrismSecondaryColors } from './palette.js';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { StructureFormat } from '@molstar/seq-molstar-poc';

// The document the sequence model builds (prism-sequence.md §3): one view with one section over a `molstar:structure`
// complex whose file the host already parsed; with no blocks the section takes every polymer unit in order, then the
// ligands. The theme (§6) is written into it: Seq*'s themes draw the colors.

/** The complex's ref in the document, and the name the link maps the 3D structure onto. */
export const ComplexRef = 'structure';
/** The section's ref. */
export const SectionRef = 'main';
/** The letters track's ref, and the pLDDT track's (a predicted model). */
export const LettersRef = 'letters';
export const PlddtRef = 'plddt';
/** The tracks' refs by what they show (the track menus are by these: `trackActions`). */
export const PlddtTrack = `${PlddtRef}-track`;
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
  /** A predicted model: its pLDDT as a track under the letters, in the AlphaFold bands (`PrismPlddtBands`). */
  predicted?: boolean;
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
  // One document for both styles: Text draws the ruler nowhere, secondary structure as underlines and the pLDDT as a
  // strip in the bands (Seq*'s lanes); Tracks draws each as a track.
  section.track({ ref: 'ruler', label: '' }).ruler();
  section.track({ ref: LettersRef, label: '' }).letters(letters);
  // Secondary structure: the model's own helices, strands and turns (blocks and arrows over a coil line).
  const ss = b.data.computed({
    ref: SecondaryRef,
    provider: 'molstar:secondary_structure',
    options: { of: ComplexRef },
  });
  section.track({ ref: SecondaryTrack, label: '' }).secondary({ data: ss, color: PrismSecondaryColors });
  if (o.predicted) {
    const plddt = b.data.computed({
      ref: PlddtRef,
      provider: 'molstar:plddt',
      options: { of: ComplexRef },
    });
    // The model's confidence as Seq*'s case 8 draws it: the scores as a line in the bands. Its thin band strip under the
    // line is left out: it sat apart from the line, a gap between.
    section
      .track({ ref: PlddtTrack, label: '', height: 2 })
      .values({ data: plddt, encoding: 'line', domain: [0, 100], color: PrismPlddtBands });
  }
  return b.getDocument();
}
