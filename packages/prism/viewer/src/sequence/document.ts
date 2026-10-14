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
}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: o.title });
  b.data.complex({
    ref: ComplexRef,
    provider: 'molstar:structure',
    source: { host_file: hostFileOf(o.key), format: o.format },
  } as never);
  const section = b.view({ ref: 'sequence', label: o.title }).section({ ref: SectionRef, style: 'compact' } as never);
  const letters = o.theme.letters === null ? {} : { color: { theme: 'scheme', params: { name: o.theme.letters } } };
  // The ruler is the tracks version's axis (the text version draws no ruler).
  section.track({ ref: 'ruler', label: '' } as never).ruler();
  section.track({ ref: LettersRef, label: '' } as never).letters(letters as never);
  // Secondary structure: the model's own helices, strands and turns (blocks and arrows over a coil line).
  const ss = b.data.computed({
    ref: SecondaryRef,
    provider: 'molstar:secondary_structure',
    options: { of: ComplexRef },
  } as never);
  section
    .track({ ref: SecondaryTrack, label: '' } as never)
    .secondary({ data: ss, color: PrismSecondaryColors } as never);
  if (o.predicted) {
    const plddt = b.data.computed({
      ref: PlddtRef,
      provider: 'molstar:plddt',
      options: { of: ComplexRef },
    } as never);
    section
      .track({ ref: `${PlddtRef}-track`, label: '' } as never)
      .values({ data: plddt, color: PrismPlddtBands } as never);
  }
  return b.getDocument();
}
