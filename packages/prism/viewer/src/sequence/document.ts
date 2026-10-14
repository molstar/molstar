/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { StructureFormat } from '@molstar/seq-molstar-poc';

// The document the sequence model builds (prism-sequence.md §3): one view with one section over a `molstar:structure`
// complex whose file the host already parsed; with no blocks the section takes every polymer unit in order, then the
// ligands. The theme (§6) is written into it: Seq*'s themes draw the colors.

/** The complex's ref in the document, and the name the link maps the 3D structure onto. */
export const ComplexRef = 'structure';
/** The section's ref. */
export const SectionRef = 'main';
/** The letters track's ref. */
export const LettersRef = 'letters';

/** The host file the document names its structure by: the parsed one, which the loaders' `parsed` hook returns. */
export const hostFileOf = (key: string) => `structure:${key}`;

/** The color scheme the letters use, by its registered name (`KernelColorSchemes`); null: not colored. */
export interface SequenceTheme {
  readonly letters: string | null;
}
export const DefaultSequenceTheme: SequenceTheme = { letters: 'zappo' };

export function sequenceDocument(o: {
  title: string;
  key: string;
  format: StructureFormat;
  theme: SequenceTheme;
}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: o.title });
  b.data.complex({
    ref: ComplexRef,
    provider: 'molstar:structure',
    source: { host_file: hostFileOf(o.key), format: o.format },
  } as never);
  const section = b.view({ ref: 'sequence', label: o.title }).section({ ref: SectionRef, style: 'compact' } as never);
  const letters = o.theme.letters === null ? {} : { color: { theme: 'scheme', params: { name: o.theme.letters } } };
  section.track({ ref: LettersRef, label: '' } as never).letters(letters as never);
  return b.getDocument();
}
