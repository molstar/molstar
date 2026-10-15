/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ColorSchemeEntry, SvsRegistryEntry, ThemeRef } from '@molstar/seq-core/format';

// The sequence model's own residue palette (prism-sequence.md §6): the residue classes of Zappo (aliphatic, aromatic,
// positive, negative, polar, the two special backbones, cysteine), in colder, softer hues, one set for a light ink and
// one for a dark. A Seq* color scheme registered from the viewer's side (`colorSchemes` of a registry entry): the
// document names it, the themes and the legend read it; Seq* knows nothing of Prism.

const groups = (
  colors: Record<'aliphatic' | 'aromatic' | 'positive' | 'negative' | 'polar' | 'special' | 'cysteine', string>,
) => [
  { letters: 'AILMV', color: colors.aliphatic, label: 'Aliphatic' },
  { letters: 'FWY', color: colors.aromatic, label: 'Aromatic' },
  { letters: 'KRH', color: colors.positive, label: 'Positive' },
  { letters: 'DE', color: colors.negative, label: 'Negative' },
  { letters: 'STNQ', color: colors.polar, label: 'Polar' },
  { letters: 'PG', color: colors.special, label: 'Proline, glycine' },
  { letters: 'C', color: colors.cysteine, label: 'Cysteine' },
];

export const PrismResidues: ColorSchemeEntry = {
  name: 'prism:residues',
  label: 'Prism residues',
  description: 'residues by class (aliphatic, aromatic, charged, polar, special) in Prism’s cool, soft hues',
  provenance: 'Zappo’s classes (Jalview’s published grouping); the hues are Prism’s own',
  colorblind: false,
  scheme: {
    kind: 'letters',
    // Seq* judges a scheme's fit by the section's first block, and a structure's first chains may be nucleic acids: all
    // the polymer alphabets fit; the letters are colored by class where they are the same letters (A, C, G, T).
    alphabets: ['protein', 'dna', 'rna'],
    // Light ink: mid tones, dark enough to read on a pale canvas.
    ...{
      groups: groups({
        aliphatic: '#6f879f',
        aromatic: '#7a6cc2',
        positive: '#3d7cc0',
        negative: '#c25d7b',
        polar: '#2f9a8d',
        special: '#b98a45',
        cysteine: '#a99a2e',
      }),
    },
    dark: {
      groups: groups({
        aliphatic: '#9db0c6',
        aromatic: '#a99be4',
        positive: '#6aa9ea',
        negative: '#ea8fa8',
        polar: '#62c5b6',
        special: '#e0b56d',
        cysteine: '#d9c866',
      }),
    },
  },
};

/** The registry entry the model adds to Seq*'s: its schemes. */
export const PrismRegistry: SvsRegistryEntry = { colorSchemes: [PrismResidues] };

/** The schemes the letters' menu offers, the default first (a name per `colorScheme`; null: not colored). */
export const LetterSchemes: readonly (string | null)[] = ['prism:residues', 'zappo', 'taylor', 'hydrophobicity', null];

// The pLDDT track's bands: AlphaFold DB's four (very low < 50 ≤ low < 70 ≤ confident < 90 ≤ very high), the official
// hues softened to the residue palette's key, still recognizably AlphaFold's blue and orange. Cells, so one set for a
// light and a dark ink. The model knows them (it set the theme), so its legend lists them (prism-sequence.md §6).

const BANDS = [
  { label: 'Very high (pLDDT > 90)', from: 90, color: '#3f6fc4' },
  { label: 'Confident (70–90)', from: 70, color: '#7ec1e3' },
  { label: 'Low (50–70)', from: 50, color: '#e8cc5a' },
  { label: 'Very low (< 50)', from: 0, color: '#e89468' },
] as const;

/** The track's color theme: a discrete palette over the value, the bands' lower bounds absolute. */
export const PrismPlddtBands: ThemeRef = {
  theme: 'field_palette',
  params: {
    field: 'value',
    palette: {
      kind: 'discrete',
      mode: 'absolute',
      colors: [...BANDS].reverse().map((b) => [b.color, b.from]),
    },
  },
};

/** The bands as the legend lists them, high to low. */
export const PlddtLegendItems: readonly { label: string; color: string }[] = BANDS.map((b) => ({
  label: b.label,
  color: b.color,
}));

// The secondary structure track's classes, softened to the residue palette's key: helices rose, strands amber, turns
// teal, the coil the line between them. One set for a light and a dark ink (they are blocks and arrows).
const SECONDARY = [
  { kind: 'helix', label: 'Helix', color: '#d1708a' },
  { kind: 'strand', label: 'Strand', color: '#e0a73f' },
  { kind: 'turn', label: 'Turn', color: '#4fa3a5' },
] as const;

/** The track's color theme: by the range's class. */
export const PrismSecondaryColors: ThemeRef = {
  theme: 'field_palette',
  params: {
    field: 'category',
    palette: {
      kind: 'categorical',
      colors: Object.fromEntries(SECONDARY.map((s) => [s.kind, s.color])),
    },
  },
};

export const SecondaryLegendItems: readonly { label: string; color: string }[] = SECONDARY.map((s) => ({
  label: s.label,
  color: s.color,
}));
