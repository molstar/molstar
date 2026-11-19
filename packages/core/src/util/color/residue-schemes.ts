/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Residue and nucleotide color schemes as data: dictionaries of one-letter codes to colors (`#rrggbb`), groups with
 * labels, and the scales they come from, each with its provenance (published colors and values, no code). Sequence
 * viewers expand them into their own themes; a Mol* color theme may read them too.
 */

export interface ResidueSchemeGroup {
  readonly label: string;
  /** One-letter codes, upper case. */
  readonly letters: string;
  readonly color: string;
}

export interface ResidueScheme {
  readonly name: string;
  readonly label: string;
  readonly alphabet: 'protein' | 'nucleotide';
  /** A color by one-letter code (upper case); a letter not listed is not colored by the scheme. */
  readonly colors: Readonly<Record<string, string>>;
  /** The groups the colors come from, when the scheme is one of groups. */
  readonly groups?: readonly ResidueSchemeGroup[];
  /** Whether the colors stay apart under the simulated color-vision deficiencies. */
  readonly cvdSafe: boolean;
  readonly provenance: string;
}

function fromGroups(groups: readonly ResidueSchemeGroup[]): Record<string, string> {
  const out: Record<string, string> = {};
  for (const g of groups) for (const l of g.letters) out[l] = g.color;
  return out;
}

function scheme(s: Omit<ResidueScheme, 'colors'> & { colors?: Record<string, string> }): ResidueScheme {
  return { ...s, colors: s.colors ?? fromGroups(s.groups ?? []) };
}

/** The twenty amino acids' one-letter codes and their three-letter names. */
export const AminoAcidNames: Readonly<Record<string, string>> = {
  A: 'ALA',
  R: 'ARG',
  N: 'ASN',
  D: 'ASP',
  C: 'CYS',
  Q: 'GLN',
  E: 'GLU',
  G: 'GLY',
  H: 'HIS',
  I: 'ILE',
  L: 'LEU',
  K: 'LYS',
  M: 'MET',
  F: 'PHE',
  P: 'PRO',
  S: 'SER',
  T: 'THR',
  W: 'TRP',
  Y: 'TYR',
  V: 'VAL',
};

export const ResidueSchemes = {
  /** Clustal X's residue groups and colors; its column thresholds are rules over a column's shares (not data here). */
  'clustal-x': scheme({
    name: 'clustal-x',
    label: 'Clustal X',
    alphabet: 'protein',
    groups: [
      { label: 'hydrophobic', letters: 'AILMFWV', color: '#80a0f0' },
      { label: 'positive', letters: 'KR', color: '#f01505' },
      { label: 'negative', letters: 'ED', color: '#c048c0' },
      { label: 'polar', letters: 'NQST', color: '#15c015' },
      { label: 'cysteine', letters: 'C', color: '#f08080' },
      { label: 'glycine', letters: 'G', color: '#f09048' },
      { label: 'proline', letters: 'P', color: '#c0c000' },
      { label: 'aromatic', letters: 'HY', color: '#15a4a4' },
    ],
    cvdSafe: false,
    provenance: 'Clustal X colour scheme (Thompson et al. 1997), as tabulated in the Jalview documentation',
  }),
  zappo: scheme({
    name: 'zappo',
    label: 'Zappo',
    alphabet: 'protein',
    groups: [
      { label: 'aliphatic', letters: 'ILVAM', color: '#ffafaf' },
      { label: 'aromatic', letters: 'FWY', color: '#ffc800' },
      { label: 'positive', letters: 'KRH', color: '#6464ff' },
      { label: 'negative', letters: 'DE', color: '#ff0000' },
      { label: 'hydrophilic', letters: 'STNQ', color: '#00ff00' },
      { label: 'conformational', letters: 'PG', color: '#ff00ff' },
      { label: 'cysteine', letters: 'C', color: '#ffff00' },
    ],
    cvdSafe: false,
    provenance: 'Zappo colour scheme, physico-chemical groups, as tabulated in the Jalview documentation',
  }),
  taylor: scheme({
    name: 'taylor',
    label: 'Taylor',
    alphabet: 'protein',
    colors: {
      A: '#ccff00',
      R: '#0000ff',
      N: '#cc00ff',
      D: '#ff0000',
      C: '#ffff00',
      Q: '#ff00cc',
      E: '#ff0066',
      G: '#ff9900',
      H: '#0066ff',
      I: '#66ff00',
      L: '#33ff00',
      K: '#6600ff',
      M: '#00ff00',
      F: '#00ff66',
      P: '#ffcc00',
      S: '#ff3300',
      T: '#ff6600',
      W: '#00ccff',
      Y: '#00ffcc',
      V: '#99ff00',
    },
    cvdSafe: false,
    provenance: 'Taylor W.R. 1997, "Residual colours: a proposal for aminochromography", Protein Eng. 10:743-746',
  }),
  helix: scheme({
    name: 'helix',
    label: 'Helix propensity',
    alphabet: 'protein',
    colors: {
      A: '#e718e7',
      R: '#6f906f',
      N: '#1be41b',
      D: '#778877',
      C: '#23dc23',
      Q: '#926d92',
      E: '#ff00ff',
      G: '#00ff00',
      H: '#758a75',
      I: '#8a758a',
      L: '#ae51ae',
      K: '#a05fa0',
      M: '#ef10ef',
      F: '#986798',
      P: '#00ff00',
      S: '#36c936',
      T: '#47b847',
      W: '#8a758a',
      Y: '#21de21',
      V: '#857a85',
    },
    cvdSafe: false,
    provenance: 'Helix propensity (Chou and Fasman 1978), green to magenta, as tabulated in the Jalview documentation',
  }),
  strand: scheme({
    name: 'strand',
    label: 'Strand propensity',
    alphabet: 'protein',
    colors: {
      A: '#5858a7',
      R: '#6b6b94',
      N: '#64649b',
      D: '#2121de',
      C: '#9d9d62',
      Q: '#8c8c73',
      E: '#0000ff',
      G: '#4949b6',
      H: '#60609f',
      I: '#ecec13',
      L: '#b2b24d',
      K: '#4747b8',
      M: '#82827d',
      F: '#c2c23d',
      P: '#2323dc',
      S: '#4949b6',
      T: '#9d9d62',
      W: '#c0c03f',
      Y: '#d3d32c',
      V: '#ffff00',
    },
    cvdSafe: false,
    provenance: 'Strand propensity (Chou and Fasman 1978), blue to yellow, as tabulated in the Jalview documentation',
  }),
  turn: scheme({
    name: 'turn',
    label: 'Turn propensity',
    alphabet: 'protein',
    colors: {
      A: '#2cd3d3',
      R: '#708f8f',
      N: '#ff0000',
      D: '#e81717',
      C: '#a85757',
      Q: '#3fc0c0',
      E: '#778888',
      G: '#ff0000',
      H: '#708f8f',
      I: '#00ffff',
      L: '#1ce3e3',
      K: '#7e8181',
      M: '#1ee1e1',
      F: '#1ee1e1',
      P: '#f60909',
      S: '#e11e1e',
      T: '#738c8c',
      W: '#738c8c',
      Y: '#9d6262',
      V: '#07f8f8',
    },
    cvdSafe: false,
    provenance: 'Turn propensity (Chou and Fasman 1978), cyan to red, as tabulated in the Jalview documentation',
  }),
  buried: scheme({
    name: 'buried',
    label: 'Buried index',
    alphabet: 'protein',
    colors: {
      A: '#00a35c',
      R: '#00fc03',
      N: '#00eb14',
      D: '#00eb14',
      C: '#0000ff',
      Q: '#00f10e',
      E: '#00f10e',
      G: '#009d62',
      H: '#00d52a',
      I: '#0054ab',
      L: '#007b84',
      K: '#00ff00',
      M: '#009768',
      F: '#008778',
      P: '#00e01f',
      S: '#00d52a',
      T: '#00db24',
      W: '#00a857',
      Y: '#00e619',
      V: '#005fa0',
    },
    cvdSafe: false,
    provenance: 'Buried index (Janin 1979), blue to green, as tabulated in the Jalview documentation',
  }),
  nucleotide: scheme({
    name: 'nucleotide',
    label: 'Nucleotide',
    alphabet: 'nucleotide',
    colors: { A: '#64f73f', C: '#ffb340', G: '#eb413c', T: '#3c88ee', U: '#3c88ee' },
    cvdSafe: false,
    provenance: 'Nucleotide colours as tabulated in the Jalview documentation',
  }),
  'purine-pyrimidine': scheme({
    name: 'purine-pyrimidine',
    label: 'Purine / pyrimidine',
    alphabet: 'nucleotide',
    groups: [
      { label: 'purine', letters: 'AGR', color: '#ff83fa' },
      { label: 'pyrimidine', letters: 'CTUY', color: '#40e0d0' },
    ],
    cvdSafe: true,
    provenance: 'Purine and pyrimidine colours as tabulated in the Jalview documentation',
  }),
  'nucleotide-cb': scheme({
    name: 'nucleotide-cb',
    label: 'Nucleotide (Okabe–Ito)',
    alphabet: 'nucleotide',
    colors: { A: '#009e73', C: '#0072b2', G: '#e69f00', T: '#d55e00', U: '#d55e00' },
    cvdSafe: true,
    provenance: 'Okabe and Ito 2008, "Color Universal Design": four of its eight colors, one a base',
  }),
} as const satisfies Readonly<Record<string, ResidueScheme>>;

export type ResidueSchemeName = keyof typeof ResidueSchemes;

/**
 * "Experimentally determined hydrophobicity scale for proteins at membrane interfaces"
 * by Wimely and White (doi:10.1038/nsb1096-842)
 * http://blanco.biomol.uci.edu/Whole_residue_HFscales.txt
 * https://www.nature.com/articles/nsb1096-842
 */
export const ResidueHydrophobicity = {
  // AA  DGwif   DGwoct  Oct-IF
  ALA: [0.17, 0.5, 0.33],
  ARG: [0.81, 1.81, 1.0],
  ASN: [0.42, 0.85, 0.43],
  ASP: [1.23, 3.64, 2.41],
  ASH: [-0.07, 0.43, 0.5],
  CYS: [-0.24, -0.02, 0.22],
  GLN: [0.58, 0.77, 0.19],
  GLU: [2.02, 3.63, 1.61],
  GLH: [-0.01, 0.11, 0.12],
  GLY: [0.01, 1.15, 1.14],
  // "His+": [  0.96,  2.33,  1.37 ],
  HIS: [0.17, 0.11, -0.06],
  ILE: [-0.31, -1.12, -0.81],
  LEU: [-0.56, -1.25, -0.69],
  LYS: [0.99, 2.8, 1.81],
  MET: [-0.23, -0.67, -0.44],
  PHE: [-1.13, -1.71, -0.58],
  PRO: [0.45, 0.14, -0.31],
  SER: [0.13, 0.46, 0.33],
  THR: [0.14, 0.25, 0.11],
  TRP: [-1.85, -2.09, -0.24],
  TYR: [-0.94, -0.71, 0.23],
  VAL: [0.07, -0.46, -0.53],

  // charmm ff
  HSD: [0.17, 0.11, -0.06],
  HSE: [0.17, 0.11, -0.06],
  HSP: [0.96, 2.33, 1.37],

  // amber ff
  HID: [0.17, 0.11, -0.06],
  HIE: [0.17, 0.11, -0.06],
  HIP: [0.96, 2.33, 1.37],
};
export const DefaultResidueHydrophobicity = [0.0, 0.0, 0.0];
