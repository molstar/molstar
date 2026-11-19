/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  AminoAcidNames,
  ResidueHydrophobicity,
  type ResidueScheme,
  ResidueSchemes,
} from '@molstar/core/util/color/residue-schemes';
import { Oklab } from '@molstar/core/util/color/spaces/oklab';
import { parseColor, toRgb } from './colors.js';
import type { ThemeNode } from './forms.js';
import { continuousLut, LutSize } from './palettes.js';

// Schemes (_spec/README.md "Themes", Schemes): named, curated, data. An entry says what it colors (its subject), what it
// needs (measures), whether it is color-blind safe and its substitute, and expands to a `letters`, `field` or `rules`
// form. Built-ins are entries like any other, in the default composition of the viewer and of `createSeqContext`.

/** What a scheme colors: letters of an alphabet, numbers, categories, secondary structure, chromosome bands. */
export type SchemeSubject = 'protein' | 'nucleotide' | 'numbers' | 'categories' | 'secondary' | 'bands';

export interface SchemeEntry {
  readonly name: string;
  readonly label: string;
  readonly subject: SchemeSubject;
  /** The measures it reads (a column's share, a reference row): offered only where each is registered. */
  readonly needs?: readonly string[];
  readonly cvdSafe: boolean;
  /** The scheme a reader's color-blind substitution draws instead (`theme.cvdSafe`). */
  readonly substitute?: string;
  readonly provenance: string;
  /**
   * The form it expands to, given the form's params (overrides, ESPript's `gaps`, `similarity`, `mode`). Literal colors
   * are the light theme's; `dark` derives (or gives) the dark theme's.
   */
  expand(params: Readonly<Record<string, unknown>>): ThemeNode;
  /** Literal colors in a dark table: `derive` (each color's lightness mapped into the dark band), `same`. */
  readonly dark: 'derive' | 'same';
}

/** The dark band a derived variant's lightness is mapped into (from the light band, Oklab L). */
export const DarkBand = { from: [0.4, 1] as const, to: [0.3, 0.62] as const };

/** A light theme's color (`#rrggbb`) in the dark band: hue kept, chroma clamped to the gamut. */
export function darkVariant(color: string): string {
  const c = parseColor(color);
  if (c === undefined) return color;
  const d = Oklab.mapLightness(toRgb(c), DarkBand.from, DarkBand.to);
  return `#${d.toString(16).padStart(6, '0')}`;
}

/** The named schemes: a registry a composition extends (`register`); the built-ins are entries of it. */
export class SchemeRegistry {
  private readonly entries = new Map<string, SchemeEntry>();
  constructor(entries: readonly SchemeEntry[] = []) {
    for (const e of entries) this.register(e);
  }
  register(entry: SchemeEntry): this {
    this.entries.set(entry.name, entry);
    return this;
  }
  get(name: string): SchemeEntry | undefined {
    return this.entries.get(name);
  }
  list(): SchemeEntry[] {
    return [...this.entries.values()];
  }
  /** The schemes that fit a subject, given the measures registered (what a picker offers). */
  fitting(subject: SchemeSubject, measures: ReadonlySet<string> = new Set()): SchemeEntry[] {
    return this.list().filter((e) => e.subject === subject && (e.needs ?? []).every((m) => measures.has(m)));
  }
}

const letters =
  (colors: Readonly<Record<string, string>>, groups?: readonly unknown[]) =>
  (params: Readonly<Record<string, unknown>>): ThemeNode => ({
    kind: 'letters',
    params: {
      colors: { ...colors, ...((params.overrides as object | undefined) ?? {}) },
      ...(groups ? { groups } : {}),
      default: 'token:missing',
    },
  });

const group = (letters: string, value: string, label: string): ThemeNode => ({
  kind: 'group',
  params: { letters, value, label },
});

/** A field form over a palette, with overrides merged into a categorical dictionary. */
const field =
  (name: string, palette: ThemeNode) =>
  (params: Readonly<Record<string, unknown>>): ThemeNode => {
    const overrides = params.overrides as Record<string, string> | undefined;
    if (overrides && palette.kind === 'categorical') {
      const colors = palette.params?.colors;
      const dict = typeof colors === 'object' && colors !== null && !Array.isArray(colors) ? colors : {};
      palette = { kind: 'categorical', params: { ...palette.params, colors: { ...dict, ...overrides } } };
    }
    return { kind: 'field', params: { name: typeof params.field === 'string' ? params.field : name, palette } };
  };

const continuous = (colors: unknown, extra: Record<string, unknown> = {}): ThemeNode => ({
  kind: 'continuous',
  params: { colors, ...extra },
});

/** The hydrophobicity scheme: Mol*'s one table (Wimley and White, its first scale), red hydrophobic to blue. */
function hydrophobicityColors(): Record<string, string> {
  const scale = AminoAcidNames;
  const values = Object.entries(scale).map(
    ([l, n]) => [l, (ResidueHydrophobicity as Record<string, number[]>)[n][0]] as const,
  );
  const lo = Math.min(...values.map((v) => v[1]));
  const hi = Math.max(...values.map((v) => v[1]));
  // low ΔG (hydrophobic) red, high blue: Mol*'s red-white-blue, interpolated in Oklab
  const lut = continuousLut([0xd73027ff, 0xf7f7f7ff, 0x4575b4ff], [0, 0.5, 1]);
  const out: Record<string, string> = {};
  for (const [l, v] of values) {
    const c = lut[Math.round(((v - lo) / (hi - lo)) * (LutSize - 1))];
    out[l] = `#${(c >>> 8).toString(16).padStart(6, '0')}`;
  }
  return out;
}

const residue = (name: keyof typeof ResidueSchemes, extra: Partial<SchemeEntry> = {}): SchemeEntry => {
  const s: ResidueScheme = ResidueSchemes[name];
  return {
    name,
    label: s.label,
    subject: s.alphabet,
    cvdSafe: s.cvdSafe,
    ...(s.alphabet === 'protein' && !s.cvdSafe ? { substitute: 'classes' } : {}),
    ...(s.alphabet === 'nucleotide' && !s.cvdSafe ? { substitute: 'nucleotide-cb' } : {}),
    provenance: s.provenance,
    // a safe scheme's colors were chosen apart: a derived variant would bring them together, so it keeps them
    dark: s.cvdSafe ? 'same' : 'derive',
    expand: letters(
      s.colors,
      s.groups?.map((g) => group(g.letters, g.color, g.label)),
    ),
    ...extra,
  };
};

const list = (
  name: string,
  label: string,
  subject: SchemeSubject,
  cvdSafe: boolean,
  provenance: string,
): SchemeEntry => ({
  name,
  label,
  subject,
  cvdSafe,
  ...(cvdSafe ? {} : { substitute: subject === 'categories' ? 'okabe-ito' : 'viridis' }),
  provenance,
  dark: 'same',
  expand: (params) =>
    subject === 'categories'
      ? {
          kind: 'field',
          params: {
            name: (params.field as string) ?? 'category',
            palette: { kind: 'categorical', params: { colors: name } },
          },
        }
      : { kind: 'field', params: { name: (params.field as string) ?? 'value', palette: continuous(name) } },
});

const ColumnMeasureNote = 'needs the rows layer’s column measures';

/** The built-in schemes (provenance per entry: published colors and rules, never GPL code). */
export const BuiltInSchemes: readonly SchemeEntry[] = [
  {
    name: 'classes',
    label: 'Residue classes',
    subject: 'protein',
    cvdSafe: false,
    provenance: 'Seq*’s residue class tokens (`seq-aa-*`)',
    dark: 'same',
    expand: letters({}, [
      group('AVLIMC', 'token:aa-hydrophobic', 'hydrophobic'),
      group('FWY', 'token:aa-aromatic', 'aromatic'),
      group('STNQ', 'token:aa-polar', 'polar'),
      group('KRH', 'token:aa-positive', 'positive'),
      group('DE', 'token:aa-negative', 'negative'),
      group('GP', 'token:aa-special', 'special'),
    ]),
  },
  residue('clustal-x', {
    needs: ['share'],
    provenance: `${ResidueSchemes['clustal-x'].provenance}; ${ColumnMeasureNote}`,
  }),
  residue('zappo'),
  residue('taylor'),
  {
    name: 'hydrophobicity',
    label: 'Hydrophobicity',
    subject: 'protein',
    cvdSafe: false,
    substitute: 'classes',
    provenance: 'Wimley and White 1996 (doi:10.1038/nsb1096-842), Mol*’s `ResidueHydrophobicity`, red to blue',
    dark: 'derive',
    expand: (p) => letters(hydrophobicityColors())(p),
  },
  residue('helix'),
  residue('strand'),
  residue('turn'),
  residue('buried'),
  {
    name: 'percent-identity',
    label: 'Percent identity',
    subject: 'protein',
    needs: ['identity'],
    cvdSafe: true,
    provenance: `Jalview’s percentage identity, blue shades by the column’s identity; ${ColumnMeasureNote}`,
    dark: 'same',
    expand: () => ({
      kind: 'rules',
      params: { else: 'token:track' },
      args: [
        { kind: 'rule', params: { when: [{ kind: 'gt', args: [{ kind: 'identity' }, 0.8] }], value: '#6464ff' } },
        { kind: 'rule', params: { when: [{ kind: 'gt', args: [{ kind: 'identity' }, 0.6] }], value: '#9999ff' } },
        { kind: 'rule', params: { when: [{ kind: 'gt', args: [{ kind: 'identity' }, 0.4] }], value: '#ccccff' } },
      ],
    }),
  },
  {
    name: 'blosum62',
    label: 'BLOSUM62',
    subject: 'protein',
    needs: ['matches-reference'],
    cvdSafe: true,
    provenance: `Jalview’s BLOSUM62 score against the column’s consensus; ${ColumnMeasureNote}`,
    dark: 'same',
    expand: () => ({
      kind: 'rules',
      params: { else: 'token:track' },
      args: [
        {
          kind: 'rule',
          params: { when: [{ kind: 'eq', args: [{ kind: 'matches-reference' }, 1] }], value: '#9a9aff' },
        },
        {
          kind: 'rule',
          params: { when: [{ kind: 'gt', args: [{ kind: 'matches-reference' }, 0] }], value: '#cccfff' },
        },
      ],
    }),
  },
  {
    name: 'espript',
    label: 'ESPript',
    subject: 'protein',
    needs: ['share'],
    cvdSafe: true,
    provenance: `ESPript (Robert and Gouet 2014): identical columns filled, similar residues by groups (an approximation of its score); ${ColumnMeasureNote}`,
    dark: 'same',
    expand: (p) => {
      const bw = p.mode === 'bw';
      const similarity = typeof p.similarity === 'number' ? p.similarity : 0.7;
      return {
        kind: 'rules',
        params: { else: 'token:track' },
        args: [
          {
            kind: 'rule',
            params: {
              when: [{ kind: 'ge', args: [{ kind: 'share', params: { gaps: p.gaps ?? 'count' } }, 1] }],
              value: bw ? '#000000' : '#d62728',
            },
          },
          {
            kind: 'rule',
            params: {
              when: [{ kind: 'ge', args: [{ kind: 'share', params: { by: 'group' } }, similarity] }],
              value: bw ? '#c8c8c8' : '#fbe3e3',
            },
          },
        ],
      };
    },
  },
  {
    name: 'nucleotide-classes',
    label: 'Nucleotides',
    subject: 'nucleotide',
    cvdSafe: false,
    substitute: 'nucleotide-cb',
    provenance: 'Seq*’s nucleotide tokens (`seq-nt-*`)',
    dark: 'same',
    expand: (p) => ({
      kind: 'letters',
      params: {
        colors: {
          A: 'token:nt-a',
          C: 'token:nt-c',
          G: 'token:nt-g',
          T: 'token:nt-t',
          U: 'token:nt-t',
          '-': 'token:missing',
          '.': 'token:missing',
          ...((p.overrides as object | undefined) ?? {}),
        },
        default: 'token:nt-other',
      },
    }),
  },
  residue('nucleotide'),
  residue('purine-pyrimidine'),
  residue('nucleotide-cb'),
  {
    name: 'ramp',
    label: 'Ramp',
    subject: 'numbers',
    cvdSafe: true,
    provenance: 'Seq*’s ramp tokens (`seq-ramp-*`), monotone in lightness',
    dark: 'same',
    expand: (p) => ({ kind: 'field', params: { name: (p.field as string) ?? 'value', palette: continuous('ramp') } }),
  },
  {
    name: 'diverging',
    label: 'Diverging',
    subject: 'numbers',
    cvdSafe: true,
    provenance: 'Seq*’s diverging tokens (`seq-div-*`), centered on 0 unless the form says `mid`',
    dark: 'same',
    expand: (p) => ({
      kind: 'field',
      params: {
        name: (p.field as string) ?? 'value',
        palette: continuous('diverging', { mid: (p.mid as number) ?? 0 }),
      },
    }),
  },
  list('viridis', 'Viridis', 'numbers', true, 'Matplotlib’s viridis (van der Walt and Smith), Mol*’s `ColorLists`'),
  list('cividis', 'Cividis', 'numbers', true, 'Nuñez, Anderton and Renslow 2018, Mol*’s `ColorLists`'),
  list('magma', 'Magma', 'numbers', true, 'Matplotlib’s magma, Mol*’s `ColorLists`'),
  list('blues', 'Blues', 'numbers', true, 'ColorBrewer Blues (Brewer), Mol*’s `ColorLists`'),
  list('rd-bu', 'Red–blue', 'numbers', true, 'ColorBrewer RdBu (Brewer), Mol*’s `ColorLists`'),
  {
    name: 'plddt',
    label: 'pLDDT',
    subject: 'numbers',
    cvdSafe: true,
    provenance: 'AlphaFold DB’s model confidence bands (Jumper et al. 2021): <50, 50–70, 70–90, ≥90',
    dark: 'same',
    expand: (p) => ({
      kind: 'field',
      params: {
        name: (p.field as string) ?? 'value',
        palette: {
          kind: 'binned',
          params: { edges: [0, 50, 70, 90, 100], colors: ['#ff7d45', '#ffdb13', '#65cbf3', '#0053d6'] },
        },
      },
    }),
  },
  {
    name: 'categories',
    label: 'Categories',
    subject: 'categories',
    cvdSafe: true,
    provenance: 'Seq*’s category tokens (`seq-cat-*`), checked for color-vision deficiencies',
    dark: 'same',
    expand: (p) => ({
      kind: 'field',
      params: {
        name: (p.field as string) ?? 'category',
        palette: { kind: 'categorical', params: { colors: 'categories' } },
      },
    }),
  },
  list('okabe-ito', 'Okabe–Ito', 'categories', true, 'Okabe and Ito 2008, Mol*’s `ColorLists`'),
  list('tableau10', 'Tableau 10', 'categories', false, 'Tableau 10, Mol*’s `ColorLists`'),
  list('category10', 'Category 10', 'categories', false, 'D3’s category10, Mol*’s `ColorLists`'),
  {
    name: 'secondary',
    label: 'Secondary structure',
    subject: 'secondary',
    cvdSafe: true,
    provenance: 'Seq*’s secondary-structure tokens (`seq-ss-*`)',
    dark: 'same',
    expand: field('kind', {
      kind: 'categorical',
      params: {
        colors: { helix: 'token:ss-helix', strand: 'token:ss-strand', turn: 'token:ss-turn' },
        unlisted: 'missing',
      },
    }),
  },
  {
    name: 'giemsa',
    label: 'Giemsa stain',
    subject: 'bands',
    cvdSafe: true,
    provenance:
      'UCSC’s cytoband stains as greys (darker for stronger stain, flipped in the dark theme); `acen` in `seq-centromere`',
    dark: 'same',
    expand: field('stain', {
      kind: 'categorical',
      params: {
        colors: {
          gneg: 'light-dark(#ffffff,#1a1d23)',
          gpos25: 'light-dark(#c8c8c8,#4a4f58)',
          gpos50: 'light-dark(#969696,#7a808a)',
          gpos75: 'light-dark(#646464,#a8aeb8)',
          gpos100: 'light-dark(#2d271e,#cfd6e1)',
          gvar: 'light-dark(#dcdcdc,#30343c)',
          stalk: 'light-dark(#787878,#8d95a1)',
          acen: 'token:centromere',
        },
        unlisted: 'missing',
      },
    }),
  },
];

/** A registry of the built-in schemes. */
export function builtInSchemes(): SchemeRegistry {
  return new SchemeRegistry(BuiltInSchemes);
}

/** The default scheme of a letters layer by its alphabet. */
export function lettersDefault(alphabet: 'protein' | 'dna' | 'rna'): string {
  return alphabet === 'protein' ? 'classes' : 'nucleotide-classes';
}

/** Every literal color a scheme's form holds (the conformance suite reads them). */
export function schemeColors(entry: SchemeEntry): string[] {
  const out: string[] = [];
  const walk = (v: unknown) => {
    if (typeof v === 'string') {
      if (v.startsWith('#') || v.startsWith('token:') || v.startsWith('light-dark(')) out.push(v);
    } else if (Array.isArray(v)) v.forEach(walk);
    else if (v && typeof v === 'object') Object.values(v).forEach(walk);
  };
  walk(entry.expand({}));
  return out;
}
