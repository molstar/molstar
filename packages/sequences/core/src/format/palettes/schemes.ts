/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Color schemes (seqspec:format#color-schemes): named schemes as data, the `colorSchemes` registry kind's kernel
// entries, and the one expansion every viewer uses (`schemeTheme`): a scheme becomes the `letter`, `field_palette` or
// `rules` theme it stands for, its dark variant when the skin is dark. What a scheme fits is data on it (alphabets, a
// subject); the core knows no biology beyond these tables.

import { type Alphabet, type Annotation, Issue, Molecules, type SeqSet } from '@molstar/seq-core/model';
import type { ThemeRef } from '../doc.js';
import type { ColorSchemeEntry, LayerTypeEntry, SvsRegistries } from '../registry.js';
import type { ResolvedSection } from '../resolve/section.js';
import { type CheckContext, isPlainObject, ptr, report, show, typeError } from '../schema/field-schema.js';
import { Color } from './color.js';
import { AminoGroupColors } from './dictionaries.js';
import { LetterDictionaryNames, letterColors } from './letters.js';
import { ColorLists } from './lists.js';

export type ColorStop = string | readonly [string, number];
export interface LetterGroup {
  readonly letters: string;
  readonly color: string;
  readonly label?: string;
}
export interface LetterColors {
  /** A letter dictionary's name, or `{letter: Color}`. */
  readonly colors?: string | Readonly<Record<string, string>>;
  /** After `colors`; a group's `label` names it in the legend. */
  readonly groups?: readonly LetterGroup[];
  /** The alphabets' other letters; null or absent: uncolored. */
  readonly default?: string | null;
}
export interface LetterScheme extends LetterColors {
  readonly kind: 'letters';
  readonly alphabets: readonly Alphabet[];
  readonly dark?: 'same' | LetterColors;
}
export interface RampScheme {
  readonly kind: 'sequential' | 'diverging';
  readonly colors: readonly ColorStop[];
  readonly dark?: 'same' | readonly ColorStop[];
}
export interface CategoricalScheme {
  readonly kind: 'categorical';
  readonly colors: readonly string[] | Readonly<Record<string, string>>;
  readonly missing?: string | null;
  readonly dark?: 'same' | readonly string[] | Readonly<Record<string, string>>;
}
/** Registered only: a theme (a `rules` table) over letters; `columns`: only where an alignment's columns are held. */
export interface ThemeScheme {
  readonly kind: 'theme';
  readonly alphabets: readonly Alphabet[];
  readonly columns: boolean;
  readonly theme: ThemeRef;
  readonly dark?: ThemeRef;
  /** What the legend lists. */
  readonly legend?: readonly { readonly label: string; readonly color: string }[];
}
export type SchemeDef = LetterScheme | RampScheme | CategoricalScheme | ThemeScheme;

/** What a layer's color is about: letters of an alphabet, numbers, or categories. */
export type SchemeSubject =
  | { readonly kind: 'letters'; readonly alphabet: Alphabet | null; readonly columns: boolean }
  | { readonly kind: 'values' }
  | { readonly kind: 'categories' };

/** The `scheme` theme's params (seqspec:interaction#color-schemes-at-runtime). */
export interface SchemeParams {
  readonly name?: string | null;
  readonly custom?: SchemeDef | null;
  readonly overrides?: Readonly<Record<string, string>>;
  readonly field?: string | null;
  readonly domain?: readonly [number | null, number | null] | null;
  readonly reverse?: boolean;
}

// --- the kernel's schemes -----------------------------------------------------------------------------------------

/** Jalview's percent identity (seqspec:format#themes): the column's consensus letter where it holds more than t. */
export function percentIdentityRules(colors: readonly string[] = ['#6464ff', '#9999ff', '#ccccff']): unknown[] {
  return [0.8, 0.6, 0.4].map((t, k) => ({
    when: [
      { measure: { rank: 'letter' }, max: 0 },
      { measure: { share: 'letter' }, min: t + 1e-9 },
    ],
    color: colors[k],
    label: `> ${t * 100}%`,
  }));
}
/** Letter groups written `letters color label; ...` (`AILMV #ffafaf aliphatic; FWY #ffc800 aromatic`). */
const groups = (alphabets: readonly Alphabet[], src: string): LetterScheme => ({
  kind: 'letters',
  alphabets,
  groups: src.split('; ').map((x) => {
    const [letters, color, ...label] = x.split(' ');
    return { letters, color, label: label.join(' ') };
  }),
});
const P: Alphabet[] = ['protein'];
const NT: Alphabet[] = ['dna', 'rna'];
const list = (kind: 'sequential' | 'diverging' | 'categorical', name: keyof typeof ColorLists) =>
  ({ kind, colors: ColorLists[name].map(Color.toHex) }) as SchemeDef;
const SCHEMES: Readonly<Record<string, SchemeDef>> = {
  clustal_x: {
    kind: 'theme',
    alphabets: P,
    columns: false,
    theme: { theme: 'clustal_x' },
    legend: Object.entries(AminoGroupColors).map(([label, c]) => ({ label, color: Color.toHex(c) })),
  },
  zappo: groups(
    P,
    'AILMV #ffafaf aliphatic; FWY #ffc800 aromatic; KRH #6464ff positive; DE #ff0000 negative; STNQ #00ff00 hydrophilic; PG #ff00ff conformational; C #ffff00 cysteine',
  ),
  taylor: { kind: 'letters', alphabets: P, colors: 'Taylor' },
  hydrophobicity: { kind: 'letters', alphabets: P, colors: 'Hydrophobicity' },
  percent_identity: {
    kind: 'theme',
    alphabets: ['protein', 'dna', 'rna', 'generic'],
    columns: true,
    theme: { theme: 'percent_identity' },
    dark: { theme: 'rules', params: { rules: percentIdentityRules(['#8080ff', '#5555b5', '#383870']) } },
  },
  nucleotide: groups(NT, 'A #64f73f A; C #ffb340 C; G #eb413c G; TU #3c88ee T, U'),
  purine_pyrimidine: groups(NT, 'AGR #ff83fa purine; CTUY #40e0d0 pyrimidine'),
  nucleotide_cb: groups(NT, 'A #009e73 A; C #0072b2 C; G #f0e442 G; TU #d55e00 T, U'),
  viridis: list('sequential', 'Viridis'),
  cividis: list('sequential', 'Cividis'),
  magma: list('sequential', 'Magma'),
  blues: list('sequential', 'Blues'),
  rd_bu: list('diverging', 'RdBu'),
  category10: list('categorical', 'Category10'),
  okabe_ito: {
    kind: 'categorical',
    colors: ['#e69f00', '#56b4e9', '#009e73', '#f0e442', '#0072b2', '#d55e00', '#cc79a7'],
  },
  tableau10: list('categorical', 'Tableau10'),
};
/** The v1 set, a row each: name | label | colorblind-safe | what it shows | provenance. */
const TABLE = `
clustal_x | Clustal X | no | residues by group where the column is conserved | Jalview's published rules (jalview.org help, Clustal X), written from them: no GPL code
zappo | Zappo | no | residues by physicochemical property | Jalview's published colors (jalview.org help, Colour Schemes), via Nightingale (MIT)
taylor | Taylor | no | each residue its own hue, similar residues near | Taylor (1997) Protein Eng. 10:743, as Jalview publishes it, via Nightingale (MIT)
hydrophobicity | Hydrophobicity | no | Kyte–Doolittle: red most hydrophobic, blue most hydrophilic | Kyte & Doolittle (1982) J. Mol. Biol. 157:105; Jalview's colors via Nightingale (MIT)
percent_identity | Percent identity | yes | the consensus residue where more than 80, 60, 40% of the column agree | Jalview's published thresholds and colors (jalview.org help, Percentage Identity); the dark variant ours
nucleotide | Nucleotide | no | A, C, G and T or U each a color | Jalview's published colors (jalview.org help, Colour Schemes), via Nightingale (MIT)
purine_pyrimidine | Purine/pyrimidine | no | purines and pyrimidines | Jalview's published colors (jalview.org help, Colour Schemes), via Nightingale (MIT)
nucleotide_cb | Nucleotide, Okabe–Ito | yes | A, C, G and T or U from the Okabe–Ito colors | Okabe & Ito (2008), Color Universal Design; the assignment ours
viridis | Viridis | yes | perceptually uniform, dark low to yellow high | van der Walt & Smith, Matplotlib (CC0), as MVS's list (Mol*, MIT)
cividis | Cividis | yes | perceptually uniform for color vision deficiency, blue to yellow | Nuñez, Anderton & Renslow (2018) PLoS ONE 13:e0199239, as MVS's list (Mol*, MIT)
magma | Magma | yes | perceptually uniform, black to cream | van der Walt & Smith, Matplotlib (CC0), as MVS's list (Mol*, MIT)
blues | Blues | yes | light low to dark blue high; the format's sequential default | ColorBrewer (Brewer, Apache-2.0), as MVS's list (Mol*, MIT)
rd_bu | Red–blue | yes | diverging at the middle: red low, blue high; the format's diverging default | ColorBrewer (Brewer, Apache-2.0), as MVS's list (Mol*, MIT)
category10 | Category 10 | no | ten categories; the format's default | D3 (ISC), as MVS's list (Mol*, MIT)
okabe_ito | Okabe–Ito | yes | seven categories distinguishable with color vision deficiency | Okabe & Ito (2008), Color Universal Design; black left out (the ink)
tableau10 | Tableau 10 | no | ten categories | Tableau 10 (Maureen Stone), as MVS's list (Mol*, MIT)`;
export const KernelColorSchemes: readonly ColorSchemeEntry[] = TABLE.trim()
  .split('\n')
  .map((line) => {
    const [name, label, cb, description, provenance] = line.split(' | ');
    return { name, label, description, provenance, colorblind: cb === 'yes', scheme: SCHEMES[name] };
  });

// --- what a scheme fits -------------------------------------------------------------------------------------------

/**
 * A layer's color subject: letters when its type's color default is `letter` (with the blocks' alphabet, or `of`'s
 * polymer's, and whether the blocks are an alignment with held columns), values when all its data are numbers,
 * categories when its color default is `field_palette`; null: no scheme fits.
 */
export function schemeSubject(
  set: SeqSet,
  section: ResolvedSection,
  type: LayerTypeEntry | undefined,
  data: readonly Annotation[],
  of?: Alphabet | null,
): SchemeSubject | null {
  const d = type?.channels.find((c) => c.kind === 'color')?.default ?? null;
  const b0 = section.blocks[0];
  if (d?.theme === 'letter') {
    const aln = b0 && b0.alignment >= 0 ? set.alignments[b0.alignment] : undefined;
    const own = b0 ? (aln?.alphabet ?? set.molecules[set.units.molecule[b0.unit]]?.alphabet ?? null) : null;
    return { kind: 'letters', alphabet: of !== undefined ? of : own, columns: !!aln && !aln.lazy && of === undefined };
  }
  if (data.length > 0 && data.every((a) => a.geometry === 'values' || a.geometry === 'matrix'))
    return { kind: 'values' };
  return d?.theme === 'field_palette' ? { kind: 'categories' } : null;
}

/** Whether `def` colors `subject`. */
export function schemeFits(def: SchemeDef, s: SchemeSubject): boolean {
  if (def.kind === 'letters' || def.kind === 'theme')
    return (
      s.kind === 'letters' &&
      s.alphabet !== null &&
      def.alphabets.includes(s.alphabet) &&
      (def.kind === 'letters' || !def.columns || s.columns)
    );
  return def.kind === 'categorical' ? s.kind === 'categories' : s.kind === 'values';
}

/** The scheme a `scheme` theme's params name or hold; undefined: none (unknown, or not given). */
export function schemeOf(reg: SvsRegistries | undefined, p: SchemeParams): SchemeDef | undefined {
  if (p.custom) return p.custom;
  return typeof p.name === 'string' ? reg?.get('colorSchemes', p.name)?.scheme : undefined;
}

/** The registered schemes that fit `subject`, in registration order (the kernel's first). */
export function schemesFor(reg: SvsRegistries, subject: SchemeSubject): ColorSchemeEntry[] {
  return reg.list('colorSchemes').filter((e) => schemeFits(e.scheme, subject));
}

// --- light and dark -----------------------------------------------------------------------------------------------

const stopColor = (s: ColorStop) => Color.decodePrimary(typeof s === 'string' ? s : s[0]);
const lum = (c: number) => 0.2126 * ((c >> 16) & 255) + 0.7152 * ((c >> 8) & 255) + 0.0722 * (c & 255);

/**
 * The ramp drawn on a dark skin when none is given: a sequential ramp lighter at its first stop is reversed, a
 * diverging one has each half reversed (its center dark, its ends bright); plain lists only.
 */
export function darkRamp(kind: 'sequential' | 'diverging', colors: readonly ColorStop[]): readonly ColorStop[] {
  if (colors.length < 2 || colors.some((c) => typeof c !== 'string')) return colors;
  if (kind === 'sequential')
    return lum(stopColor(colors[0])) > lum(stopColor(colors.at(-1) as ColorStop)) ? [...colors].reverse() : colors;
  const h = colors.length >> 1;
  const mid = colors.length % 2 === 1 ? [colors[h]] : [];
  return [...colors.slice(0, h).reverse(), ...mid, ...colors.slice(colors.length - h).reverse()];
}

/** The scheme's variant for the skin, `dark` resolved: given, "same", or derived. */
export function schemeVariant(def: SchemeDef, dark: boolean): SchemeDef {
  if (!dark || def.dark === 'same') return def;
  if (def.kind === 'theme') return def.dark ? { ...def, theme: def.dark } : def;
  if (def.kind === 'letters') return def.dark ? { ...def, ...def.dark } : def;
  if (def.kind === 'categorical') return def.dark ? { ...def, colors: def.dark as CategoricalScheme['colors'] } : def;
  return { ...def, colors: (def.dark as readonly ColorStop[] | undefined) ?? darkRamp(def.kind, def.colors) };
}

/** A letters scheme's colors, `{letter: Color}` upper case: its dictionary, its groups, its default, `overrides`. */
export function schemeLetters(
  def: LetterScheme,
  overrides: Readonly<Record<string, string>> = {},
): Record<string, string> {
  const out: Record<string, string> = {};
  const named = typeof def.colors === 'string' ? letterColors(def.colors, null) : undefined;
  for (const [k, c] of Object.entries(named ?? {})) out[k.toUpperCase()] = Color.toHex(c);
  if (isPlainObject(def.colors)) for (const [k, c] of Object.entries(def.colors)) out[k.toUpperCase()] = String(c);
  for (const grp of def.groups ?? []) for (const k of grp.letters) out[k.toUpperCase()] = grp.color;
  if (def.default != null)
    for (const a of def.alphabets)
      for (let c = 33; c < 127; c++) {
        const k = String.fromCharCode(c);
        if ((c < 97 || c > 122) && Molecules.isLetter(a, c) && !Object.hasOwn(out, k)) out[k] = def.default;
      }
  for (const [k, c] of Object.entries(overrides)) out[k.toUpperCase()] = c;
  return out;
}

/**
 * The theme a `scheme` theme stands for, its variant for the skin: `letter` with the scheme's letters, `field_palette`
 * with its ramp or categories, or a registered table's theme; null: no scheme.
 */
export function schemeTheme(reg: SvsRegistries | undefined, p: SchemeParams, o: { dark: boolean }): ThemeRef | null {
  const def = schemeOf(reg, p);
  if (!def) return null;
  const v = schemeVariant(def, o.dark);
  if (v.kind === 'theme') return v.theme;
  if (v.kind === 'letters') return { theme: 'letter', params: { colors: schemeLetters(v, p.overrides) } };
  const field = p.field ?? null;
  if (v.kind === 'categorical') {
    const palette = { kind: 'categorical', colors: v.colors, missing_color: v.missing ?? null };
    return { theme: 'field_palette', params: { field, palette } };
  }
  const palette = {
    kind: 'continuous',
    colors: v.colors,
    reverse: p.reverse ?? false,
    value_domain: p.domain ?? [null, null],
    // A value beyond the domain (a range narrower than the data) takes the ramp's end color.
    underflow_color: 'auto',
    overflow_color: 'auto',
  };
  return { theme: 'field_palette', params: { field, palette } };
}

/** What a legend lists for a scheme beside its theme's own: a letters scheme's labeled groups, a table's legend. */
export function schemeLegend(
  reg: SvsRegistries | undefined,
  p: SchemeParams,
  dark: boolean,
): { label: string; color: string }[] {
  const def = schemeOf(reg, p);
  if (!def) return [];
  const v = schemeVariant(def, dark);
  if (v.kind === 'theme') return [...(v.legend ?? [])];
  if (v.kind !== 'letters') return [];
  const letters = schemeLetters(v, p.overrides);
  const labeled = (v.groups ?? []).flatMap((x) =>
    x.label ? [{ label: x.label, color: letters[x.letters[0].toUpperCase()] ?? x.color }] : [],
  );
  if (labeled.length > 0) return labeled;
  // No labeled groups (Taylor's hue per residue, a scale's colors): the letters, those of one color together, so the
  // legend still lists the scheme (and its picker stays in reach).
  const byColor = new Map<string, string[]>();
  for (const k of Object.keys(letters).sort()) {
    if (k < 'A' || k > 'Z' || k.length !== 1) continue;
    const c = letters[k].toLowerCase();
    const same = byColor.get(c);
    if (same) same.push(k);
    else byColor.set(c, [k]);
  }
  return Array.from(byColor, ([color, ks]) => ({ label: ks.join(' '), color }));
}

// --- validation ---------------------------------------------------------------------------------------------------

const ALPHABETS: readonly Alphabet[] = ['protein', 'dna', 'rna', 'generic'];
const lettersOf = (a: Alphabet) =>
  Array.from({ length: 94 }, (_, i) => i + 33)
    .filter((c) => (c < 97 || c > 122) && Molecules.isLetter(a, c))
    .map((c) => String.fromCharCode(c))
    .join('');
/** `key` is a letter of one of `alphabets` (either case); else reported with the alphabets' letters. */
function checkLetter(key: string, alphabets: readonly Alphabet[], path: string, ctx: CheckContext) {
  const c = key.length === 1 ? key.toUpperCase().charCodeAt(0) : -1;
  if (alphabets.some((a) => Molecules.isLetter(a, c))) return;
  const which = alphabets.map((a) => `${a} (${lettersOf(a)})`).join(', ');
  report(ctx, 'svs.scheme.letter', path, `${show(key)} is not a letter of ${which || 'any alphabet given'}`);
}
function checkColor(v: unknown, path: string, ctx: CheckContext) {
  if (!Color.isValid(v)) typeError(ctx, path, 'a Color: #rgb, #rrggbb, a color name, or a/b', v, 'svs.color.invalid');
}
function checkKeys(o: Record<string, unknown>, known: readonly string[], path: string, ctx: CheckContext) {
  for (const k of Object.keys(o))
    if (!known.includes(k))
      report(ctx, 'svs.params.unknown', ptr(path, k), `unknown key "${k}"; known: ${known.join(', ')}`);
}
function checkLetterColors(
  o: Record<string, unknown>,
  alphabets: readonly Alphabet[],
  path: string,
  ctx: CheckContext,
) {
  const names = [...LetterDictionaryNames, 'ElementSymbol'];
  if (typeof o.colors === 'string' && !names.includes(o.colors))
    typeError(ctx, ptr(path, 'colors'), `a letter dictionary (${names.join(', ')}) or {letter: Color}`, o.colors);
  else if (isPlainObject(o.colors))
    for (const [k, c] of Object.entries(o.colors)) {
      checkLetter(k, alphabets, ptr(ptr(path, 'colors'), k), ctx);
      checkColor(c, ptr(ptr(path, 'colors'), k), ctx);
    }
  else if (o.colors !== undefined)
    typeError(ctx, ptr(path, 'colors'), 'a letter dictionary name or {letter: Color}', o.colors);
  if (o.groups !== undefined && !Array.isArray(o.groups))
    typeError(ctx, ptr(path, 'groups'), 'a list of {letters, color, label}', o.groups);
  (Array.isArray(o.groups) ? o.groups : []).forEach((grp, i) => {
    const at = ptr(ptr(path, 'groups'), i);
    if (!isPlainObject(grp) || typeof grp.letters !== 'string' || grp.letters.length === 0) {
      typeError(ctx, at, '{letters: "DE", color, label?}', grp);
      return;
    }
    checkKeys(grp, ['letters', 'color', 'label'], at, ctx);
    for (const k of grp.letters) checkLetter(k, alphabets, ptr(at, 'letters'), ctx);
    checkColor(grp.color, ptr(at, 'color'), ctx);
    if (grp.label !== undefined && typeof grp.label !== 'string') typeError(ctx, ptr(at, 'label'), 'str', grp.label);
  });
  if (o.default != null) checkColor(o.default, ptr(path, 'default'), ctx);
}
function checkStops(v: unknown, path: string, ctx: CheckContext) {
  if (!Array.isArray(v) || v.length < 2) return typeError(ctx, path, 'at least two colors, or [Color, t] stops', v);
  v.forEach((s, i) => {
    if (Array.isArray(s) && s.length === 2 && typeof s[1] === 'number') checkColor(s[0], ptr(ptr(path, i), 0), ctx);
    else checkColor(s, ptr(path, i), ctx);
  });
}
function checkCategories(v: unknown, path: string, ctx: CheckContext) {
  if (Array.isArray(v) && v.length > 0) for (let i = 0; i < v.length; i++) checkColor(v[i], ptr(path, i), ctx);
  else if (isPlainObject(v)) for (const [k, c] of Object.entries(v)) checkColor(c, ptr(path, k), ctx);
  else typeError(ctx, path, 'a list of colors, or {value: Color}', v);
}

/** An inline scheme (`scheme`'s `custom`): its shape, colors and letters, each fault with what to write instead. */
export function checkSchemeDef(v: unknown, path: string, ctx: CheckContext) {
  const kinds = ['letters', 'sequential', 'diverging', 'categorical'];
  if (!isPlainObject(v) || !kinds.includes(v.kind as string))
    return typeError(
      ctx,
      isPlainObject(v) ? ptr(path, 'kind') : path,
      `a scheme of kind ${kinds.join(', ')}`,
      isPlainObject(v) ? v.kind : v,
    );
  const dark = ptr(path, 'dark');
  if (v.kind === 'letters') {
    checkKeys(v, ['kind', 'alphabets', 'colors', 'groups', 'default', 'dark'], path, ctx);
    const alphabets = Array.isArray(v.alphabets) ? (v.alphabets as Alphabet[]) : [];
    if (alphabets.length === 0 || alphabets.some((a) => !ALPHABETS.includes(a)))
      return typeError(ctx, ptr(path, 'alphabets'), `a list of ${ALPHABETS.join(', ')}`, v.alphabets);
    checkLetterColors(v, alphabets, path, ctx);
    if (isPlainObject(v.dark)) {
      checkKeys(v.dark, ['colors', 'groups', 'default'], dark, ctx);
      checkLetterColors(v.dark, alphabets, dark, ctx);
    } else if (v.dark !== undefined && v.dark !== 'same')
      typeError(ctx, dark, '"same" or {colors, groups, default}', v.dark);
    return;
  }
  const ramp = v.kind !== 'categorical';
  checkKeys(v, ramp ? ['kind', 'colors', 'dark'] : ['kind', 'colors', 'missing', 'dark'], path, ctx);
  (ramp ? checkStops : checkCategories)(v.colors, ptr(path, 'colors'), ctx);
  if (!ramp && v.missing != null) checkColor(v.missing, ptr(path, 'missing'), ctx);
  if (v.dark !== undefined && v.dark !== 'same') (ramp ? checkStops : checkCategories)(v.dark, dark, ctx);
}

/**
 * `scheme`'s params across fields (`ThemeEntry.check`): exactly one of `name` (registered; an extension's unknown here
 * is drawn as the layer's default) and `custom`; `overrides` only on a letters scheme, its letters the scheme's.
 */
export function checkSchemeParams(p: Readonly<Record<string, unknown>>, path: string, ctx: CheckContext) {
  const named = typeof p.name === 'string';
  if (named && p.custom != null)
    return report(ctx, 'svs.params.exclusive', path, "give a scheme's name or a custom scheme, not both");
  if (!named && p.custom == null)
    return report(
      ctx,
      'svs.params.missing',
      ptr(path, 'name'),
      'missing "name" (a color scheme) or "custom" (an inline scheme)',
    );
  if (named && !ctx.registries.get('colorSchemes', p.name as string)) {
    const known = ctx.registries.list('colorSchemes').map((e) => e.name);
    const why = `unknown color scheme ${show(p.name)}; registered: ${known.join(', ')}`;
    if (/:/.test(p.name as string)) ctx.issues.push(Issue.warning('svs.render.theme', ptr(path, 'name'), why));
    else report(ctx, 'svs.scheme.unknown', ptr(path, 'name'), why);
    return;
  }
  const def = schemeOf(ctx.registries, p as SchemeParams);
  const overrides = isPlainObject(p.overrides) ? Object.keys(p.overrides) : [];
  if (overrides.length === 0 || !def || (def.kind === 'letters' && !Array.isArray(def.alphabets))) return;
  if (def.kind !== 'letters')
    return report(
      ctx,
      'svs.params.type',
      ptr(path, 'overrides'),
      'overrides recolor letters; for values or categories give a custom scheme',
    );
  for (const k of overrides) checkLetter(k, def.alphabets, ptr(ptr(path, 'overrides'), k), ctx);
}
