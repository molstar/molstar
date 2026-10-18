/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The `espript` theme (seqspec:format#themes): ESPript's black-and-grey shading of an alignment, as a named rule table.
 * A residue whose letter is every residue of its column is `identical`; else one whose letter's group (`groups`, letters
 * read as similar) holds at least `similarity` of the column's residues is `similar`; nothing elsewhere. Shares are of
 * the column's residues (gaps are not counted), as `rules`' `share` reads them. An approximation of ESPript's
 * similarity score (Risler's matrix), by groups.
 */

/** ESPript-like groups of similar residues (the default of the theme's `groups`). */
export const ESPRIPT_GROUPS: readonly string[] = ['ILMV', 'FWY', 'KRH', 'DE', 'ST', 'NQ', 'AG'];

/** The rules the theme stands for, from its params (`identical`, `similar`, `similarity`, `groups`). */
export function espriptRules(p: Readonly<Record<string, unknown>> = {}): unknown[] {
  const identical = typeof p.identical === 'string' ? p.identical : '#000000';
  const similar = typeof p.similar === 'string' ? p.similar : '#bdbdbd';
  const t = typeof p.similarity === 'number' ? p.similarity : 0.7;
  const groups = Array.isArray(p.groups) ? (p.groups as string[]) : ESPRIPT_GROUPS;
  const out: unknown[] = [
    { when: [{ measure: { share: 'letter' }, min: 1 - 1e-9 }], color: identical, label: 'identical' },
  ];
  for (const g of groups) {
    const letters = g.toUpperCase();
    if (!letters) continue;
    const all = [...letters, ...letters.toLowerCase()];
    out.push({ when: [{ measure: { share: letters }, min: t - 1e-9, in: all }], color: similar, label: 'similar' });
  }
  return out;
}
