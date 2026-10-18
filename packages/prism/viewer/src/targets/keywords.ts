/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CommonProteinCaps, PolymerNames } from '@molstar/model/model/structure/model/types';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import type { Keyword, PersonKeyword } from './types.js';

// Ported from Stories Untold's `any-molql` (keywords.ts).

/**
 * The keywords as MolScript, for joining a keyword with other selections. The expressions are Mol*'s own
 * (`structure-selection-query` of the plugin, as MVS's static selectors use them), copied: that module imports the
 * plugin, which targets must not. `ligand` and `branched` are without the "plus connected" atoms MVS adds when it
 * draws them.
 */
const ammp = MS.ammp;
const polymerTest = (subtypes: string) =>
  MS.core.logic.and([
    MS.core.rel.eq([ammp('entityType'), 'polymer']),
    MS.core.str.match([MS.re(subtypes, 'i'), ammp('entitySubtype')]),
  ]);
const atomistic = MS.core.rel.eq([ammp('objectPrimitive'), 'atomistic']);
const nonPolymerResidue = MS.core.str.match([
  MS.re('non-polymer|(amino|carboxy) terminus|peptide-like', 'i'),
  ammp('chemCompType'),
]);
const group = (tests: Record<string, Expression>) => MS.struct.modifier.union([MS.struct.generator.atomGroups(tests)]);

const branched = group({
  'entity-test': MS.core.logic.or([
    MS.core.rel.eq([ammp('entityType'), 'branched']),
    MS.core.logic.and([
      MS.core.rel.eq([ammp('entityType'), 'non-polymer']),
      MS.core.str.match([MS.re('oligosaccharide', 'i'), ammp('entitySubtype')]),
    ]),
  ]),
});

const ligand = MS.struct.modifier.union([
  MS.struct.modifier.exceptBy({
    0: MS.struct.modifier.union([
      MS.struct.combinator.merge([
        group({
          'entity-test': MS.core.logic.and([
            MS.core.logic.or([
              MS.core.rel.eq([ammp('entityType'), 'non-polymer']),
              MS.core.rel.neq([ammp('entityPrdId'), '']),
            ]),
            MS.core.logic.not([MS.core.str.match([MS.re('(oligosaccharide|lipid|ion)', 'i'), ammp('entitySubtype')])]),
          ]),
          'chain-test': atomistic,
          'residue-test': MS.core.logic.not([MS.core.str.match([MS.re('saccharide', 'i'), ammp('chemCompType')])]),
        }),
        group({
          'entity-test': MS.core.rel.eq([ammp('entityType'), 'polymer']),
          'chain-test': atomistic,
          'residue-test': nonPolymerResidue,
        }),
      ]),
    ]),
    by: MS.struct.combinator.merge([
      group({
        'entity-test': MS.core.rel.eq([ammp('entityType'), 'polymer']),
        'chain-test': atomistic,
        'residue-test': MS.core.set.has([MS.set(...PolymerNames), ammp('label_comp_id')]),
      }),
      MS.struct.generator.atomGroups({
        'chain-test': atomistic,
        'residue-test': MS.core.set.has([MS.set(...CommonProteinCaps), ammp('label_comp_id')]),
      }),
    ]),
  }),
]);

export const KEYWORD_EXPRESSIONS: Readonly<Record<Keyword, Expression>> = {
  all: MS.struct.generator.all(),
  polymer: group({
    'entity-test': polymerTest('(polypeptide|cyclic-pseudo-peptide|peptide-like|nucleotide|peptide nucleic acid)'),
  }),
  protein: group({ 'entity-test': polymerTest('(polypeptide|cyclic-pseudo-peptide|peptide-like)') }),
  nucleic: group({ 'entity-test': polymerTest('(nucleotide|peptide nucleic acid)') }),
  branched,
  ligand,
  ion: group({ 'entity-test': MS.core.rel.eq([ammp('entitySubtype'), 'ion']) }),
  water: group({ 'entity-test': MS.core.rel.eq([ammp('entityType'), 'water']) }),
  coarse: group({ 'chain-test': MS.core.set.has([MS.set('sphere', 'gaussian'), ammp('objectPrimitive')]) }),
};

export const KEYWORDS = Object.keys(KEYWORD_EXPRESSIONS) as Keyword[];

/** The keywords that are the person's (resolved by the viewer, not here). */
export const PERSON_KEYWORDS: readonly PersonKeyword[] = ['selection', 'focus'];

export function isKeyword(text: string): text is Keyword {
  return Object.hasOwn(KEYWORD_EXPRESSIONS, text);
}

export function isPersonKeyword(text: string): text is PersonKeyword {
  return text === 'selection' || text === 'focus';
}

/**
 * A word read as a keyword: any case, surrounding spaces, and the plural a writer puts on a thing ("ligands", "ions",
 * "waters"). Undefined when it is none.
 */
export function readKeyword(text: string): Keyword | PersonKeyword | undefined {
  const word = text.trim().toLowerCase();
  if (isKeyword(word) || isPersonKeyword(word)) return word;
  if (word.endsWith('s')) {
    const singular = word.slice(0, -1);
    if (singular !== 'all' && singular !== 'coarse' && (isKeyword(singular) || isPersonKeyword(singular)))
      return singular;
  }
  return undefined;
}
