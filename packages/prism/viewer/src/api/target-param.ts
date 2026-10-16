/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { StructureElement } from '@molstar/model/model/structure';
import { defaultQueryLanguage, type QueryLanguageEntry, queryLanguages } from '../targets/languages.js';
import type { QueryLanguage, Target } from '../targets/types.js';

// The `target` param (targets.md §1): a keyword, an address list, a query text or a selector (schema items, `{ molql }`),
// and, from code in the page (the sequence link), loci directly (`{ loci }`, never stored). A plain value param: its
// command form shows a text field (`commands/operations.ts`), its tool field is the same string (agent.md §2).

/** What an operation's `target` takes. */
export type TargetValue = Target | { readonly loci: StructureElement.Loci };

export function TargetParam(description = 'What to act on: ligand, A:120-130, resn HEM, …'): PD.Value<TargetValue> {
  return PD.Value<TargetValue>('', { description });
}

/**
 * The query language of a target's text (targets.md §1): one of those this viewer reads (its spec's `languages`), the
 * first unless said; hidden where it reads none. The languages are registered when a viewer is made, after the areas'
 * params are declared, so the default, the options and whether it shows are read when asked (a form, a tool's schema,
 * a call's defaults).
 */
export function LanguageParam(): PD.Select<QueryLanguage> {
  let options: readonly (readonly [QueryLanguage, string])[] = [];
  let of: readonly QueryLanguageEntry[] | undefined;
  return Object.defineProperties(PD.Select<QueryLanguage>('pymol', []), {
    defaultValue: { get: defaultQueryLanguage, enumerable: true },
    options: {
      get() {
        const languages = queryLanguages();
        if (languages !== of) {
          of = languages;
          options = languages.map((l) => [l.id, l.name] as const);
        }
        return options;
      },
      enumerable: true,
    },
    isHidden: { get: () => queryLanguages().length === 0, enumerable: true },
  });
}

/** Whether a value names loci directly. */
export function isLociTarget(t: unknown): t is { readonly loci: StructureElement.Loci } {
  return !!t && typeof t === 'object' && 'loci' in t && !!(t as { loci: unknown }).loci;
}
