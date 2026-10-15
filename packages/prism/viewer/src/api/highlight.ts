/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import { ApiArea, ApiResult } from './area.js';
import { StructureParam } from './params.js';
import { isFailure } from './scope.js';
import { isLociTarget, LanguageParam, TargetParam } from './target-param.js';

// `highlight` (operations.md §2, plan 2.3): the hover mark, for embedding hosts and the sequence link. View state: no
// history step (sequence.md §7.3). It runs at pointer rate from the sequence, so loci from code in the page go
// straight to the highlight manager (no resolving, no summary), and the result says little.

export class HighlightApi extends ApiArea {
  readonly id = 'highlight';

  readonly set = this.operation({
    id: 'highlight.set',
    title: 'Highlight',
    surface: { kind: 'ui-only', reason: 'Hover and the sequence link light things up; nothing to run on purpose.' },
    params: {
      target: TargetParam('What to highlight: ligand, A:120-130, resn HEM, …'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    history: false,
    async run({ target, structure, language }, ctx) {
      const highlights = ctx.plugin.managers.interactivity.lociHighlights;
      if (isLociTarget(target)) {
        if (!StructureElement.Loci.is(target.loci) || StructureElement.Loci.isEmpty(target.loci)) {
          return ApiResult.fail('not-found', 'There is nothing to highlight.');
        }
        highlights.highlightOnly({ loci: target.loci }, false);
        return ApiResult.ok('Highlighted.', { view: true });
      }
      const resolved = await ctx.targets.resolve(target, {
        ...(structure ? { structure } : {}),
        ...(language ? { language } : {}),
      });
      if (isFailure(resolved)) return resolved;
      if (resolved.empty) return ApiResult.fail('not-found', resolved.summary, 'Nothing was changed.');
      let first = true;
      for (const p of resolved.parts) {
        if (p.selected.atoms === 0) continue;
        if (first) highlights.highlightOnly({ loci: p.loci }, false);
        else highlights.highlight({ loci: p.loci }, false);
        first = false;
      }
      return ApiResult.ok(`Highlighted ${resolved.summary}.`, { view: true });
    },
  });

  readonly clear = this.operation({
    id: 'highlight.clear',
    title: 'Clear the highlight',
    surface: { kind: 'ui-only', reason: 'Hover and the sequence link light things up; nothing to run on purpose.' },
    params: {},
    history: false,
    async run(_, ctx) {
      ctx.plugin.managers.interactivity.lociHighlights.clearHighlights();
      return ApiResult.ok('Highlight cleared.', { view: true });
    },
  });
}
