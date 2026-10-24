/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureElement } from '@molstar/model/model/structure';
import type { Mark, MarkTone } from '../sources/marks.js';
import { ApiArea, ApiResult } from './area.js';
import { StructureParam } from './params.js';
import { isFailure, loadedStructures } from './scope.js';
import { isLociTarget, LanguageParam, TargetParam } from './target-param.js';

// `highlight` (api.md §2): the hover mark, for embedding hosts and the sequence link; with `sequence`,
// lit in the sequence too (a view other than the sequence pointing at residues: a plot's cell). And the marks
// (marks.md): an owner's named residue ranges in a tone, in the 3D view and the sequence until cleared. View state: no
// history step (sequence.md §7.3). It runs at pointer rate from the sequence, so loci from code in the page go
// straight to the highlight manager (no resolving, no summary), and the result says little.

/** One mark as `highlight.mark` takes it. */
export interface MarkInput {
  readonly id: string;
  readonly label?: string;
  readonly tone: MarkTone;
  /** Its residues: a target (`A:92-190`, `ligand`). */
  readonly target: string;
}

const HighlightOwner = 'highlight';

export class HighlightApi extends ApiArea {
  readonly id = 'highlight';

  readonly set = this.operation({
    id: 'highlight.set',
    title: 'Highlight',
    surface: { kind: 'ui-only', reason: 'Hover and the sequence link light things up; nothing to run on purpose.' },
    params: {
      target: TargetParam('What to highlight: ligand, A:120-130, resn HEM, …'),
      structure: StructureParam(),
      language: LanguageParam(this.ctx.languages),
      sequence: PD.Boolean(false, {
        isHidden: true,
        description: 'Light it in the sequence too (a view other than the sequence pointing at it)',
      }),
    },
    history: false,
    async run({ target, structure, language, sequence }, ctx) {
      const highlights = ctx.plugin.managers.interactivity.lociHighlights;
      if (isLociTarget(target)) {
        if (!StructureElement.Loci.is(target.loci) || StructureElement.Loci.isEmpty(target.loci)) {
          return ApiResult.fail('not-found', 'There is nothing to highlight.');
        }
        highlights.highlightOnly({ loci: target.loci }, false);
        const root = target.loci.structure.root;
        const ref = sequence
          ? loadedStructures(ctx).find((s) => s.cell.obj?.data.root === root)?.cell.transform.ref
          : undefined;
        ctx.sources.marks.setHover(HighlightOwner, ref ? { structure: ref, loci: target.loci } : undefined);
        return ApiResult.ok('Highlighted.', { view: true });
      }
      const resolved = await ctx.targets.resolve(target, {
        ...(structure ? { structure } : {}),
        ...(language ? { language } : {}),
      });
      if (isFailure(resolved)) return resolved;
      if (resolved.empty) return ApiResult.fail('not-found', resolved.summary, 'Nothing was changed.');
      // The sequence follows this call: a hover it lit before is not left behind.
      if (!sequence) ctx.sources.marks.setHover(HighlightOwner, undefined);
      let first = true;
      for (const p of resolved.parts) {
        if (p.selected.atoms === 0) continue;
        if (first) highlights.highlightOnly({ loci: p.loci }, false);
        else highlights.highlight({ loci: p.loci }, false);
        if (first && sequence) ctx.sources.marks.setHover(HighlightOwner, { structure: p.structure.ref, loci: p.loci });
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
      ctx.sources.marks.setHover(HighlightOwner, undefined);
      return ApiResult.ok('Highlight cleared.', { view: true });
    },
  });

  /**
   * An owner's marks (marks.md): residue ranges in a tone, in the 3D view and the sequence, replacing the ones it had,
   * until `highlight.unmark`. What a view shows of its own state on the scene (the PAE plot's box): view state.
   */
  readonly mark = this.operation({
    id: 'highlight.mark',
    title: 'Mark residue ranges',
    surface: {
      kind: 'ui-only',
      reason: "A view's own marks on the scene (a plot's box); the agent highlights or selects.",
    },
    params: {
      owner: PD.Text('', { description: 'Who marks (its marks are replaced and cleared together)' }),
      marks: PD.Value<readonly MarkInput[]>([], {
        description: 'The marks: id, label, tone (highlight, accent), target',
      }),
      structure: StructureParam(),
      language: LanguageParam(this.ctx.languages),
      dimRest: PD.Boolean(false, { description: 'Grey the rest of the structure while the marks are shown' }),
      focusSequence: PD.Boolean(false, {
        description: 'Bring the marks into view in the sequence, all of them together (it zooms out to show them)',
      }),
    },
    history: false,
    async run({ owner, marks, structure, language, dimRest, focusSequence }, ctx) {
      if (!owner) return ApiResult.fail('invalid', 'Who marks?', 'Give an owner.');
      const out: Omit<Mark, 'owner'>[] = [];
      for (const m of marks) {
        const resolved = await ctx.targets.resolve(m.target, {
          ...(structure ? { structure } : {}),
          ...(language ? { language } : {}),
        });
        if (isFailure(resolved)) return resolved;
        if (resolved.empty) return ApiResult.fail('not-found', `${m.id}: ${resolved.summary}`, 'Nothing was marked.');
        for (const p of resolved.parts) {
          if (p.selected.atoms === 0) continue;
          out.push({
            id: m.id,
            label: m.label ?? m.target,
            tone: m.tone,
            structure: p.structure.ref,
            loci: p.loci,
            dimRest,
            ...(focusSequence ? { focusSequence } : {}),
          });
        }
      }
      ctx.sources.marks.set(owner, out);
      return ApiResult.ok(out.length ? `Marked ${out.map((m) => m.label).join('; ')}.` : 'Nothing to mark.', {
        view: true,
      });
    },
  });

  readonly unmark = this.operation({
    id: 'highlight.unmark',
    title: "Clear an owner's marks",
    surface: { kind: 'ui-only', reason: "A view's own marks on the scene (a plot's box)." },
    params: { owner: PD.Text('', { description: 'Whose marks' }) },
    history: false,
    async run({ owner }, ctx) {
      ctx.sources.marks.clear(owner);
      return ApiResult.ok('Marks cleared.', { view: true });
    },
  });
}
