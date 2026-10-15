/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { Selected } from '../targets/types.js';
import { summarize } from '../targets/index.js';
import { ApiArea, ApiResult } from './area.js';
import { StructureParam } from './params.js';
import { isFailure } from './scope.js';
import { LanguageParam, TargetParam } from './target-param.js';

// `focus` (operations.md §2, plan 2.3): Mol*'s focus (`managers.structure.focus`), the one part of a structure the
// views dwell on, and, where the host has a canvas, the camera moves to it. View state: no history step
// (sequence.md §7.3). The focus is the plugin's own state (state.md §3), so a headless viewer has one too.

/** What `focus.get` gives back. */
export interface FocusValue {
  readonly structure: string;
  readonly label: string;
  readonly selected: Selected;
}

export class FocusApi extends ApiArea {
  readonly id = 'focus';

  readonly set = this.operation({
    id: 'focus.set',
    title: 'Focus',
    params: {
      target: TargetParam('What to focus on: ligand, A:120-130, resn HEM, …'),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    history: false,
    async run({ target, structure, language }, ctx) {
      const resolved = await ctx.targets.resolve(target, {
        ...(structure ? { structure } : {}),
        ...(language ? { language } : {}),
      });
      if (isFailure(resolved)) return resolved;
      if (resolved.empty) return ApiResult.fail('not-found', resolved.summary, 'Nothing was changed.');
      const parts = resolved.parts.filter((p) => p.selected.atoms > 0);
      const part = parts[0];
      const canvas = ctx.modules.has('canvas');
      const current = ctx.plugin.managers.structure.focus.current?.loci;
      if (!canvas && current && StructureElement.Loci.areEqual(current, part.loci)) {
        return ApiResult.already(`Already focused on ${part.selected.line}.`);
      }
      const focus = ctx.plugin.managers.structure.focus;
      focus.setFromLoci(part.loci);
      if (canvas) ctx.plugin.managers.camera.focusLoci(part.loci);
      return ApiResult.ok(
        `Focused on ${part.selected.line}.`,
        { view: true },
        parts.length > 1
          ? {
              warnings: [
                `The target is on ${parts.length} structures; the focus is one part: ${part.structure.label}.`,
              ],
            }
          : {},
      );
    },
  });

  readonly clear = this.operation({
    id: 'focus.clear',
    title: 'Clear the focus',
    params: {},
    history: false,
    async run(_, ctx) {
      const focus = ctx.plugin.managers.structure.focus;
      if (!focus.current) return ApiResult.already('Nothing is focused.');
      focus.clear();
      return ApiResult.ok('Focus cleared.', { view: true });
    },
  });

  readonly get = this.operation({
    id: 'focus.get',
    title: 'Describe the focus',
    params: {},
    read: true,
    async run(_, ctx) {
      const loci = ctx.plugin.managers.structure.focus.current?.loci;
      if (!loci || StructureElement.Loci.isEmpty(loci))
        return ApiResult.ok('Nothing is focused.', {}, { value: undefined });
      const cell = ctx.plugin.managers.structure.hierarchy.current.structures.find(
        (s) => s.cell.obj?.data.root === loci.structure.root,
      );
      const selected = summarize(loci);
      const value: FocusValue = {
        structure: cell?.cell.transform.ref ?? '',
        label: ctx.plugin.managers.structure.focus.current?.label ?? selected.line,
        selected,
      };
      return ApiResult.ok(`Focused on ${selected.line}.`, {}, { value });
    },
  });
}
