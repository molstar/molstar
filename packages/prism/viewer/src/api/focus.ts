/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureElement } from '@molstar/model/model/structure';
import type { FocusEntry } from '@molstar/plugin/state/manager/structure/focus';
import type { Selected } from '../targets/types.js';
import { summarize } from '../targets/index.js';
import { ApiArea, ApiResult } from './area.js';
import { StructureParam } from './params.js';
import { isFailure, type StructureId, structureId } from './scope.js';
import { LanguageParam, TargetParam } from './target-param.js';
import { duration as cameraDuration } from './camera.js';
import type { ApiAreaContext } from './area.js';

// `focus` (api.md §2): Mol*'s focus (`managers.structure.focus`), the one part of a structure the
// views dwell on, and, where the host has a canvas, the camera moves to it. View state: no history step
// (sequence.md §7.3). The focus is the plugin's own state (state.md §3), so a headless viewer has one too.

/** What `focus.get` gives back. */
export interface FocusValue {
  readonly structure: StructureId;
  readonly label: string;
  readonly selected: Selected;
}

/** One of the recent foci (`focus.history`): a target `focus.set` takes back, with the structure it is on. */
export interface FocusHistoryEntry {
  readonly label: string;
  /** The structure (its `ref` is what `focus.set`'s `structure` takes). */
  readonly structure: StructureId;
  /** The focus as addresses (`A:120-130, A:135`): residue-level, so an atom-level focus comes back as its residues. */
  readonly target: string;
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
      if (canvas) ctx.plugin.managers.camera.focusLoci(part.loci, { durationMs: cameraDuration(ctx) });
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
        structure: cell ? structureId(ctx, cell) : { ref: '', label: '' },
        label: ctx.plugin.managers.structure.focus.current?.label ?? selected.line,
        selected,
      };
      return ApiResult.ok(`Focused on ${selected.line}.`, {}, { value });
    },
  });

  // The focus history (parity.md §3): Mol*'s own (`focus.history`, newest first, eight at most), as v5's focus picker
  // lists it. Focusing again moves an entry to the top; view state like the rest of `focus`.

  readonly history = this.operation({
    id: 'focus.history',
    title: 'List the recent foci',
    params: {},
    read: true,
    async run(_, ctx) {
      const entries = recentFoci(ctx).map(({ entry, structure }): FocusHistoryEntry => {
        const selected = summarize(entry.loci, { maxRangesPerChain: Number.MAX_SAFE_INTEGER });
        return { label: entry.label, structure, target: selected.addresses.join(', '), selected };
      });
      if (!entries.length) return ApiResult.ok('Nothing was focused yet.', {}, { value: entries });
      const list = entries.map((e, i) => `${i + 1}. ${e.label}`).join('; ');
      return ApiResult.ok(`${entries.length} recent foci, newest first: ${list}.`, {}, { value: entries });
    },
  });

  readonly again = this.operation({
    id: 'focus.again',
    title: 'Focus again',
    params: (ctx) => ({
      entry: PD.Select<string>(defaultEntry(ctx), historyOptions(ctx), {
        description: 'One of the recent foci (`focus.history`), by its place in the list, newest first: 1, 2, …',
      }),
    }),
    history: false,
    when: (ctx) => (recentFoci(ctx).length ? true : 'Nothing was focused yet.'),
    async run({ entry }, ctx) {
      const list = recentFoci(ctx);
      const picked = list[Number(entry) - 1];
      if (!picked)
        return ApiResult.fail('not-found', `There is no recent focus ${entry}.`, `There are ${list.length}.`);
      const focus = ctx.plugin.managers.structure.focus;
      const canvas = ctx.modules.has('canvas');
      const current = focus.current?.loci;
      if (!canvas && current && StructureElement.Loci.areEqual(current, picked.entry.loci)) {
        return ApiResult.already(`Already focused on ${picked.entry.label}.`);
      }
      // The entry as Mol* keeps it (its label and category), not re-resolved: it moves to the top of the history.
      focus.set(picked.entry);
      if (canvas) ctx.plugin.managers.camera.focusLoci(picked.entry.loci, { durationMs: cameraDuration(ctx) });
      return ApiResult.ok(`Focused on ${picked.entry.label}.`, { view: true });
    },
  });
}

/** The plugin's focus history with each entry's structure; entries whose structure is gone are left out. */
function recentFoci(ctx: ApiAreaContext): { entry: FocusEntry; structure: StructureId }[] {
  const structures = ctx.plugin.managers.structure.hierarchy.current.structures;
  const out: { entry: FocusEntry; structure: StructureId }[] = [];
  for (const entry of ctx.plugin.managers.structure.focus.history) {
    const cell = structures.find((s) => s.cell.obj?.data.root === entry.loci.structure.root);
    if (cell) out.push({ entry, structure: structureId(ctx, cell) });
  }
  return out;
}

function historyOptions(ctx: ApiAreaContext): [string, string][] {
  const list = recentFoci(ctx);
  if (!list.length) return [['1', 'Nothing was focused yet']];
  const several = new Set(list.map((e) => e.structure.ref)).size > 1;
  return list.map(({ entry, structure }, i) => {
    const on = several ? ` (${structure.label})` : '';
    return [`${i + 1}`, `${entry.label}${on}`];
  });
}

/** The newest entry that is not the current focus (the current one is usually first). */
function defaultEntry(ctx: ApiAreaContext): string {
  const list = recentFoci(ctx);
  const current = ctx.plugin.managers.structure.focus.current?.loci;
  const i = list.findIndex((e) => !current || !StructureElement.Loci.areEqual(current, e.entry.loci));
  return `${i < 0 ? 1 : i + 1}`;
}
