/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Loci } from '@molstar/model/model/loci';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { Selected } from '../targets/types.js';
import { summarize } from '../targets/index.js';
import { type ApiContext, type ApiFailure, ApiArea, ApiResult } from './area.js';
import { EveryStructureParam, StructureParam } from './params.js';
import { isFailure, labelOf, someStructures } from './scope.js';
import { LanguageParam, TargetParam } from './target-param.js';
import type { ResolveOptions } from './targets.js';

// `select` (operations.md §2, plan 2.3): the selection of Mol*'s structure selection manager, as view state: no history
// step (sequence.md §7.3, `history: false`). Every edit is computed per structure as the selection it should end as
// (`next`), compared with the current one ("already so" when it is the same), and written through the manager, so
// the markers, the sequence and the labels follow as for a click. Targets are exact: Mol*'s picking granularity does
// not widen them (it is for the pointer, `select.granularity`). The selection mode and granularity are the plugin's
// own state (state.md §3).

const Modes = [
  ['set', 'Replace the selection'],
  ['add', 'Add to the selection'],
  ['remove', 'Remove from the selection'],
  ['intersect', 'Keep what is in both'],
] as const satisfies readonly [string, string][];
type Mode = (typeof Modes)[number][0];

const ExpandTo = [
  ['residues', 'Whole residues'],
  ['chains', 'Whole chains'],
  ['entities', 'Whole entities'],
  ['molecules', 'Whole molecules (entities)'],
] as const satisfies readonly [string, string][];
type ExpandTarget = (typeof ExpandTo)[number][0];

/** What `select.get` gives back for one structure. */
export interface SelectionInStructure {
  readonly structure: { readonly ref: string; readonly label: string };
  readonly selected: Selected;
}

/** One structure's selection: now, and as it should be. */
interface Edit {
  readonly current: StructureElement.Loci;
  readonly next: StructureElement.Loci;
}

export class SelectApi extends ApiArea {
  readonly id = 'select';

  readonly set = this.operation({
    id: 'select.set',
    title: 'Select',
    params: {
      target: TargetParam('What to select: ligand, A:120-130, resn HEM, …'),
      mode: PD.Select<Mode>('set', Modes as unknown as [Mode, string][], {
        description: 'Replace the selection, add to it, remove from it or keep what is in both',
      }),
      structure: StructureParam(),
      language: LanguageParam(),
    },
    history: false,
    async run({ target, mode, structure, language }, ctx) {
      const resolved = await ctx.targets.resolve(target, resolveOptions(structure, language));
      if (isFailure(resolved)) return resolved;
      if (resolved.empty) return ApiResult.fail('not-found', resolved.summary, 'Nothing was changed.');
      const selection = ctx.plugin.managers.structure.selection;
      const edits = resolved.parts.map((p) => {
        const current = currentOf(selection, p.loci.structure);
        return { current, next: combine(mode, current, p.loci) };
      });
      if (mode === 'set' && !structure) {
        // The selection becomes the target: structures it does not reach lose theirs.
        for (const s of ctx.plugin.managers.structure.hierarchy.current.structures) {
          const data = s.cell.obj?.data;
          if (!data || resolved.parts.some((p) => p.loci.structure.root === data.root)) continue;
          const current = currentOf(selection, data);
          if (!StructureElement.Loci.isEmpty(current)) edits.push({ current, next: StructureElement.Loci.none(data) });
        }
      }
      if (!write(ctx.plugin, edits)) {
        return ApiResult.already(
          mode === 'remove'
            ? 'None of that is selected.'
            : mode === 'intersect'
              ? 'The selection is already within that.'
              : 'That is already selected.',
        );
      }
      const did =
        mode === 'set'
          ? `Selected ${resolved.summary}.`
          : mode === 'add'
            ? `Added ${resolved.summary} to the selection.`
            : mode === 'remove'
              ? `Removed ${resolved.summary} from the selection.`
              : `Narrowed the selection to its overlap with ${resolved.summary}.`;
      return ApiResult.ok(did, { selection: true });
    },
  });

  readonly all = this.operation({
    id: 'select.all',
    title: 'Select everything',
    params: { structure: EveryStructureParam() },
    history: false,
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ structure }, ctx) {
      const targets = structuresOf(ctx, structure);
      if (isFailure(targets)) return targets;
      const edits = targets.map((s) => ({ current: currentOf(selection(ctx), s), next: StructureElement.Loci.all(s) }));
      if (!write(ctx.plugin, edits)) return ApiResult.already('Everything is already selected.');
      return ApiResult.ok(`Selected ${describe(edits)}.`, { selection: true });
    },
  });

  readonly clear = this.operation({
    id: 'select.clear',
    title: 'Clear the selection',
    params: {},
    history: false,
    async run(_, ctx) {
      if (selection(ctx).elementCount() === 0) return ApiResult.already('Nothing is selected.');
      ctx.plugin.managers.interactivity.lociSelects.deselectAll();
      return ApiResult.ok('Selection cleared.', { selection: true });
    },
  });

  readonly invert = this.operation({
    id: 'select.invert',
    title: 'Invert the selection',
    params: { structure: EveryStructureParam() },
    history: false,
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ structure }, ctx) {
      const targets = structuresOf(ctx, structure);
      if (isFailure(targets)) return targets;
      const edits = targets.map((s) => {
        const current = currentOf(selection(ctx), s);
        return { current, next: StructureElement.Loci.subtract(StructureElement.Loci.all(s), current) };
      });
      if (!write(ctx.plugin, edits)) return ApiResult.already('The structure is empty.');
      return ApiResult.ok(`Inverted the selection: ${describe(edits)} are selected.`, { selection: true });
    },
  });

  readonly around = this.operation({
    id: 'select.around',
    title: 'Select the surroundings',
    params: {
      radius: PD.Numeric(5, { min: 0, max: 100, step: 0.5 }, { description: 'Distance from the selection, in Å' }),
      asWholeResidues: PD.Boolean(true, { description: 'Include the residues that reach into the radius whole' }),
      structure: EveryStructureParam(),
    },
    history: false,
    when: nothingSelected,
    async run({ radius, asWholeResidues, structure }, ctx) {
      const targets = selectedIn(ctx, structure);
      if (isFailure(targets)) return targets;
      const edits = targets.map((current) => {
        const around = StructureElement.Loci.extendToRadius(current, radius);
        return { current, next: asWholeResidues ? StructureElement.Loci.extendToWholeResidues(around) : around };
      });
      if (!write(ctx.plugin, edits))
        return ApiResult.already(`The selection already includes what is within ${radius} Å.`);
      return ApiResult.ok(`Selection grown by ${radius} Å: ${describe(edits)}.`, { selection: true });
    },
  });

  readonly expand = this.operation({
    id: 'select.expand',
    title: 'Expand the selection',
    params: {
      to: PD.Select<ExpandTarget>('residues', ExpandTo as unknown as [ExpandTarget, string][], {
        description: 'Extend the selection to whole residues, chains or entities',
      }),
      structure: EveryStructureParam(),
    },
    history: false,
    when: nothingSelected,
    async run({ to, structure }, ctx) {
      const targets = selectedIn(ctx, structure);
      if (isFailure(targets)) return targets;
      const edits = targets.map((current) => ({ current, next: expand(current, to) }));
      if (!write(ctx.plugin, edits)) return ApiResult.already(`The selection already covers whole ${to}.`);
      return ApiResult.ok(`Selection expanded to whole ${to}: ${describe(edits)}.`, { selection: true });
    },
  });

  readonly mode = this.operation({
    id: 'select.mode',
    title: 'Selection mode',
    params: { on: PD.Boolean(true, { description: 'Clicks in the canvas select instead of focus' }) },
    history: false,
    async run({ on }, ctx) {
      if (ctx.plugin.selectionMode === on) {
        return ApiResult.already(`Selection mode is already ${on ? 'on' : 'off'}.`);
      }
      ctx.plugin.selectionMode = on;
      return ApiResult.ok(`Selection mode ${on ? 'on' : 'off'}.`, { view: true });
    },
  });

  readonly granularity = this.operation({
    id: 'select.granularity',
    title: 'Picking level',
    params: {
      granularity: PD.Select<Loci.Granularity>('residue', Loci.GranularityOptions, {
        description: 'What a click in the canvas selects: an atom, a residue, a chain, …',
      }),
    },
    history: false,
    async run({ granularity }, ctx) {
      const interactivity = ctx.plugin.managers.interactivity;
      if (interactivity.props.granularity === granularity) {
        return ApiResult.already(`The picking level is already ${granularity}.`);
      }
      interactivity.setProps({ granularity });
      return ApiResult.ok(`Picking level set to ${granularity}.`, { view: true });
    },
  });

  readonly get = this.operation({
    id: 'select.get',
    title: 'Describe the selection',
    params: { structure: EveryStructureParam() },
    read: true,
    async run({ structure }, ctx) {
      const targets = someStructures(ctx, structure || undefined);
      if (isFailure(targets)) return targets;
      const value: SelectionInStructure[] = [];
      for (const s of targets) {
        const data = s.cell.obj?.data;
        if (!data) continue;
        const current = currentOf(selection(ctx), data);
        if (StructureElement.Loci.isEmpty(current)) continue;
        value.push({ structure: { ref: s.cell.transform.ref, label: labelOf(ctx, s) }, selected: summarize(current) });
      }
      if (!value.length) return ApiResult.ok('Nothing is selected.', {}, { value });
      const many = targets.length > 1;
      const text = value.map((v) => (many ? `${v.structure.label}: ${v.selected.line}` : v.selected.line));
      return ApiResult.ok(`Selected: ${text.join('; ')}.`, {}, { value });
    },
  });
}

function resolveOptions(structure: string, language: ResolveOptions['language']): ResolveOptions {
  return { ...(structure ? { structure } : {}), ...(language ? { language } : {}) };
}

const selection = (ctx: ApiContext) => ctx.plugin.managers.structure.selection;

/** A structure's selection (an empty one where it has none). */
function currentOf(
  manager: PluginContext['managers']['structure']['selection'],
  structure: Structure,
): StructureElement.Loci {
  const loci = manager.getLoci(structure);
  return StructureElement.Loci.is(loci) ? loci : StructureElement.Loci.none(structure);
}

/** The data of the named structure, or of every loaded one. */
function structuresOf(ctx: ApiContext, name: string): Structure[] | ApiFailure {
  const list = someStructures(ctx, name || undefined);
  if (isFailure(list)) return list;
  return list.flatMap((s) => (s.cell.obj ? [s.cell.obj.data] : []));
}

/** The selections of the named structure, or of every one that has some. */
function selectedIn(ctx: ApiContext, name: string): StructureElement.Loci[] | ApiFailure {
  const list = structuresOf(ctx, name);
  if (isFailure(list)) return list;
  const found = list.map((s) => currentOf(selection(ctx), s)).filter((l) => !StructureElement.Loci.isEmpty(l));
  return found.length ? found : ApiResult.fail('blocked', 'Nothing is selected.', 'Select something first.');
}

function nothingSelected(ctx: { plugin: PluginContext }): true | string {
  return ctx.plugin.managers.structure.selection.elementCount() > 0 ? true : 'Nothing is selected.';
}

function combine(mode: Mode, current: StructureElement.Loci, target: StructureElement.Loci): StructureElement.Loci {
  switch (mode) {
    case 'add':
      return StructureElement.Loci.union(current, target);
    case 'remove':
      return StructureElement.Loci.subtract(current, target);
    case 'intersect':
      return StructureElement.Loci.intersect(current, target);
    default:
      return target;
  }
}

function expand(loci: StructureElement.Loci, to: ExpandTarget): StructureElement.Loci {
  switch (to) {
    case 'residues':
      return StructureElement.Loci.extendToWholeResidues(loci);
    case 'chains':
      return StructureElement.Loci.extendToWholeChains(loci);
    default:
      // A molecule is an entity here (`project`'s "molecules", targets.md §2).
      return StructureElement.Loci.extendToWholeEntities(loci);
  }
}

/** The same elements (sizes are equal and one holds the other; the order of units does not matter). */
function same(a: StructureElement.Loci, b: StructureElement.Loci): boolean {
  const size = StructureElement.Loci.size(a);
  if (size !== StructureElement.Loci.size(b)) return false;
  return size === 0 || StructureElement.Loci.isSubset(a, b);
}

/**
 * Writes the edits that change something through the selection manager (without Mol*'s picking granularity: the loci
 * are exact). Returns whether anything changed.
 */
function write(plugin: PluginContext, edits: readonly Edit[]): boolean {
  const manager = plugin.managers.structure.selection;
  let changed = false;
  for (const { current, next } of edits) {
    if (same(current, next)) continue;
    changed = true;
    // An empty loci is no selection to set: what was selected is removed.
    if (StructureElement.Loci.isEmpty(next)) manager.fromLoci('remove', current, false);
    else manager.fromLoci('set', next, false);
  }
  return changed;
}

/** What the edits leave selected, in words. */
function describe(edits: readonly Edit[]): string {
  const lines = edits.filter((e) => !StructureElement.Loci.isEmpty(e.next)).map((e) => summarize(e.next).line);
  return lines.length ? lines.join('; ') : 'nothing';
}
