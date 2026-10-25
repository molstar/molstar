/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureElement } from '@molstar/model/model/structure';
import {
  structureAreEqual,
  structureAreIntersecting,
  structureIntersect,
  structureSubtract,
  structureUnion,
} from '@molstar/model/model/structure/query/utils/structure-set';
import type { StructureComponentParams } from '@molstar/plugin/state/helpers/structure-component';
import type { StructureComponentRef, StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import { compileTarget, fromMolqlJson, PERSON_KEYWORDS, readKeyword, type Target } from '../targets/index.js';
import { type ApiAreaContext, type ApiContext, ApiArea, ApiResult } from './area.js';
import { type ComponentParams, expressionOf, ownParams, reapplyContents, setOwnParams } from './contents.js';
import { ComponentParam, EveryStructureParam, InStructure, StructureParam } from './params.js';
import {
  componentName,
  inStructures,
  isFailure,
  labelOf,
  type StructureId,
  loadedStructures,
  someStructures,
} from './scope.js';
import { isLociTarget, LanguageParam, TargetParam, type TargetValue } from './target-param.js';

// `component` (api.md §2): the components of a structure, which are the unit the representations hang
// off. `add` makes one of what a target selects (with a representation), `modify` changes what one holds by a target,
// `remove`, `rename`, `visibility`; `list` reads them (the component manager's options are `look`'s). The scene areas share a
// shape (see `structure.ts`): a `structure` field naming one by label or ref, a component by label or ref, registry
// ids validated with the valid ones in the hint, "already so" when nothing would change. `add` follows decision 30:
// a component of exactly that target on the structure already exists, unless `copy` asks for another.

/** A representation of a component, as `list` gives it. */
export interface RepresentationInfo {
  readonly ref: string;
  /** The representation type id (`cartoon`, `ball-and-stick`, …). */
  readonly type: string;
  readonly visible: boolean;
}

/** A component, as `list` gives it. */
export interface ComponentInfo {
  readonly ref: string;
  readonly label: string;
  readonly atoms: number;
  readonly visible: boolean;
  readonly representations: readonly RepresentationInfo[];
}

/** The components of one structure. */
export interface StructureComponents {
  readonly structure: { readonly ref: string; readonly label: string };
  readonly components: readonly ComponentInfo[];
}

// The selection's words (api.md §1 "Words"): add, remove, intersect.
const ModeParam = () =>
  PD.Select(
    'add',
    [
      ['add', 'Add'],
      ['remove', 'Remove'],
      ['intersect', 'Keep only'],
    ] as const,
    { description: 'add: the target to the component; remove: take it out; intersect: keep what both hold' },
  );

export class ComponentApi extends ApiArea {
  readonly id = 'component';

  readonly add = this.operation({
    id: 'component.add',
    title: 'Add a component',
    params: {
      target: TargetParam('What the component holds: ligand, A:120-130, resn HEM, …'),
      label: PD.Text('', { description: 'Its name; empty: the target as written' }),
      representation: PD.Text('', {
        description:
          "A representation type id ('cartoon', 'ball-and-stick', …), 'none' for none; empty: ball-and-stick",
      }),
      structure: StructureParam(),
      language: LanguageParam(this.ctx.languages),
    },
    // Identity (decision 30): a component of exactly that target (or with that label) on the structure. The fast
    // path knows only the targets it can compile without the viewer (not `selection`, `focus` or loci); `run` checks
    // the rest once the target is resolved.
    exists(values, ctx) {
      const target = values.target;
      if (isLociTarget(target)) return undefined;
      const compiled = compileTarget(target as Target, { language: values.language, languages: ctx.languages });
      if (!compiled.ok) return undefined;
      const list = values.structure ? someStructures(ctx, values.structure) : loadedStructures(ctx);
      if (isFailure(list) || (!values.structure && list.length !== 1)) return undefined;
      const key = JSON.stringify(compiled.molql);
      const found = list.map((s) => findExisting(ctx, s, key));
      if (found.some((c) => !c)) return undefined;
      return `${labelOf(ctx, list[0])} already has ${list.length === 1 ? `the component '${componentName(found[0]!)}'` : 'those components'}.`;
    },
    async run({ target, label, representation, structure, language }, ctx) {
      const registry = ctx.plugin.representation.structure.registry;
      if (representation && representation !== 'none' && !registry.has(representation)) {
        return ApiResult.fail(
          'invalid',
          `There is no representation '${representation}'.`,
          `Representations: ${registry.types.map((t) => t[0]).join(', ')}, or none.`,
        );
      }
      const resolved = await ctx.targets.resolve(target, {
        ...(structure ? { structure } : {}),
        ...(language ? { language } : {}),
      });
      if (!resolved.ok) return resolved;
      const text = targetText(target);
      const parts = resolved.parts.filter((p) => p.selected.atoms > 0);
      if (!parts.length) {
        return ApiResult.fail(
          'not-found',
          `'${text}' selects nothing in ${resolved.parts.map((p) => p.structure.label).join(', ')}.`,
          'Check the target, or name the structure.',
        );
      }

      // The representation type for each part: the one asked for, which has to apply to what is selected; else
      // ball-and-stick (else the first that applies).
      const types: string[] = [];
      for (const p of parts) {
        const sub = StructureElement.Loci.toStructure(p.loci);
        if (representation && representation !== 'none') {
          if (!registry.get(representation).isApplicable(sub)) {
            return ApiResult.fail(
              'invalid',
              `'${representation}' does not apply to ${p.selected.line} in ${p.structure.label}.`,
              `Representations that do: ${
                registry
                  .getApplicableTypes(sub)
                  .map((t) => t[0])
                  .join(', ') || 'none'
              }.`,
            );
          }
          types.push(representation);
        } else if (representation === 'none') {
          types.push('none');
        } else {
          types.push(
            registry.has('ball-and-stick') && registry.get('ball-and-stick').isApplicable(sub)
              ? 'ball-and-stick'
              : (registry.getApplicableTypes(sub)[0]?.[0] ?? 'none'),
          );
        }
      }

      const person = isPerson(target);
      const stored = person || isLociTarget(target) ? undefined : fromMolqlJson(resolved.stored);
      const created: { ref: string; label: string; structure: StructureId; type: string; atoms: number }[] = [];
      const existing: string[] = [];
      // The structures by ref (a `find` a part was parts × structures), and one data transaction for all the parts:
      // the viewer syncs once at its end, not twice a structure over every structure.
      const byRef = new Map(loadedStructures(ctx).map((x) => [x.cell.transform.ref, x]));
      await ctx.plugin.dataTransaction(async () => {
        for (let i = 0; i < parts.length; i++) {
          const part = parts[i];
          const s = byRef.get(part.structure.ref);
          if (!s) continue;
          // What it holds is kept as the target's own MolQL, so the component follows a model change; the person's
          // selection and focus, and loci, are only addresses: kept as the atoms they were.
          const expression: Expression =
            stored?.ok === true
              ? stored.molql
              : StructureElement.Bundle.toExpression(StructureElement.Bundle.fromLoci(part.loci));
          if (!ctx.copy) {
            const found = findExisting(ctx, s, JSON.stringify(expression));
            if (found) {
              existing.push(`${part.structure.label} already has '${componentName(found)}'`);
              continue;
            }
          }
          const component = await ctx.plugin.builders.structure.tryCreateComponentFromExpression(
            s.cell,
            expression,
            UUID.create22(),
            { label: label || text },
          );
          if (!component) continue;
          const type = types[i];
          if (type !== 'none') await addRepresentation(ctx, component.ref, type);
          created.push({
            ref: component.ref,
            label: label || text,
            structure: part.structure,
            type,
            atoms: part.selected.atoms,
          });
        }
      });

      if (!created.length) {
        return existing.length
          ? ApiResult.already(`${existing.join('; ')}.`)
          : ApiResult.fail(
              'not-found',
              `'${text}' selects nothing in ${parts.map((p) => p.structure.label).join(', ')}.`,
            );
      }
      // What the person hid (entities, chains) stays out of the new components too (decision 18).
      const touched = new Set(created.map((c) => c.structure.ref));
      await reapplyContents(
        ctx,
        loadedStructures(ctx).filter((x) => touched.has(x.cell.transform.ref)),
      );
      const did = created
        .map(
          (c) =>
            `Added '${c.label}' (${atomsText(c.atoms)}) to ${c.structure}${c.type === 'none' ? ' without a representation' : ` as ${c.type}`}`,
        )
        .join('; ');
      return ApiResult.ok(
        `${did}.`,
        { nodes: created.map((c) => c.ref) },
        { value: created, ...(existing.length ? { warnings: existing } : {}) },
      );
    },
  });

  readonly remove = this.operation({
    id: 'component.remove',
    title: 'Remove a component',
    params: { component: ComponentParam(), structure: StructureParam(InStructure) },
    approval: 'destructive',
    async run({ component, structure }, ctx) {
      const found = inStructures(ctx, component, structure);
      if (isFailure(found)) return found;
      const name = componentName(found.component);
      // The component and the representations under it.
      await ctx.plugin.managers.structure.hierarchy.remove([found.component], false);
      return ApiResult.ok(`Removed '${name}' from ${labelOf(ctx, found.structure)}.`, {
        nodes: [found.component.cell.transform.ref],
      });
    },
  });

  readonly rename = this.operation({
    id: 'component.rename',
    title: 'Rename a component',
    params: {
      component: ComponentParam(),
      label: PD.Text('', { description: 'The new name' }),
      structure: StructureParam(InStructure),
    },
    async run({ component, label, structure }, ctx) {
      const name = label.trim();
      if (!name) return ApiResult.fail('invalid', 'The new name is empty.', 'Give the component a name.');
      const found = inStructures(ctx, component, structure);
      if (isFailure(found)) return found;
      const c = found.component;
      const old = componentName(c);
      if (old === name) return ApiResult.already(`The component is already named '${name}'.`);
      // The manager's `updateLabel`, awaited (it does not wait for its commit, and the step must see the result); its
      // own params when its structure hides contents.
      const ref = c.cell.transform.ref;
      const params = (ownParams(ctx, ref) ?? c.cell.params?.values) as StructureComponentParams;
      await setOwnParams(ctx, found.structure, ref, {
        type: params.type,
        nullIfEmpty: params.nullIfEmpty,
        label: name,
      } as ComponentParams);
      return ApiResult.ok(`Renamed '${old}' to '${name}'.`, { nodes: [ref] });
    },
  });

  readonly visibility = this.operation({
    id: 'component.visibility',
    title: 'Show or hide a component',
    params: { component: ComponentParam(), visible: PD.Boolean(true), structure: StructureParam(InStructure) },
    async run({ component, visible, structure }, ctx) {
      const found = inStructures(ctx, component, structure);
      if (isFailure(found)) return found;
      const c = found.component;
      const name = componentName(c);
      if (!!c.cell.state.isHidden === !visible) {
        return ApiResult.already(`'${name}' is already ${visible ? 'shown' : 'hidden'}.`);
      }
      // `toggleVisibility` toggles from the first component's state; it differs from the wish here.
      ctx.plugin.managers.structure.component.toggleVisibility([c]);
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} '${name}'.`, { nodes: [c.cell.transform.ref] });
    },
  });

  readonly modify = this.operation({
    id: 'component.modify',
    title: 'Change what a component holds',
    params: {
      component: ComponentParam(),
      target: TargetParam('What to add, take out or keep: ligand, A:120-130, …'),
      mode: ModeParam(),
      structure: StructureParam(InStructure),
      language: LanguageParam(this.ctx.languages),
    },
    async run({ component, target, mode, structure, language }, ctx) {
      const found = inStructures(ctx, component, structure);
      if (isFailure(found)) return found;
      const { component: c, structure: s } = found;
      const name = componentName(c);
      const current = c.cell.obj?.data;
      const parent = s.cell.obj?.data;
      if (!current || !parent) return ApiResult.fail('blocked', `'${name}' is not built yet.`);
      const resolved = await ctx.targets.resolve(target, {
        structure: s.cell.transform.ref,
        ...(language ? { language } : {}),
      });
      if (!resolved.ok) return resolved;
      const part = resolved.parts.find((p) => p.structure.ref === s.cell.transform.ref && p.selected.atoms > 0);
      if (!part) {
        return ApiResult.fail(
          'not-found',
          `'${targetText(target)}' selects nothing in ${labelOf(ctx, s)}.`,
          'Check the target.',
        );
      }

      // Mol*'s "modify by selection" (`StructureComponentManager.modifyComponent`, private): the set operation on the
      // component's atoms, then the component holds exactly those atoms.
      const by = StructureElement.Loci.toStructure(part.loci);
      const overlaps = structureAreIntersecting(current, by);
      if (mode === 'remove' && !overlaps) {
        return ApiResult.already(`'${name}' holds none of ${part.selected.line}.`);
      }
      if (mode === 'intersect' && !overlaps) {
        return ApiResult.fail(
          'invalid',
          `'${name}' and ${part.selected.line} do not overlap, so nothing would be left.`,
          "Use mode 'add' to add it, or remove the component.",
        );
      }
      const modified =
        mode === 'add'
          ? structureUnion(parent, [current, by])
          : mode === 'intersect'
            ? structureIntersect(current, by)
            : structureSubtract(current, by);
      if (structureAreEqual(modified, current)) {
        return ApiResult.already(
          mode === 'add' ? `'${name}' already holds ${part.selected.line}.` : `'${name}' is already so.`,
        );
      }
      const ref = c.cell.transform.ref;
      const own = ownParams(ctx, ref);
      if (own) {
        // Its structure hides contents: the change goes to its own params, as an expression over them, so what is
        // hidden stays out now and comes back with the rest when shown (contents.ts).
        const mine = expressionOf(own)!;
        const byExpression = StructureElement.Bundle.toExpression(StructureElement.Bundle.fromLoci(part.loci));
        const expression =
          mode === 'add'
            ? MS.struct.combinator.merge([mine, byExpression])
            : mode === 'intersect'
              ? MS.struct.modifier.intersectBy({ 0: mine, by: byExpression })
              : MS.struct.modifier.exceptBy({ 0: mine, by: byExpression });
        await setOwnParams(ctx, s, ref, { ...own, type: { name: 'expression', params: expression }, label: name });
        const now = ctx.plugin.state.data.cells.get(ref)?.obj?.data as { elementCount?: number } | undefined;
        return ApiResult.ok(
          `'${name}' now holds ${atomsText(now?.elementCount ?? 0)} shown (was ${atomsText(current.elementCount)}).`,
          { nodes: [ref] },
          { value: { atoms: now?.elementCount ?? 0 } },
        );
      }
      if (modified.elementCount === 0) {
        await ctx.plugin.managers.structure.hierarchy.remove([c], false);
        return ApiResult.ok(`'${name}' held nothing after that, so it was removed.`, {
          nodes: [c.cell.transform.ref],
        });
      }
      const params: StructureComponentParams = {
        type: { name: 'bundle', params: StructureElement.Bundle.fromSubStructure(parent, modified) },
        nullIfEmpty: true,
        label: name,
      };
      await ctx.plugin.state.data.build().to(c.cell).update(params).commit();
      return ApiResult.ok(
        `'${name}' now holds ${atomsText(modified.elementCount)} (was ${atomsText(current.elementCount)}).`,
        { nodes: [c.cell.transform.ref] },
        { value: { atoms: modified.elementCount } },
      );
    },
  });

  readonly list = this.operation({
    id: 'component.list',
    title: 'List the components',
    params: { structure: EveryStructureParam() },
    read: true,
    async run({ structure }, ctx) {
      const list = someStructures(ctx, structure || undefined);
      if (isFailure(list)) return list;
      const value: StructureComponents[] = list.map((s) => ({
        structure: { ref: s.cell.transform.ref, label: labelOf(ctx, s) },
        components: s.components.map(componentInfo),
      }));
      if (value.every((v) => !v.components.length)) return ApiResult.ok('There are no components.', {}, { value });
      const lines = value.map((v) => {
        const cs = v.components.map(
          (c) =>
            `${c.label} (${c.ref}${c.visible ? '' : ', hidden'})${c.representations.length ? ` [${c.representations.map((r) => r.type).join(', ')}]` : ''}`,
        );
        return `${v.structure.label}: ${cs.length ? cs.join('; ') : 'no components'}`;
      });
      return ApiResult.ok(lines.join('\n'), {}, { value });
    },
  });
}

function componentInfo(c: StructureComponentRef): ComponentInfo {
  return {
    ref: c.cell.transform.ref,
    label: componentName(c),
    atoms: c.cell.obj?.data.elementCount ?? 0,
    visible: !c.cell.state.isHidden,
    representations: c.representations.map((r) => ({
      ref: r.cell.transform.ref,
      type: (r.cell.transform.params as { type?: { name?: string } } | undefined)?.type?.name ?? '',
      visible: !r.cell.state.isHidden,
    })),
  };
}

/** The target as the person wrote it (the default name of a component). */
function targetText(target: TargetValue): string {
  if (isLociTarget(target)) return 'selection';
  return typeof target === 'string' ? target.trim() : JSON.stringify(target);
}

/** `selection` and `focus`: the person's, only addresses once resolved. */
function isPerson(target: TargetValue): boolean {
  if (typeof target !== 'string') return false;
  const word = readKeyword(target);
  return !!word && (PERSON_KEYWORDS as readonly string[]).includes(word);
}

function atomsText(n: number): string {
  return `${n.toLocaleString('en')} atom${n === 1 ? '' : 's'}`;
}

/** A component of the same target (decision 30: identity is the structure and the target, not the name). */
function findExisting(ctx: ApiAreaContext, s: StructureRef, key: string): StructureComponentRef | undefined {
  return s.components.find((c) => {
    // Its own params when contents are hidden (its cell's leave out what is hidden).
    const params = ownParams(ctx, c.cell.transform.ref) ?? (c.cell.params?.values as ComponentParams | undefined);
    const type = params?.type;
    return type?.name === 'expression' && JSON.stringify(type.params) === key;
  });
}

/**
 * A representation on a component with the component manager's options (hydrogens, quality, lighting, …): the manager's
 * own `addRepresentation`, on the component's hierarchy ref.
 */
async function addRepresentation(ctx: ApiContext, ref: string, type: string): Promise<void> {
  const plugin = ctx.plugin;
  const component = plugin.managers.structure.hierarchy.current.refs.get(ref);
  if (component?.kind === 'structure-component') {
    await plugin.managers.structure.component.addRepresentation([component], type);
  } else {
    await plugin.builders.structure.representation.addRepresentation(ref, {
      type: plugin.representation.structure.registry.get(type),
    });
  }
}
