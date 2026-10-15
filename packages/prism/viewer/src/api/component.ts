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
import { StructureComponentManager } from '@molstar/plugin/state/manager/structure/component';
import type { StructureComponentRef, StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type { Expression } from '@molstar/query-language/language/expression';
import { compileTarget, fromMolqlJson, PERSON_KEYWORDS, readKeyword, type Target } from '../targets/index.js';
import { type ApiContext, type ApiFailure, ApiArea, ApiResult } from './area.js';
import { reapplyContents } from './contents.js';
import { isFailure, labelOf, loadedStructures, someStructures } from './scope.js';
import { isLociTarget, LanguageParam, TargetParam, type TargetValue } from './target-param.js';

// `component` (operations.md §2, plan 2.2): the components of a structure, which are the unit the representations hang
// off. `add` makes one of what a target selects (with a representation), `modify` changes what one holds by a target,
// `remove`, `rename`, `visibility`, and the component manager's `options`; `list` reads them. The scene areas share a
// shape (see `structure.ts`): a `structure` field naming one by label or ref, a component by label or ref, registry
// ids validated with the valid ones in the hint, "already so" when nothing would change. `add` follows decision 30:
// a component of exactly that target on the structure already exists, unless `copy` asks for another.

const StructureParam = () =>
  PD.Text('', { description: 'A structure by label or state ref; empty: the one in view (or all, where it says)' });

const ComponentParam = () =>
  PD.Text('', { description: 'A component by label or state ref (within the structure, when it is named)' });

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

/** How the component manager's options are named here (the manager's names where they were long or odd). */
const OptionKeys = {
  hydrogens: 'hydrogens',
  quality: 'visualQuality',
  ignoreLight: 'ignoreLight',
  solidSurface: 'solidSurface',
  material: 'materialStyle',
  clip: 'clipObjects',
  interior: 'interior',
  animation: 'animation',
} as const satisfies Record<string, keyof StructureComponentManager.Options>;

const ModeParam = () =>
  PD.Select(
    'union',
    [
      ['union', 'Add'],
      ['subtract', 'Remove'],
      ['intersect', 'Keep only'],
    ] as const,
    { description: 'union: add the target to the component; subtract: take it out; intersect: keep what both hold' },
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
      language: LanguageParam(),
    },
    // Identity (decision 30): a component of exactly that target (or with that label) on the structure. The fast
    // path knows only the targets it can compile without the viewer (not `selection`, `focus` or loci); `run` checks
    // the rest once the target is resolved.
    exists(values, ctx) {
      const target = values.target;
      if (isLociTarget(target)) return undefined;
      const compiled = compileTarget(target as Target, { language: values.language });
      if (!compiled.ok) return undefined;
      const list = values.structure ? someStructures(ctx, values.structure) : loadedStructures(ctx);
      if (isFailure(list) || (!values.structure && list.length !== 1)) return undefined;
      const key = JSON.stringify(compiled.molql);
      const found = list.map((s) => findExisting(s, key, values.label));
      if (found.some((c) => !c)) return undefined;
      return `${labelOf(ctx, list[0])} already has ${list.length === 1 ? `the component '${nameOf(found[0]!)}'` : 'those components'}.`;
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
      const created: { ref: string; label: string; structure: string; type: string; atoms: number }[] = [];
      const existing: string[] = [];
      for (let i = 0; i < parts.length; i++) {
        const part = parts[i];
        const s = loadedStructures(ctx).find((x) => x.cell.transform.ref === part.structure.ref);
        if (!s) continue;
        // What it holds is kept as the target's own MolQL, so the component follows a model change; the person's
        // selection and focus, and loci, are only addresses: kept as the atoms they were.
        const expression: Expression =
          stored?.ok === true
            ? stored.molql
            : StructureElement.Bundle.toExpression(StructureElement.Bundle.fromLoci(part.loci));
        if (!ctx.copy) {
          const found = findExisting(s, JSON.stringify(expression), label);
          if (found) {
            existing.push(`${part.structure.label} already has '${nameOf(found)}'`);
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
          structure: part.structure.label,
          type,
          atoms: part.selected.atoms,
        });
      }

      if (!created.length) {
        return existing.length
          ? ApiResult.already(`${existing.join('; ')}.`)
          : ApiResult.fail(
              'not-found',
              `'${text}' selects nothing in ${parts.map((p) => p.structure.label).join(', ')}.`,
            );
      }
      // What the person hid (entities, chains) stays out of the new components too (decision 18).
      const touched = new Set(created.map((c) => c.structure));
      await reapplyContents(
        ctx,
        loadedStructures(ctx).filter((x) => touched.has(labelOf(ctx, x))),
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
    params: { component: ComponentParam(), structure: StructureParam() },
    async run({ component, structure }, ctx) {
      const found = findComponent(ctx, component, structure);
      if (isFailure(found)) return found;
      const name = nameOf(found.component);
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
      structure: StructureParam(),
    },
    async run({ component, label, structure }, ctx) {
      const name = label.trim();
      if (!name) return ApiResult.fail('invalid', 'The new name is empty.', 'Give the component a name.');
      const found = findComponent(ctx, component, structure);
      if (isFailure(found)) return found;
      const c = found.component;
      if (nameOf(c) === name) return ApiResult.already(`The component is already named '${name}'.`);
      // The manager's `updateLabel`, awaited: it does not wait for its commit, and the history step must see the result.
      const params = c.cell.params?.values as StructureComponentParams;
      await ctx.plugin.state.data
        .build()
        .to(c.cell)
        .update({ type: params.type, nullIfEmpty: params.nullIfEmpty, label: name } satisfies StructureComponentParams)
        .commit();
      return ApiResult.ok(`Renamed '${nameOf(c)}' to '${name}'.`, { nodes: [c.cell.transform.ref] });
    },
  });

  readonly visibility = this.operation({
    id: 'component.visibility',
    title: 'Show or hide a component',
    params: { component: ComponentParam(), visible: PD.Boolean(true), structure: StructureParam() },
    async run({ component, visible, structure }, ctx) {
      const found = findComponent(ctx, component, structure);
      if (isFailure(found)) return found;
      const c = found.component;
      const name = nameOf(c);
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
      structure: StructureParam(),
      language: LanguageParam(),
    },
    async run({ component, target, mode, structure, language }, ctx) {
      const found = findComponent(ctx, component, structure);
      if (isFailure(found)) return found;
      const { component: c, structure: s } = found;
      const name = nameOf(c);
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
      if (mode === 'subtract' && !overlaps) {
        return ApiResult.already(`'${name}' holds none of ${part.selected.line}.`);
      }
      if (mode === 'intersect' && !overlaps) {
        return ApiResult.fail(
          'invalid',
          `'${name}' and ${part.selected.line} do not overlap, so nothing would be left.`,
          "Use mode 'union' to add it, or remove the component.",
        );
      }
      const modified =
        mode === 'union'
          ? structureUnion(parent, [current, by])
          : mode === 'intersect'
            ? structureIntersect(current, by)
            : structureSubtract(current, by);
      if (structureAreEqual(modified, current)) {
        return ApiResult.already(
          mode === 'union' ? `'${name}' already holds ${part.selected.line}.` : `'${name}' is already so.`,
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

  readonly options = this.operation({
    id: 'component.options',
    title: 'Set how representations are made',
    // Defaults are the current values, so what is not given stays as it is.
    params: (ctx) => optionParams(ctx.plugin.managers.structure.component.state.options),
    async run(values, ctx) {
      const manager = ctx.plugin.managers.structure.component;
      const current = manager.state.options;
      const next: StructureComponentManager.Options = { ...current };
      const changed: string[] = [];
      for (const [name, key] of Object.entries(OptionKeys) as [keyof typeof OptionKeys, keyof typeof current][]) {
        const param = StructureComponentManager.OptionsParams[key] as PD.Any;
        if (PD.isParamEqual(param, current[key], values[name])) continue;
        (next as Record<string, unknown>)[key] = values[name];
        changed.push(name);
      }
      if (!changed.length) return ApiResult.already('Those options are already set.');
      await manager.setOptions(next);
      return ApiResult.ok(`Set ${changed.join(', ')}.`, {}, { value: changed });
    },
  });

  readonly list = this.operation({
    id: 'component.list',
    title: 'List the components',
    params: { structure: StructureParam() },
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

/** The options as params whose defaults are the current values. */
function optionParams(current: StructureComponentManager.Options) {
  const P = StructureComponentManager.OptionsParams;
  const withDefault = <T extends PD.Any>(p: T, value: unknown): T => ({ ...p, defaultValue: value }) as T;
  return {
    hydrogens: withDefault(P.hydrogens, current.hydrogens),
    quality: withDefault(P.visualQuality, current.visualQuality),
    ignoreLight: withDefault(P.ignoreLight, current.ignoreLight),
    solidSurface: withDefault(P.solidSurface, current.solidSurface),
    material: withDefault(P.materialStyle, current.materialStyle),
    clip: withDefault(P.clipObjects, current.clipObjects),
    interior: withDefault(P.interior, current.interior),
    animation: withDefault(P.animation, current.animation),
  };
}

/** The component's name for people. */
function nameOf(c: StructureComponentRef): string {
  const params = c.cell.params?.values as StructureComponentParams | undefined;
  return params?.label || c.cell.obj?.label || c.key || c.cell.transform.ref;
}

function componentInfo(c: StructureComponentRef): ComponentInfo {
  return {
    ref: c.cell.transform.ref,
    label: nameOf(c),
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

/** A component of exactly that MolQL (as JSON), or with that name when one is given. */
function findExisting(s: StructureRef, key: string, label: string): StructureComponentRef | undefined {
  return s.components.find((c) => {
    if (label && nameOf(c) === label.trim()) return true;
    const type = (c.cell.params?.values as StructureComponentParams | undefined)?.type;
    return type?.name === 'expression' && JSON.stringify(type.params) === key;
  });
}

/** One component by state ref or name (the case may differ), in the named structure or in all of them. */
function findComponent(
  ctx: ApiContext,
  name: string,
  structure: string,
): { structure: StructureRef; component: StructureComponentRef } | ApiFailure {
  const structures = someStructures(ctx, structure || undefined);
  if (isFailure(structures)) return structures;
  const all = structures.flatMap((s) => s.components.map((component) => ({ structure: s, component })));
  const describe = (x: (typeof all)[number]) =>
    `${nameOf(x.component)} (${x.component.cell.transform.ref}${structures.length > 1 ? `, in ${labelOf(ctx, x.structure)}` : ''})`;
  const list = all.length ? `Components: ${all.map(describe).join(', ')}.` : 'There are no components.';
  if (!name) {
    if (all.length === 1) return all[0];
    return ApiResult.fail('ambiguous', 'Which component?', list);
  }
  const byRef = all.filter((x) => x.component.cell.transform.ref === name);
  const matches = byRef.length
    ? byRef
    : (() => {
        const exact = all.filter((x) => nameOf(x.component) === name);
        return exact.length ? exact : all.filter((x) => nameOf(x.component).toLowerCase() === name.toLowerCase());
      })();
  if (matches.length === 1) return matches[0];
  if (!matches.length) return ApiResult.fail('not-found', `There is no component '${name}'.`, list);
  return ApiResult.fail(
    'ambiguous',
    `Which '${name}'? ${matches.length} components have that name.`,
    `Give its state ref${structure ? '' : ' or name the structure'}: ${matches.map(describe).join(', ')}.`,
  );
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
