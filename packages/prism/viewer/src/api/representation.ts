/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Color } from '@molstar/core/util/color';
import { decodeColor } from '@molstar/core/util/color/utils';
import { Material } from '@molstar/core/util/material';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { setSubtreeVisibility } from '@molstar/plugin/behavior/static/state';
import { PluginCommands } from '@molstar/plugin/commands';
import { isSurfaceRepresentationType } from '@molstar/plugin/state/helpers/structure-representation-params';
import type {
  StructureComponentRef,
  StructureRef,
  StructureRepresentationRef,
} from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { setStructureClipping } from '@molstar/plugin/state/helpers/structure-clipping';
import { setStructureEmissive } from '@molstar/plugin/state/helpers/structure-emissive';
import { setStructureOverpaint } from '@molstar/plugin/state/helpers/structure-overpaint';
import { setStructureTransparency } from '@molstar/plugin/state/helpers/structure-transparency';
import { type ApiAreaContext, type ApiContext, ApiArea, type ApiFailure, ApiResult } from './area.js';
import {
  ComponentParam,
  EveryStructureParam,
  InStructure,
  mergeParams,
  StructureParam,
  OptionalParam,
  NodeParam,
} from './params.js';
import { type DataPalette, dataPaletteOf, setDataPalette } from '../sources/palette.js';
import { nodesFor, updateEach } from './nodes.js';
import {
  componentName,
  findComponent,
  inStructures,
  isFailure,
  labelOf,
  loadedStructures,
  someStructures,
} from './scope.js';
import {
  readStyleLayers,
  StyleLayerCells,
  type StyleLayerInfo,
  type StyleLayerKind,
  StyleLayerKinds,
  styleLayerCellsOf,
  styleLayerId,
  styleLayerKey,
} from '../sources/style-layers.js';
import { setRepresentation } from './set-representation.js';
import { LanguageParam, TargetParam } from './target-param.js';

// `representation` (api.md §2): what a component shows and how it looks. `set` and `remove` act on a
// component's representations (by type); the look operations take a target (targets.md §1–2) and act, per structure,
// on the representations of the components the target covers: `color`, `transparency`, `emissive` and `clipping` are
// Mol*'s layers (overpaint, transparency, emissive, clipping on a representation), which hold the target's elements;
// `color` with a theme, `size` and `material` change the representation's own params. Transparency is the overlay from
// opaque (0) to hidden (1), keeping the geometry (README decision 18). Desired state: a call that would change nothing
// says so (`already`, no history step).

// The layers (the plugin's `state/helpers/structure-*.ts`): which transformer holds them under a representation, and
// the tag those helpers give its cells (they keep it private); `sources/style-layers.ts` reads them for people.
const Layers = StyleLayerCells;

type LayerKind = StyleLayerKind;

/** The Mol* clipping groups, by name. */
const ClippingNames = Object.keys(Clipping.Groups.Names) as Clipping.Groups.Names[];

export class RepresentationApi extends ApiArea {
  readonly id = 'representation';

  // accessibility.md: the colors the categorical themes (chains, entities, …) give by default: Prism's, Mol*'s classic
  // or colorblind-safe; the representations on the old default follow (viewer-spec.md §2.8).
  readonly defaultPalette = this.operation({
    id: 'representation.defaultPalette',
    title: 'Set the default palette',
    params: {
      palette: PD.Select<DataPalette>(
        'colorblind',
        [
          ['prism', 'Prism'],
          ['classic', 'Classic (Mol*)'],
          ['colorblind', 'Colorblind-safe'],
        ],
        { description: 'The default colors of chains, entities, polymers and other categories' },
      ),
    },
    async run({ palette }, ctx) {
      if (dataPaletteOf(ctx.plugin) === palette) return ApiResult.already(`The default palette is already ${palette}.`);
      const nodes = (await setDataPalette(ctx.plugin, palette)) ?? [];
      return ApiResult.ok(
        palette === 'colorblind'
          ? 'Chains and entities now use colorblind-safe colors.'
          : palette === 'classic'
            ? 'Chains and entities now use Mol*’s classic colors.'
            : 'Chains and entities now use Prism’s colors.',
        { nodes },
      );
    },
  });

  // The inspector's form: the whole params of the representations a node selector matches (`#ref` for a
  // row's); equal params: already.
  readonly update = this.operation({
    id: 'representation.update',
    title: 'Update a representation',
    params: {
      node: NodeParam('Representations by selector (component[label=Ligand] > representation) or a state ref'),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The representation's params (type, color theme, size theme)" },
      ),
    },
    surface: { kind: 'ui-only', reason: "The inspector's form; tools use set, color and size" },
    async run({ node, params }, ctx) {
      const found = representationsOf(ctx, node);
      if (isFailure(found)) return found;
      return updateEach(ctx, found, params, 'Representation');
    },
  });

  // A representation row's eye (the scene tree): the matched cells and what hangs off them (their layers).
  readonly visibility = this.operation({
    id: 'representation.visibility',
    title: 'Show or hide a representation',
    params: {
      node: NodeParam('Representations by selector (representation[type=cartoon]) or a state ref'),
      visible: PD.Boolean(true),
    },
    async run({ node, visible }, ctx) {
      const found = representationsOf(ctx, node);
      if (isFailure(found)) return found;
      const todo = found.filter((r) => !!r.cell.state.isHidden === visible);
      if (!todo.length) {
        return ApiResult.already(
          `${representationsName(found)} ${found.length > 1 ? 'are' : 'is'} already ${visible ? 'shown' : 'hidden'}.`,
        );
      }
      for (const r of todo) setSubtreeVisibility(ctx.plugin.state.data, refOf(r), !visible);
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${representationsName(todo)}.`, { nodes: todo.map(refOf) });
    },
  });

  // api.md §2: on a component, the representation of a type: added if absent, its params updated if present.
  readonly set = this.operation({
    id: 'representation.set',
    title: 'Show a component as…',
    params: (ctx) => {
      const types = ctx.plugin.representation.structure.registry.types;
      return {
        component: ComponentParam(),
        type: PD.Select(types[0]?.[0] ?? '', types, {
          description: 'A representation type: cartoon, ball-and-stick, gaussian-surface, …',
        }),
        params: PD.Value<Record<string, unknown>>({}, { description: "Some of the type's params, by name" }),
        structure: StructureParam(InStructure),
      };
    },
    async run({ component, type, params, structure }, ctx) {
      const found = inStructures(ctx, component, structure);
      if (isFailure(found)) return found;
      const c = found.component;
      return setRepresentation(ctx, {
        name: componentName(c),
        type,
        patch: params,
        registry: ctx.plugin.representation.structure.registry,
        themes: ctx.plugin.representation.structure.themes,
        data: c.cell.obj?.data,
        existing: c.representations,
        // Added the way the component manager adds one (its options: hydrogens, quality, material, …), then the params.
        async create(provider) {
          const selector = await ctx.plugin.builders.structure.representation.addRepresentation(c.cell, {
            type: provider as never,
            typeParams: mergeParams(managerTypeParams(ctx, type), params ?? {}) as never,
          });
          return [selector ? selector.ref : c.cell.transform.ref];
        },
      });
    },
  });

  // api.md §2: removes that type (or all) from the component; or the representations a node selector matches
  // (a tree row's by `#ref`, which another of the same type may sit beside).
  readonly remove = this.operation({
    id: 'representation.remove',
    title: 'Remove a representation',
    params: {
      component: ComponentParam(),
      type: PD.Text('', { description: "A representation type; empty: the component's only one" }),
      all: OptionalParam(PD.Boolean(false, { description: 'Every representation of the component' })),
      structure: StructureParam(InStructure),
      node: NodeParam('Representations by selector or a state ref (instead of a component and type)'),
    },
    approval: 'destructive',
    async run({ component, type, all: every, structure, node }, ctx) {
      if (node) {
        const found = representationsOf(ctx, node);
        if (isFailure(found)) return found;
        await ctx.plugin.managers.structure.hierarchy.remove(found, false);
        return ApiResult.ok(`Removed ${representationsName(found)}.`, { nodes: found.map(refOf) });
      }
      const found = inStructures(ctx, component, structure);
      if (isFailure(found)) return found;
      const c = found.component;
      if (type && !ctx.plugin.representation.structure.registry.has(type)) {
        return ApiResult.fail('invalid', `There is no representation type '${type}'.`, `Types: ${typeIds(ctx)}.`);
      }
      // An empty remove is the only one (api.md §1 "Words"); `all` clears.
      if (!type && !every && c.representations.length > 1) {
        return ApiResult.fail(
          'ambiguous',
          `Which representation of ${componentName(c)}? It has ${c.representations.length}.`,
          `Give a type (${[...new Set(c.representations.map(typeOf))].join(', ')}), a node, or all: true.`,
        );
      }
      const doomed = c.representations.filter((r) => !type || typeOf(r) === type);
      if (!doomed.length) {
        return ApiResult.already(
          type ? `${componentName(c)} does not show ${type}.` : `${componentName(c)} shows nothing already.`,
        );
      }
      await ctx.plugin.managers.structure.hierarchy.remove(doomed, false);
      return ApiResult.ok(`Removed ${[...new Set(doomed.map(typeOf))].join(', ')} from ${componentName(c)}.`, {
        nodes: doomed.map((r) => r.cell.transform.ref),
      });
    },
  });

  // api.md §2: a color on what the target selects (the overpaint layer), or a color theme on the
  // representations of the components it covers.
  readonly color = this.operation({
    id: 'representation.color',
    title: 'Color',
    params: (ctx) => {
      const themes = ctx.plugin.representation.structure.themes.colorThemeRegistry.types;
      return {
        target: TargetParam(),
        color: PD.Text('', {
          description: "A CSS color: red, #ff8800, rgb(255, 136, 0); with theme 'uniform', the theme's color",
        }),
        theme: PD.Select('', [['', 'None'], ...themes.map((t) => [t[0], t[1]] as [string, string])], {
          description: 'A color theme: element-symbol, chain-id, uncertainty, …',
        }),
        structure: StructureParam(),
        language: LanguageParam(this.ctx.languages),
      };
    },
    async run({ target, color, theme, structure, language }, ctx) {
      // The uniform theme takes its color (a color alone paints the overpaint layer).
      const uniform = theme === 'uniform' && !!color;
      if (!!color === !!theme && !uniform) {
        return ApiResult.fail(
          'invalid',
          color ? 'Give a color or a theme, not both.' : 'Give a color or a theme.',
          'color: a CSS color painted on the target; theme: a color theme for its representations.',
        );
      }
      if (theme) {
        const registry = ctx.plugin.representation.structure.themes.colorThemeRegistry;
        if (!registry.has(theme)) {
          return ApiResult.fail(
            'invalid',
            `There is no color theme '${theme}'.`,
            `Themes: ${registry.types.map((t) => t[0]).join(', ')}.`,
          );
        }
      }
      const painted = color ? decodeColor(color) : undefined;
      if (color && painted === undefined) {
        return ApiResult.fail(
          'invalid',
          `'${color}' is not a color.`,
          'Give a CSS color: a name (red), #rrggbb, or rgb(r, g, b).',
        );
      }

      const where = await resolveCovered(ctx, target, structure, language);
      if (isFailure(where)) return where;

      if (theme) return applyTheme(ctx, where, { color: theme, ...(uniform ? { value: painted } : {}) });
      return applyLayer(ctx, where, 'color', `Colored ${where.summary} ${color}.`, (part, getLoci) =>
        setStructureOverpaint(ctx.plugin, part.components, painted!, getLoci),
      );
    },
  });

  // api.md §2, decision 18: transparency from opaque (0) to hidden (1), keeping the geometry.
  readonly transparency = this.operation({
    id: 'representation.transparency',
    title: 'Set transparency',
    params: {
      target: TargetParam(),
      alpha: PD.Numeric(
        0.5,
        { min: 0, max: 1, step: 0.01 },
        { description: '0: opaque, 1: hidden (the geometry stays)' },
      ),
      structure: StructureParam(),
      language: LanguageParam(this.ctx.languages),
    },
    async run({ target, alpha, structure, language }, ctx) {
      const where = await resolveCovered(ctx, target, structure, language);
      if (isFailure(where)) return where;
      return applyLayer(
        ctx,
        where,
        'transparency',
        alpha === 0 ? `Made ${where.summary} opaque.` : `Set the transparency of ${where.summary} to ${alpha}.`,
        (part, getLoci) => setStructureTransparency(ctx.plugin, part.components, alpha, getLoci),
        alpha === 0,
      );
    },
  });

  // api.md §2: the emissive layer.
  readonly emissive = this.operation({
    id: 'representation.emissive',
    title: 'Set emissive strength',
    params: {
      target: TargetParam(),
      strength: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { description: '0: none, 1: the full color glows' }),
      structure: StructureParam(),
      language: LanguageParam(this.ctx.languages),
    },
    async run({ target, strength, structure, language }, ctx) {
      const where = await resolveCovered(ctx, target, structure, language);
      if (isFailure(where)) return where;
      return applyLayer(
        ctx,
        where,
        'emissive',
        strength === 0
          ? `Made ${where.summary} not glow.`
          : `Set the emissive strength of ${where.summary} to ${strength}.`,
        (part, getLoci) => setStructureEmissive(ctx.plugin, part.components, strength, getLoci),
        strength === 0,
      );
    },
  });

  // api.md §2: the clipping layer. Mol*'s layer names the clip objects (canvas clipping, groups one to six) the
  // target is excluded from: by default all of them, so the target stays whole while the rest is clipped.
  readonly clipping = this.operation({
    id: 'representation.clipping',
    title: 'Exclude from clipping',
    params: {
      target: TargetParam(),
      groups: PD.MultiSelect<Clipping.Groups.Names>(
        [...ClippingNames],
        ClippingNames.map((n) => [n, n] as [Clipping.Groups.Names, string]),
        { description: 'The clip objects the target is excluded from (clipped by none of these); default all' },
      ),
      structure: StructureParam(),
      language: LanguageParam(this.ctx.languages),
    },
    async run({ target, groups, structure, language }, ctx) {
      const names = groups ?? ClippingNames;
      const unknown = names.filter((n) => !Clipping.Groups.isName(n));
      if (unknown.length) {
        return ApiResult.fail(
          'invalid',
          `There are no clip groups ${unknown.join(', ')}.`,
          `Groups: ${ClippingNames.join(', ')}.`,
        );
      }
      const where = await resolveCovered(ctx, target, structure, language);
      if (isFailure(where)) return where;
      const flags = Clipping.Groups.fromNames([...names]);
      return applyLayer(
        ctx,
        where,
        'clipping',
        `Excluded ${where.summary} from clip ${names.length === 1 ? 'object' : 'objects'} ${names.join(', ') || '(none)'}.`,
        (part, getLoci) => setStructureClipping(ctx.plugin, part.components, flags, getLoci),
      );
    },
  });

  // api.md §2: the size theme and the size factor on the representations the target covers (all, when empty).
  readonly size = this.operation({
    id: 'representation.size',
    title: 'Set size',
    params: (ctx) => {
      const themes = ctx.plugin.representation.structure.themes.sizeThemeRegistry.types;
      return {
        target: TargetParam('What to resize: ligand, A:120-130, …; empty: everything'),
        // Each absent: kept (api.md §1); one of them is needed.
        theme: OptionalParam(
          PD.Select(
            themes[0]?.[0] ?? 'uniform',
            themes.map((t) => [t[0], t[1]] as [string, string]),
            {
              description: 'A size theme: uniform, physical, …',
            },
          ),
        ),
        factor: OptionalParam(
          PD.Numeric(
            1,
            { min: 0.05, max: 10, step: 0.05 },
            { description: 'A multiple of the default size (1: the default, 2: twice as big)' },
          ),
        ),
        structure: EveryStructureParam(),
        language: LanguageParam(this.ctx.languages),
      };
    },
    async run({ target, theme, factor, structure, language }, ctx) {
      if (!theme && !factor) {
        return ApiResult.fail(
          'invalid',
          'Give a size theme or a factor.',
          'theme: a size theme; factor: a multiple of the default size.',
        );
      }
      if (theme) {
        const registry = ctx.plugin.representation.structure.themes.sizeThemeRegistry;
        if (!registry.has(theme)) {
          return ApiResult.fail(
            'invalid',
            `There is no size theme '${theme}'.`,
            `Themes: ${registry.types.map((t) => t[0]).join(', ')}.`,
          );
        }
      }
      const where = target
        ? await resolveCovered(ctx, target, structure, language)
        : await allCovered(ctx, structure || undefined);
      if (isFailure(where)) return where;

      const before = signature(ctx, where.reprs);
      const warnings: string[] = [];
      await ctx.plugin.dataTransaction(async () => {
        if (theme) await changeTheme(ctx, where, { size: theme }, warnings);
        if (factor) await changeSizeFactor(ctx, where, factor, warnings);
      });
      const touched = touchedText(where);
      if (before === signature(ctx, where.reprs)) {
        return warnings.length
          ? ApiResult.fail('blocked', warnings.join(' '))
          : ApiResult.already(`${touched.text} already has that size.`);
      }
      if (touched.warning) warnings.push(touched.warning);
      const did = [theme ? `the ${theme} size theme` : '', factor ? `size ×${factor}` : '']
        .filter(Boolean)
        .join(' and ');
      return ApiResult.ok(
        `Set ${did} on ${touched.text}.`,
        { nodes: where.reprs.map(refOf) },
        warnings.length ? { warnings } : {},
      );
    },
  });

  // api.md §2: the material params (metalness, roughness, bumpiness) of representations.
  readonly material = this.operation({
    id: 'representation.material',
    title: 'Set material',
    params: {
      material: Material.getParam({ isFlat: true }),
      component: ComponentParam('A component by label or state ref; empty: all of the structure'),
      structure: EveryStructureParam(),
    },
    async run({ material, component, structure }, ctx) {
      const targets = someStructures(ctx, structure || undefined);
      if (isFailure(targets)) return targets;
      const reprs: StructureRepresentationRef[] = [];
      let matched = 0;
      for (const s of targets) {
        if (component) {
          // A component only some of the structures have is fine; only a name that is unclear in one is a failure.
          const c = findComponent(ctx, [s], component);
          if (isFailure(c)) {
            if (c.code !== 'not-found') return c;
            continue;
          }
          matched++;
          reprs.push(...c.component.representations);
        } else {
          for (const c of s.components) reprs.push(...c.representations);
        }
      }
      if (component && !matched) {
        const none = findComponent(ctx, targets, component);
        return isFailure(none) ? none : ApiResult.fail('not-found', `There is no component '${component}'.`);
      }
      if (!reprs.length) return ApiResult.fail('not-found', 'There are no representations to change.');

      const patch = { ...Material.Zero, ...material };
      const changing = reprs.filter((r) => {
        const p = (r.cell.transform.params as RepresentationParams).type.params;
        return p.material && !Material.areEqual(p.material as Material, patch);
      });
      const supporting = reprs.filter(
        (r) => 'material' in (r.cell.transform.params as RepresentationParams).type.params,
      );
      if (!supporting.length) return ApiResult.fail('blocked', 'None of these representations has a material.');
      if (!changing.length) return ApiResult.already(`Their material is already ${Material.toString(patch)}.`);
      // One state update for all of them (an update each was a commit each: representations × the tree).
      const state = ctx.plugin.state.data;
      const tree = state.build();
      for (const r of changing) {
        const old = r.cell.transform.params as RepresentationParams;
        tree.to(r.cell.transform.ref).update({
          ...old,
          type: { name: old.type.name, params: { ...old.type.params, material: { ...patch } } },
        });
      }
      await ctx.plugin.dataTransaction(() => PluginCommands.State.Update(ctx.plugin, { state, tree }));
      return ApiResult.ok(
        `Set the material of ${changing.length} representation${changing.length === 1 ? '' : 's'} (${Material.toString(patch)}).`,
        {
          nodes: changing.map(refOf),
        },
      );
    },
  });

  // api.md §2: clears a layer: every representation of the structure (or all) loses what it holds; `all`: every
  // kind at once (the Selection tab's Clear all), one step.
  readonly reset = this.operation({
    id: 'representation.reset',
    title: 'Reset a style layer',
    params: {
      layer: PD.Select<LayerKind | 'all'>('color', [
        ['color', 'Color'],
        ['transparency', 'Transparency'],
        ['emissive', 'Emissive'],
        ['substance', 'Material'],
        ['clipping', 'Clipping'],
        ['all', 'All layers'],
      ]),
      structure: EveryStructureParam(),
    },
    async run({ layer, structure }, ctx) {
      if (layer !== 'all' && !(layer in Layers)) {
        return ApiResult.fail(
          'invalid',
          `There is no layer '${layer}'.`,
          `Layers: ${Object.keys(Layers).join(', ')}, all.`,
        );
      }
      const targets = someStructures(ctx, structure || undefined);
      if (isFailure(targets)) return targets;
      const reprs = targets.flatMap((s) => s.components.flatMap((c) => c.representations));
      const kinds = layer === 'all' ? StyleLayerKinds : [layer];
      const cells = kinds.flatMap((k) => layerCells(ctx, k, reprs));
      const names = targets.map((x) => labelOf(ctx, x)).join(', ');
      const what = layer === 'all' ? 'style layers' : Layers[layer].name;
      if (!cells.length) return ApiResult.already(`${names} has no ${what} to reset.`);
      const update = ctx.plugin.state.data.build();
      for (const c of cells) update.delete(c.ref);
      await update.commit({ doNotUpdateCurrent: true });
      return ApiResult.ok(`Reset the ${what} of ${names}.`, {
        nodes: [...new Set(cells.map((c) => c.parent))],
      });
    },
  });

  // The applied layers, one per value of a kind on a structure (`sources/style-layers.ts`): what the Selection tab
  // lists, and the ids `layer-update` and `layer-remove` take.
  readonly layers = this.operation({
    id: 'representation.layers',
    title: 'List the style layers',
    params: { structure: EveryStructureParam() },
    read: true,
    async run({ structure }, ctx) {
      const targets = someStructures(ctx, structure || undefined);
      if (isFailure(targets)) return targets;
      const refs = new Set(targets.map((t) => t.cell.transform.ref));
      const list = readStyleLayers(ctx.plugin).filter((l) => refs.has(l.structure.ref));
      if (!list.length) {
        return ApiResult.ok(
          `${targets.map((t) => labelOf(ctx, t)).join(', ')} has no style layers.`,
          {},
          { value: list },
        );
      }
      return ApiResult.ok(list.map((l) => `${l.id}: ${layerLine(l)}`).join('\n'), {}, { value: list });
    },
  });

  // One layer's value changed where it is (a color, a transparency, an emissive strength): the Selection tab's edit.
  readonly layerUpdate = this.operation({
    id: 'representation.layer-update',
    title: 'Change a style layer',
    params: (ctx) => ({
      layer: LayerParam(ctx, (l) => l.kind === 'color' || l.kind === 'transparency' || l.kind === 'emissive'),
      value: PD.Text('', {
        description: 'The new value: a CSS color for a color layer; 0 to 1 for transparency (1: hidden) and emissive',
      }),
    }),
    when: (ctx) => (ctx.sources.layers.value.length ? true : 'No style layers are applied.'),
    async run({ layer, value }, ctx) {
      const found = findLayer(ctx, layer);
      if (isFailure(found)) return found;
      const { kind } = found;
      let patch: { color: Color } | { value: number };
      if (kind === 'color') {
        if (found.clear)
          return ApiResult.fail('invalid', 'A color reset has no color to change.', 'Remove it instead.');
        const color = decodeColor(value);
        if (color === undefined) {
          return ApiResult.fail(
            'invalid',
            `'${value}' is not a color.`,
            'Give a CSS color: a name (red), #rrggbb, or rgb(r, g, b).',
          );
        }
        patch = { color };
      } else if (kind === 'transparency' || kind === 'emissive') {
        const n = Number(value);
        if (!value.trim() || !Number.isFinite(n) || n < 0 || n > 1) {
          return ApiResult.fail('invalid', `'${value}' is not a number from 0 to 1.`);
        }
        patch = { value: n };
      } else {
        return ApiResult.fail(
          'invalid',
          `A ${Layers[kind].name} layer is not changed in place.`,
          'Remove it and apply another to the target.',
        );
      }
      const next = styleLayerKey(kind, patch);
      const id = styleLayerId(kind, found.structure.ref, next);
      if (next === found.value) return ApiResult.already(`It is already ${next}.`, { value: { id } });
      const touched = await editLayerCells(ctx, found, (l) => ({ ...l, ...patch }));
      return ApiResult.ok(
        `Changed the ${Layers[kind].name} of ${found.where} (${found.structure.label}) from ${found.value} to ${next}.`,
        { nodes: touched },
        { value: { id } },
      );
    },
  });

  // One layer removed from every representation of its structure (the rest of the kind stays): the Selection tab's ✕.
  readonly layerRemove = this.operation({
    id: 'representation.layer-remove',
    title: 'Remove a style layer',
    params: (ctx) => ({ layer: LayerParam(ctx) }),
    when: (ctx) => (ctx.sources.layers.value.length ? true : 'No style layers are applied.'),
    async run({ layer }, ctx) {
      const found = findLayer(ctx, layer);
      if (isFailure(found)) return found;
      const touched = await editLayerCells(ctx, found, () => undefined);
      return ApiResult.ok(`Removed ${layerLine(found)}.`, { nodes: touched });
    },
  });
}

// --- the layers by id -------------------------------------------------------------------------------------------

/** A layer for people: "color #e0523a on A:10-17 of 1CBS". */
function layerLine(l: StyleLayerInfo): string {
  const what =
    l.kind === 'color'
      ? l.clear
        ? 'the color reset'
        : `color ${l.value}`
      : l.kind === 'substance'
        ? l.clear
          ? 'the material reset'
          : `material ${l.value}`
        : l.kind === 'clipping'
          ? `the clipping exclusion (${l.value})`
          : `${Layers[l.kind].name} ${l.value}`;
  return `${what} on ${l.where} of ${l.structure.label}`;
}

/** A layer by id, its choices the layers applied now (those `keep` keeps). Left out, it is the only layer there is. */
function LayerParam(ctx: ApiAreaContext, keep: (l: StyleLayerInfo) => boolean = () => true) {
  const list = ctx.sources.layers.value.filter(keep);
  return PD.Select<string>('', [['', 'Choose one'], ...list.map((l) => [l.id, layerLine(l)] as [string, string])], {
    description: 'A style layer by id, as representation.layers lists them; empty: the only one',
  });
}

function findLayer(ctx: ApiContext, id: string): StyleLayerInfo | ApiFailure {
  const list = readStyleLayers(ctx.plugin);
  // No id: the only layer, else the call does not say which (an operation that removes or recolors never guesses).
  if (!id && list.length === 1) return list[0];
  if (!id && list.length > 1) {
    return ApiResult.fail(
      'ambiguous',
      `Which layer? ${list.length} are applied.`,
      `Name one: ${list.map((l) => l.id).join(', ')}.`,
    );
  }
  const found = list.find((l) => l.id === id);
  if (found) return found;
  return ApiResult.fail(
    'not-found',
    `There is no style layer '${id}'.`,
    list.length ? `Layers: ${list.map((l) => l.id).join(', ')}.` : 'No style layers are applied.',
  );
}

/**
 * Each cell that holds the layer gets it changed (`edit` gives the new one) or removed (`edit` gives undefined); a cell
 * left without layers is deleted. One tree update; the representations touched are returned.
 */
async function editLayerCells(
  ctx: ApiContext,
  layer: StyleLayerInfo,
  edit: (l: Record<string, unknown>) => Record<string, unknown> | undefined,
): Promise<string[]> {
  const state = ctx.plugin.state.data;
  const update = state.build();
  const touched = new Set<string>();
  for (const ref of layer.cells) {
    const transform = state.tree.transforms.get(ref);
    if (!transform) continue;
    const params = transform.params as { layers: Record<string, unknown>[] };
    const next: Record<string, unknown>[] = [];
    for (const l of params.layers) {
      if (styleLayerKey(layer.kind, l) !== layer.value) next.push(l);
      else {
        const edited = edit(l);
        if (edited) next.push(edited);
      }
    }
    if (next.length) update.to(ref).update({ ...params, layers: next });
    else update.delete(ref);
    touched.add(transform.parent);
  }
  await update.commit({ doNotUpdateCurrent: true });
  return [...touched];
}

// --- what the operations act on ---------------------------------------------------------------------------------

interface RepresentationParams {
  type: { name: string; params: Record<string, any> };
  colorTheme: { name: string; params: Record<string, any> };
  sizeTheme: { name: string; params: Record<string, any> };
}

const typeOf = (r: StructureRepresentationRef) => (r.cell.transform.params as RepresentationParams).type.name;
const refOf = (r: StructureRepresentationRef) => r.cell.transform.ref;

/** The scene's structure representations a node selector matches (api.md §1 "Nodes"). */
function representationsOf(ctx: ApiContext, node: string): StructureRepresentationRef[] | ApiFailure {
  const byRef = new Map<string, StructureRepresentationRef>();
  for (const s of loadedStructures(ctx)) {
    for (const c of s.components) for (const r of c.representations) byRef.set(refOf(r), r);
  }
  return nodesFor(ctx.plugin.state.data, node, 'representation', (cell) => byRef.get(cell.transform.ref));
}

/** "cartoon of Polymer", "3 representations". */
function representationsName(list: readonly StructureRepresentationRef[]): string {
  if (list.length === 1) return `${typeOf(list[0])} of ${componentName(list[0].component)}`;
  return `${list.length} representations`;
}

function typeIds(ctx: ApiContext): string {
  return ctx.plugin.representation.structure.registry.types.map((t) => t[0]).join(', ');
}

/** A structure's representations of one target: the structure, its components that the target covers. */
interface CoveredPart {
  readonly structure: StructureRef;
  readonly components: StructureComponentRef[];
  readonly loci?: StructureElement.Loci;
}

interface Covered {
  readonly summary: string;
  readonly parts: CoveredPart[];
  readonly reprs: StructureRepresentationRef[];
}

/** The target resolved per structure, to the components with representations that it covers (targets.md §1–2). */
async function resolveCovered(
  ctx: ApiContext,
  target: unknown,
  structure: string,
  language: 'pymol' | 'vmd' | 'jmol' | 'molscript',
): Promise<Covered | ApiFailure> {
  const resolved = await ctx.targets.resolve(target, {
    ...(structure ? { structure } : {}),
    language,
  });
  if (!resolved.ok) return resolved;
  if (resolved.empty) {
    return ApiResult.fail('not-found', resolved.summary, 'Check the target, or name another structure.');
  }
  const structures = loadedStructures(ctx);
  const parts: CoveredPart[] = [];
  for (const p of resolved.parts) {
    if (!p.selected.atoms) continue;
    const s = structures.find((x) => x.cell.transform.ref === p.structure.ref);
    if (!s) continue;
    const components = s.components.filter((c) => {
      const data = c.cell.obj?.data;
      return (
        !!data &&
        c.representations.length > 0 &&
        !StructureElement.Loci.isEmpty(StructureElement.Loci.remap(p.loci, data))
      );
    });
    if (components.length) parts.push({ structure: s, components, loci: p.loci });
  }
  if (!parts.length) {
    return ApiResult.fail('not-found', `${resolved.summary} No representation shows it.`, 'Show it first.');
  }
  return {
    summary: resolved.summary,
    parts,
    reprs: parts.flatMap((p) => p.components.flatMap((c) => c.representations)),
  };
}

/** Every representation of the named structure (or all). */
async function allCovered(ctx: ApiContext, structure: string | undefined): Promise<Covered | ApiFailure> {
  const targets = someStructures(ctx, structure);
  if (isFailure(targets)) return targets;
  const parts = targets
    .map((s) => ({ structure: s, components: s.components.filter((c) => c.representations.length > 0) }))
    .filter((p) => p.components.length);
  if (!parts.length)
    return ApiResult.fail('not-found', 'There are no representations to change.', 'Show something first.');
  return {
    summary: 'everything',
    parts,
    reprs: parts.flatMap((p) => p.components.flatMap((c) => c.representations)),
  };
}

/**
 * What a theme or a size changes, as text: the representations of the whole components a target touches (unlike a
 * layer, which holds just the target). A target that covers its components only in part names them, with a warning.
 */
function touchedText(where: Covered): { text: string; warning?: string } {
  const partial = where.parts.some((p) =>
    p.components.some((c) => {
      const data = c.cell.obj?.data;
      return (
        !!p.loci && !!data && StructureElement.Loci.size(StructureElement.Loci.remap(p.loci, data)) < data.elementCount
      );
    }),
  );
  if (!partial) return { text: where.summary };
  const names = [...new Set(where.parts.flatMap((p) => p.components.map(componentName)))].join(', ');
  return {
    text: `the whole of ${names} (holding ${where.summary})`,
    warning: `This applies to whole components: ${names} changed, not only ${where.summary}.`,
  };
}

// --- layers -----------------------------------------------------------------------------------------------------

interface LayerCell {
  readonly ref: string;
  readonly parent: string;
  readonly params: unknown;
}

/** The layer cells of `kind` under these representations. */
function layerCells(ctx: ApiContext, kind: LayerKind, reprs: readonly StructureRepresentationRef[]): LayerCell[] {
  const out: LayerCell[] = [];
  for (const r of reprs) {
    for (const cell of styleLayerCellsOf(ctx.plugin, kind, r.cell.transform.ref)) {
      out.push({ ref: cell.transform.ref, parent: r.cell.transform.ref, params: cell.transform.params });
    }
  }
  return out;
}

/**
 * What a change would touch, as text: the representations' params and their layer cells' (a state update of equal
 * params still makes a new version, so versions say nothing). Equal before and after: nothing changed.
 */
function signature(ctx: ApiContext, reprs: readonly StructureRepresentationRef[]): string {
  const state = ctx.plugin.state.data;
  const parts: string[] = [];
  for (const r of reprs) {
    parts.push(`${r.cell.transform.ref}=${JSON.stringify(state.cells.get(r.cell.transform.ref)?.transform.params)}`);
  }
  for (const kind of Object.keys(Layers) as LayerKind[]) {
    for (const c of layerCells(ctx, kind, reprs)) parts.push(`${c.ref}=${JSON.stringify(c.params)}`);
  }
  return parts.join('\n');
}

/**
 * Applies a layer per structure; the loci of a part are remapped onto the structure each representation is built on.
 * `skipEmpty`: setting the neutral value (0) where no layer exists changes nothing, so those components are left out.
 */
async function applyLayer(
  ctx: ApiContext,
  where: Covered,
  kind: LayerKind,
  did: string,
  set: (part: CoveredPart, getLoci: (s: Structure) => Promise<StructureElement.Loci>) => Promise<void>,
  skipEmpty = false,
) {
  const parts = skipEmpty
    ? where.parts
        .map((p) => ({
          ...p,
          components: p.components.filter((c) => layerCells(ctx, kind, c.representations).length > 0),
        }))
        .filter((p) => p.components.length)
    : where.parts;
  const reprs = parts.flatMap((p) => p.components.flatMap((c) => c.representations));
  if (!reprs.length) return ApiResult.already(`${where.summary} has no ${Layers[kind].name} to clear.`);

  const before = signature(ctx, reprs);
  await ctx.plugin.dataTransaction(async () => {
    for (const part of parts) {
      const loci = part.loci!;
      await set(part, async (s) => StructureElement.Loci.remap(loci, s));
    }
  });
  if (before === signature(ctx, reprs))
    return ApiResult.already(`${where.summary} already has that ${Layers[kind].name}.`);
  return ApiResult.ok(did, { nodes: reprs.map(refOf) });
}

// --- themes and sizes -------------------------------------------------------------------------------------------

/** A color or size theme by name; `value`: the uniform color theme's color. */
type ThemeChange = { color?: string; size?: string; value?: Color };

/** Sets a color or size theme on the representations; `warnings` names what a theme does not apply to. */
async function changeTheme(ctx: ApiContext, where: Covered, change: ThemeChange, warnings: string[]): Promise<void> {
  const themes = ctx.plugin.representation.structure.themes;
  const manager = ctx.plugin.managers.structure.component;
  const name = (change.color ?? change.size)!;
  const registry = change.color ? themes.colorThemeRegistry : themes.sizeThemeRegistry;
  const provider = registry.get(name);
  for (const part of where.parts) {
    const applicable = part.components.filter((c) => {
      const data = c.cell.obj?.data;
      return !!data && provider.isApplicable({ structure: data });
    });
    const skipped = part.components.length - applicable.length;
    if (skipped)
      warnings.push(
        `${provider.label} does not apply to ${skipped} component${skipped === 1 ? '' : 's'} of ${labelOf(ctx, part.structure)}.`,
      );
    const stale = applicable.filter((c) =>
      c.representations.some((r) => {
        const p = r.cell.transform.params as RepresentationParams;
        if ((change.color ? p.colorTheme.name : p.sizeTheme.name) !== name) return true;
        return change.value !== undefined && p.colorTheme.params?.value !== change.value;
      }),
    );
    const { value, ...theme } = change;
    if (stale.length)
      await manager.updateRepresentationsTheme(
        stale,
        (value === undefined ? theme : { ...theme, colorParams: { value } }) as never,
      );
  }
}

async function applyTheme(ctx: ApiContext, where: Covered, change: ThemeChange) {
  const before = signature(ctx, where.reprs);
  const warnings: string[] = [];
  await ctx.plugin.dataTransaction(() => changeTheme(ctx, where, change, warnings));
  const name = change.color ?? change.size;
  const touched = touchedText(where);
  if (before === signature(ctx, where.reprs)) {
    return warnings.length
      ? ApiResult.fail('blocked', warnings.join(' '))
      : ApiResult.already(`${touched.text} already uses the ${name} ${change.color ? 'color' : 'size'} theme.`);
  }
  if (touched.warning) warnings.push(touched.warning);
  return ApiResult.ok(
    `Colored ${touched.text} by ${name}.`,
    { nodes: where.reprs.map(refOf) },
    warnings.length ? { warnings } : {},
  );
}

/** The size factor: a multiple of each representation type's own default `sizeFactor`; types without one are skipped. */
async function changeSizeFactor(ctx: ApiContext, where: Covered, factor: number, warnings: string[]) {
  const themes = ctx.plugin.representation.structure.themes;
  const registry = ctx.plugin.representation.structure.registry;
  const update = ctx.plugin.state.data.build();
  let skipped = 0;
  let any = false;
  for (const part of where.parts) {
    for (const c of part.components) {
      const data = c.cell.obj?.data;
      for (const r of c.representations) {
        const old = r.cell.transform.params as RepresentationParams;
        if (!data || !('sizeFactor' in old.type.params)) {
          skipped++;
          continue;
        }
        const defaults = PD.getDefaultValues(registry.get(old.type.name).getParams(themes, data)) as Record<
          string,
          any
        >;
        const want = round((defaults.sizeFactor as number) * factor);
        if (Math.abs(old.type.params.sizeFactor - want) < 1e-6) continue;
        any = true;
        update
          .to(r.cell.transform.ref)
          .update({ ...old, type: { name: old.type.name, params: { ...old.type.params, sizeFactor: want } } });
      }
    }
  }
  if (skipped) warnings.push(`${skipped} representation${skipped === 1 ? ' has' : 's have'} no size factor.`);
  if (any) await update.commit();
}

const round = (x: number) => Math.round(x * 1e6) / 1e6;

// --- params -----------------------------------------------------------------------------------------------------

/** The type params the component manager adds a representation with (its options: hydrogens, quality, …). */
function managerTypeParams(ctx: ApiContext, type: string): Record<string, any> {
  const {
    hydrogens,
    visualQuality: quality,
    ignoreLight,
    materialStyle: material,
    clipObjects: clip,
    solidSurface,
    interior,
    animation,
  } = ctx.plugin.managers.structure.component.state.options;
  return {
    ignoreHydrogens: hydrogens !== 'all',
    ignoreHydrogensVariant: hydrogens === 'only-polar' ? 'non-polar' : 'all',
    quality,
    ignoreLight,
    material,
    clip,
    interior,
    animation,
    ...(isSurfaceRepresentationType(type) && { solidInterior: solidSurface }),
  };
}
