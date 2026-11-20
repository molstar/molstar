/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `legend` behavior (_spec/README.md "Legends", "Commands"): editing a theme is one patch, previewed (a scheme
// hovered in a menu), and the legend's gestures as commands: a layer's channel (`layer.channel`, several layers in one
// patch), a scale node's params (`scale.set`: every layer naming it follows), one value's color (`theme.value-color`),
// a legend entry's items selected or highlighted (`legend.select`, `legend.hover`), and the reader's color-blind
// substitution (`theme.cvd-safe`, session state, never the spec); a value's items hidden (`layer.filter`: a patch of the
// layer's `data`, a query that leaves them out; a theme never hides).

/**
 * A layer's data split: what the track is (the annotation and a leading run of item forms: its base) and the values of
 * `field` the user's filter, a trailing `(except …)`, leaves out (_spec/README.md "Themes", Filters; core's `splitDataForms`
 * is the one parser, so what the viewer calls the base is what the layers read as it).
 */
export function dataFilter(data: string, field: string): { readonly base: string; readonly hidden: readonly string[] } {
  const p = parseAddress(data).pattern;
  const segs = p?.segments ?? [];
  const ids: string[] = [];
  for (const s of segs) {
    if (typeof s !== 'string') break;
    ids.push(s);
  }
  const split = splitDataForms(segs.slice(ids.length));
  if (!p || !ids.length || split.unread) return { base: data, hidden: [] };
  const base = printAddress([...ids, ...split.lead]);
  const of = split.filter?.params?.of;
  const v = isNode(of) && of.kind === 'item' ? of.params?.[field] : undefined;
  const hidden = v === undefined ? [] : (Array.isArray(v) ? v : [v]).map(String);
  return { base, hidden };
}

/** The address of a value's items on a layer's data: the base's forms and `(item :<field> v)`, one query. */
export function itemsAddress(data: string, field: string, value: string | readonly string[]): string | undefined {
  const p = parseAddress(data).pattern;
  const segs = p?.segments ?? [];
  const ids: string[] = [];
  for (const s of segs) {
    if (typeof s !== 'string') break;
    ids.push(s);
  }
  const split = splitDataForms(segs.slice(ids.length));
  if (!p || !ids.length || split.unread) return undefined;
  const one: Node = { kind: 'item', params: { [field]: typeof value === 'string' ? value : [...value] } };
  // item forms alone: one item holding every param (they are and-ed); an `intersect` or `merge` in the base: an `intersect` of all
  const plain = split.lead.every((f) => f.kind === 'item');
  const form: Node = plain
    ? { kind: 'item', params: Object.assign({}, ...split.lead.map((f) => f.params ?? {}), one.params) }
    : { kind: 'intersect', params: { of: [...split.lead, one] } };
  return printAddress([...ids, form]);
}

import type { Op } from '@molstar/core/repr/patch';
import { isNode, type Node, parseAddress, printAddress, resolveAddress, type Value } from '@molstar/core/repr/tree';
import { splitDataForms } from '@molstar/seqx-core/view';
import {
  builtInSchemes,
  colorHex,
  namedList,
  type SchemeEntry,
  type SchemeRegistry,
  TokenLists,
} from '@molstar/seqx-core/theme/index';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CommandArgs, CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { AddressModel, EditOutcome, SpecModel } from '@molstar/seqx-viewer/models';
import type { EditOrigin } from '@molstar/seqx-viewer/ports';
import type { LegendModel } from '@molstar/seqx-viewer/legends';

export interface ThemingContext {
  /** The layer (`layer.channel`, `theme.value-color`, `legend.*`) or the scale (`scale.set`). */
  readonly address?: string;
  /** `layer.channel`, `layer.filter` cleared (no `entry`): several layers at once, one patch. */
  readonly addresses?: readonly string[];
  /** `layer.channel`: the channel (`color` when absent). */
  readonly channel?: string;
  /** `layer.channel`: its form or value; absent: unset (the layer type's default again), or `scheme`'s. */
  readonly value?: Value;
  /**
   * `layer.channel`'s select param (one layer's color: a fitting scheme by name), as a header's submenu item passes it:
   * read when `value` is absent.
   */
  readonly scheme?: string;
  /** `scale.set`: the params set (`null` at a domain's end reopens it; `undefined` leaves a param). */
  readonly scale?: Readonly<Record<string, Value | undefined>>;
  /** `theme.value-color`, `legend.*`: the value (its text). */
  readonly entry?: string;
  /** `theme.value-color`: its color; absent: the value's entry removed. */
  readonly color?: string;
  /** `legend.select`: toggle the value's items in or out of the selection (Cmd/Ctrl). */
  readonly toggle?: boolean;
  /** `theme.cvd-safe`: on or off (toggles when absent). */
  readonly on?: boolean;
  readonly origin?: EditOrigin;
}

export const ThemingCommandIds = {
  channel: 'layer.channel',
  scale: 'scale.set',
  valueColor: 'theme.value-color',
  select: 'legend.select',
  hover: 'legend.hover',
  cvdSafe: 'theme.cvd-safe',
  filter: 'layer.filter',
} as const;

export interface ThemingBehavior {
  readonly commands: readonly CommandDefinition<ThemingContext>[];
  /** The schemes that fit a layer (what its menu and its legend's title offer), the reader's substitution said. */
  fittingSchemes(layer: string): SchemeEntry[];
  /** A layer's channel previewed (a scheme hovered in a menu): a gesture's preview, ended by `layer.channel` or Esc. */
  preview(layer: string, channel: string, value: Value | undefined): boolean;
  cancelPreview(): void;
  dispose(): void;
}

const ids = (address: string) => address.split('/').filter(Boolean);

/** The field a layer's color reads (its legend's): a `field` form's name, a scheme's `field`, else none. */
export function legendField(layer: Node): string | undefined {
  const c = layer.params?.color;
  if (!isNode(c)) return undefined;
  if (c.kind === 'field') return typeof c.params?.name === 'string' ? c.params.name : 'value';
  if (c.kind === 'scheme' && typeof c.params?.field === 'string') return c.params.field;
  return undefined;
}

export function themingBehavior(
  spec: SpecModel,
  addresses: AddressModel,
  legends: LegendModel,
  schemes: SchemeRegistry = builtInSchemes(),
): ThemingBehavior {
  let disposed = false;
  let previewing = false;
  /** The values toggled into the selection, by layer: the selection is their `or`. */
  const toggled = new Map<string, string[]>();
  const node = (address: string | undefined, kind: string): Node | undefined => {
    if (address === undefined) return undefined;
    const n = resolveAddress(spec.root.value, ids(address));
    return n?.kind === kind ? n : undefined;
  };
  const said =
    (done: string) =>
    (o: EditOutcome): CommandResult =>
      o.ok
        ? { message: o.changed ? done : 'Already so', undoable: o.changed }
        : { message: o.issues.find((i) => i.severity === 'error')?.message ?? 'Refused' };
  const edit = (label: string, ops: () => Op[], ctx: ThemingContext, done: string) => {
    if (previewing) cancelPreview();
    const o = spec.patch(label, (_root, revision) => ({ base: revision, ops: ops() }), {
      ...(ctx.origin ? { origin: ctx.origin } : {}),
    });
    return o instanceof Promise ? o.then(said(done)) : said(done)(o);
  };
  const layers = (ctx: ThemingContext) => ctx.addresses ?? (ctx.address !== undefined ? [ctx.address] : []);
  const channelOp = (layer: string, channel: string, value: Value | undefined): Op =>
    value === undefined
      ? { op: 'unset', at: ids(layer), key: channel }
      : { op: 'set', at: ids(layer), key: channel, value };
  const live = (f: (ctx: ThemingContext) => true | string) => (ctx: ThemingContext) =>
    disposed ? 'The viewer is retired' : f(ctx);
  /** `layer.channel`'s value: the context's, else its select param's scheme (the args', then the context's). */
  const channelValue = (args: CommandArgs, ctx: ThemingContext): Value | undefined => {
    if (ctx.value !== undefined) return ctx.value;
    const name = typeof args.scheme === 'string' ? args.scheme : ctx.scheme;
    return name === undefined ? undefined : ({ kind: 'scheme', params: { name } } as unknown as Value);
  };
  /**
   * `layer.channel`'s choices (_spec/README.md "Legends", the legend behavior's commands: the command owns them): on one
   * layer's color, its fitting schemes as one select param, `scheme`, its default the scheme the layer has (none: no
   * default, optional). What a header's menu (a submenu), a palette and an agent read. A literal, not `PD.Select`: the
   * module stays out of the bundle.
   */
  const channelParams = (ctx: ThemingContext): PD.Params => {
    const ls = layers(ctx);
    const fit = ls.length === 1 && (ctx.channel ?? 'color') === 'color' ? fittingSchemes(ls[0]) : [];
    if (!fit.length) return {};
    const c = node(ls[0], 'layer')?.params?.color;
    const scheme = {
      type: 'select',
      isOptional: true,
      defaultValue: isNode(c) && c.kind === 'scheme' ? c.params?.name : undefined,
      options: fit.map((s) => [s.name, s.label]),
    } as PD.Select<string>;
    return { scheme };
  };

  /** The address of a value's items on a layer's data: its base's forms and `(item :<field> v)` in one query. */
  const itemsOf = (layer: string, value: string): string | undefined => {
    const n = node(layer, 'layer');
    const data = n?.params?.data;
    const field = n && legendField(n);
    if (typeof data !== 'string' || !field) return undefined;
    return itemsAddress(data, field, value);
  };

  /** A layer's color with one value's color set (or removed): a dictionary entry of its palette, or of its scale's. */
  const valueColorOps = (layer: string, value: string, color: string | undefined): Op[] | string => {
    const n = node(layer, 'layer');
    if (!n) return 'No layer there';
    const c = n.params?.color;
    // a field form's palette, or none (the layer's default over its category, written out); a scheme, rules or any
    // other form is not one value's color to set: refused, never replaced
    if (c !== undefined && !(isNode(c) && c.kind === 'field'))
      return 'The layer’s color is not by a field: a value’s color is set on a field’s palette only';
    const field: Node = isNode(c) ? c : { kind: 'field', params: { name: 'category' } };
    const scale = field.params?.scale;
    if (typeof scale === 'string') {
      const s = node(scale, 'scale');
      if (!s) return 'No scale there';
      const palette = s.params?.palette;
      return [
        {
          op: 'set',
          at: ids(scale),
          key: 'palette',
          value: withEntry(isNode(palette) ? palette : undefined, value, color),
        },
      ];
    }
    const palette = field.params?.palette;
    const next: Node = {
      ...field,
      params: { ...field.params, palette: withEntry(isNode(palette) ? palette : undefined, value, color) },
    };
    return [{ op: 'set', at: ids(layer), key: 'color', value: next as unknown as Value }];
  };

  /** A layer's data with a value's items hidden or shown again (`value` absent: every filter cleared). */
  const filterOps = (layer: string, value: string | undefined): Op[] | string => {
    const n = node(layer, 'layer');
    const data = n?.params?.data;
    if (!n || typeof data !== 'string') return 'No layer with data there';
    const field = legendField(n) ?? 'category';
    const { base, hidden } = dataFilter(data, field);
    if (value === undefined) {
      if (base === data) return 'Not filtered';
      return [{ op: 'set', at: ids(layer), key: 'data', value: base }];
    }
    const next = hidden.includes(value) ? hidden.filter((h) => h !== value) : [...hidden, value];
    const to = next.length
      ? `${base}/${printAddress([{ kind: 'except', params: { of: { kind: 'item', params: { [field]: next.length === 1 ? next[0] : next } } } }]).slice(1)}`
      : base;
    return [{ op: 'set', at: ids(layer), key: 'data', value: to }];
  };

  const cancelPreview = () => {
    if (!previewing) return;
    previewing = false;
    spec.cancelPreview();
  };
  function fittingSchemes(layer: string): SchemeEntry[] {
    const n = node(layer, 'layer');
    if (!n) return [];
    const type = n.params?.type;
    const data = n.params?.data;
    const target = typeof data === 'string' ? node(data, 'sequence') : undefined;
    const alphabet = type === 'letters' ? (target?.params?.alphabet ?? letterAlphabet(spec, layer)) : undefined;
    const subject =
      alphabet === 'protein' ? 'protein' : alphabet ? 'nucleotide' : type === 'values' ? 'numbers' : 'categories';
    return schemes.fitting(subject);
  }
  function previewChannel(layer: string, channel: string, value: Value | undefined): boolean {
    if (disposed || !node(layer, 'layer')) return false;
    const ok = spec.setPreview({ ops: [channelOp(layer, channel, value)] }).ok;
    previewing = ok;
    return ok;
  }

  const commands: CommandDefinition<ThemingContext>[] = [
    {
      id: ThemingCommandIds.channel,
      title: 'Theme a layer’s channel',
      category: 'Edit',
      when: live((ctx) => (layers(ctx).length && layers(ctx).every((a) => node(a, 'layer')) ? true : 'No layer there')),
      params: (ctx) => channelParams(ctx),
      // a choice hovered in a menu: the gesture's preview (ended by the choice's patch, or by Esc)
      preview: (args, ctx) => {
        const ls = layers(ctx);
        return ls.length === 1 && previewChannel(ls[0], ctx.channel ?? 'color', channelValue(args, ctx));
      },
      run: (args, ctx) => {
        const value = channelValue(args, ctx);
        return edit(
          'Theme a layer',
          () => layers(ctx).map((a) => channelOp(a, ctx.channel ?? 'color', value)),
          ctx,
          value === undefined ? 'Default theme' : 'Themed',
        );
      },
    },
    {
      id: ThemingCommandIds.scale,
      title: 'Set a scale',
      category: 'Edit',
      when: live((ctx) => (node(ctx.address, 'scale') ? true : 'No scale there')),
      run: (_, ctx) =>
        edit(
          'Set a scale',
          () =>
            Object.entries(ctx.scale ?? {})
              .filter(([, v]) => v !== undefined)
              .map(([key, value]) => ({ op: 'set', at: ids(ctx.address!), key, value: value as Value }) as Op),
          ctx,
          'Scale set',
        ),
    },
    {
      id: ThemingCommandIds.valueColor,
      title: 'Color one value',
      category: 'Edit',
      when: live((ctx) =>
        ctx.entry === undefined ? 'No value' : node(ctx.address, 'layer') ? true : 'No layer there',
      ),
      run: (_, ctx) => {
        const ops = valueColorOps(ctx.address!, ctx.entry!, ctx.color);
        if (typeof ops === 'string') return { message: ops };
        return edit('Color one value', () => ops, ctx, ctx.color === undefined ? 'Default color' : 'Colored');
      },
    },
    {
      id: ThemingCommandIds.filter,
      title: 'Hide a value’s items',
      category: 'Legend',
      // a value hidden on one layer (`address`); no value: every listed layer's filter cleared in one patch (`addresses`:
      // a track's "n filtered" chip)
      when: live((ctx) =>
        layers(ctx).length && layers(ctx).every((a) => node(a, 'layer')) && (ctx.entry === undefined || !ctx.addresses)
          ? true
          : 'No layer there',
      ),
      run: (_, ctx) => {
        let ops: Op[] | string = 'Not filtered';
        for (const a of layers(ctx)) {
          const o = filterOps(a, ctx.entry);
          if (typeof o !== 'string') ops = typeof ops === 'string' ? o : [...ops, ...o];
          else if (layers(ctx).length === 1) ops = o;
        }
        if (typeof ops === 'string') return { message: ops };
        const all = ops;
        return edit(
          ctx.entry === undefined ? 'Clear a filter' : 'Filter a layer',
          () => all,
          ctx,
          ctx.entry === undefined ? 'Filter cleared' : 'Filtered',
        );
      },
    },
    {
      id: ThemingCommandIds.select,
      title: 'Select a legend entry’s items',
      category: 'Legend',
      when: live((ctx) =>
        ctx.entry !== undefined && ctx.address !== undefined && itemsOf(ctx.address, ctx.entry)
          ? true
          : 'No entry there',
      ),
      run: (_, ctx) => {
        const layer = ctx.address!;
        const prev = toggled.get(layer) ?? [];
        // the grammar's click rule: clicking the only selected value clears the selection
        let values: string[];
        if (ctx.toggle)
          values = prev.includes(ctx.entry!) ? prev.filter((v) => v !== ctx.entry) : [...prev, ctx.entry!];
        else values = prev.length === 1 && prev[0] === ctx.entry ? [] : [ctx.entry!];
        toggled.set(layer, values);
        if (!values.length) {
          addresses.setSelection(undefined);
          return { message: 'Cleared' };
        }
        const n = node(layer, 'layer')!;
        const data = n.params?.data;
        const address =
          typeof data === 'string'
            ? itemsAddress(data, legendField(n)!, values.length === 1 ? values[0] : values)
            : undefined;
        if (!address) return { message: 'No entry there' };
        return { message: addresses.select(address).message };
      },
    },
    {
      id: ThemingCommandIds.hover,
      title: 'Highlight a legend entry’s items',
      category: 'Legend',
      when: live(() => true),
      run: (_, ctx) => {
        if (ctx.entry === undefined || ctx.address === undefined) {
          addresses.clearHighlight('legend');
          return { message: 'Cleared' };
        }
        const a = itemsOf(ctx.address, ctx.entry);
        return { message: a ? addresses.highlightAddress(a, 'legend').message : 'No entry there' };
      },
    },
    {
      id: ThemingCommandIds.cvdSafe,
      title: 'Color-blind safe colors',
      category: 'View',
      when: live(() => true),
      run: (_, ctx) => {
        const on = ctx.on ?? !legends.cvdSafe.value;
        legends.setCvdSafe(on);
        return { message: on ? 'Color-blind safe colors' : 'The figure’s colors' };
      },
    },
  ];

  return {
    commands,
    fittingSchemes,
    preview: previewChannel,
    cancelPreview,
    dispose() {
      disposed = true;
      cancelPreview();
    },
  };
}

/** A letters layer's alphabet: its section's system's. */
function letterAlphabet(spec: SpecModel, layer: string): string | undefined {
  const section = ids(layer).slice(0, -2);
  const s = resolveAddress(spec.root.value, section);
  const of = s?.params?.of;
  if (typeof of !== 'string') return undefined;
  const sys = resolveAddress(spec.root.value, ids(of));
  const a = sys?.params?.alphabet;
  return typeof a === 'string' ? a : undefined;
}

/** A categorical palette with one value's color set or removed; a default palette's other values keep theirs. */
function withEntry(palette: Node | undefined, value: string, color: string | undefined): Value {
  const params = (palette?.kind === 'categorical' ? palette.params : undefined) ?? {};
  const colors = params.colors;
  const dict: Record<string, Value> =
    colors && typeof colors === 'object' && !Array.isArray(colors) && !isNode(colors as Value)
      ? { ...(colors as Record<string, Value>) }
      : {};
  // listed values keep the color their place gave them: written into the dictionary
  const listed = Array.isArray(params.values) ? (params.values as string[]) : [];
  const list = typeof colors === 'string' ? listColors(colors) : Array.isArray(colors) ? (colors as Value[]) : [];
  listed.forEach((v, i) => {
    if (!(v in dict) && list.length) dict[v] = list[i % list.length];
  });
  if (color === undefined) delete dict[value];
  else dict[value] = color;
  // a list or a named list becomes the `others` the dictionary hashes into: every other value keeps its color
  const others = typeof colors === 'string' || Array.isArray(colors) ? colors : (params.others ?? 'categories');
  const { values: _, ...rest } = params as Record<string, Value>;
  return { kind: 'categorical', params: { ...rest, colors: dict, others } } as unknown as Value;
}

/** A named list's colors as a spec writes them: its tokens (`token:cat-0`), else Mol*'s colors as hex. */
function listColors(name: string): Value[] {
  const tokens = TokenLists[name];
  if (tokens) return tokens.map((t) => `token:${t.slice(4)}`);
  return (namedList(name, { theme: 'light', colors: {}, version: 0 }) ?? []).map((c) => colorHex(c).slice(0, 7));
}
