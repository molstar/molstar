/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `legend` behavior (_spec/README.md "Legends", "Commands"): editing a theme is one patch, previewed (a scheme
// hovered in a menu), and the legend's gestures as commands: a layer's channel (`layer.channel`, several layers in one
// patch), a scale node's params (`scale.set`: every layer naming it follows), one value's color (`theme.value-color`),
// a legend entry's items selected or highlighted (`legend.select`, `legend.hover`), and the reader's color-blind
// substitution (`theme.cvd-safe`, session state, never the spec).

import type { Op } from '@molstar/core/repr/patch';
import { isNode, type Node, resolveAddress, type Value } from '@molstar/core/repr/tree';
import {
  builtInSchemes,
  colorHex,
  namedList,
  type SchemeEntry,
  type SchemeRegistry,
  TokenLists,
} from '@molstar/seqx-core/theme/index';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { AddressModel, EditOutcome, SpecModel } from '@molstar/seqx-viewer/models';
import type { EditOrigin } from '@molstar/seqx-viewer/ports';
import type { LegendModel } from '@molstar/seqx-viewer/legends';

export interface ThemingContext {
  /** The layer (`layer.channel`, `theme.value-color`, `legend.*`) or the scale (`scale.set`). */
  readonly address?: string;
  /** `layer.channel`: several layers at once, one patch. */
  readonly addresses?: readonly string[];
  /** `layer.channel`: the channel (`color` when absent). */
  readonly channel?: string;
  /** `layer.channel`: its form or value; absent: unset (the layer type's default again). */
  readonly value?: Value;
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
const quote = (v: string) => (/^[A-Za-z_][A-Za-z0-9_.-]*$/.test(v) ? v : JSON.stringify(v));

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

  /** The address of a value's items on a layer's data: `<data>/(item :<field> v)`. */
  const itemsOf = (layer: string, value: string): string | undefined => {
    const n = node(layer, 'layer');
    const data = n?.params?.data;
    const field = n && legendField(n);
    if (typeof data !== 'string' || !field) return undefined;
    return `${data}/(item :${field} ${quote(value)})`;
  };

  /** A layer's color with one value's color set (or removed): a dictionary entry of its palette, or of its scale's. */
  const valueColorOps = (layer: string, value: string, color: string | undefined): Op[] | string => {
    const n = node(layer, 'layer');
    if (!n) return 'No layer there';
    const c = n.params?.color;
    // the default: the layer's theme written out over the field it reads by default (its category), others kept
    const field: Node = isNode(c) && c.kind === 'field' ? c : { kind: 'field', params: { name: 'category' } };
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

  const cancelPreview = () => {
    if (!previewing) return;
    previewing = false;
    spec.cancelPreview();
  };

  const commands: CommandDefinition<ThemingContext>[] = [
    {
      id: ThemingCommandIds.channel,
      title: 'Theme a layer’s channel',
      category: 'Edit',
      when: live((ctx) => (layers(ctx).length && layers(ctx).every((a) => node(a, 'layer')) ? true : 'No layer there')),
      run: (_, ctx) =>
        edit(
          'Theme a layer',
          () => layers(ctx).map((a) => channelOp(a, ctx.channel ?? 'color', ctx.value)),
          ctx,
          ctx.value === undefined ? 'Default theme' : 'Themed',
        ),
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
        const one = itemsOf(layer, values[0])!;
        const address =
          values.length === 1
            ? one
            : `${one.slice(0, one.lastIndexOf('/('))}/(or :of [${values.map((v) => `(item :${legendField(node(layer, 'layer')!)} ${quote(v)})`).join(' ')}])`;
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
    fittingSchemes(layer) {
      const n = node(layer, 'layer');
      if (!n) return [];
      const type = n.params?.type;
      const data = n.params?.data;
      const target = typeof data === 'string' ? node(data, 'sequence') : undefined;
      const alphabet = type === 'letters' ? (target?.params?.alphabet ?? letterAlphabet(spec, layer)) : undefined;
      const subject =
        alphabet === 'protein' ? 'protein' : alphabet ? 'nucleotide' : type === 'values' ? 'numbers' : 'categories';
      return schemes.fitting(subject);
    },
    preview(layer, channel, value) {
      if (disposed || !node(layer, 'layer')) return false;
      const ok = spec.setPreview({ ops: [channelOp(layer, channel, value)] }).ok;
      previewing = ok;
      return ok;
    },
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
