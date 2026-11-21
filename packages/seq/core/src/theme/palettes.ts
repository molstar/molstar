/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color/color';
import {
  ColorListBlues,
  ColorListCategory10,
  ColorListCividis,
  ColorListGreens,
  ColorListMagma,
  ColorListOkabeIto,
  ColorListRedBlue,
  ColorListTableau10,
  ColorListViridis,
} from '@molstar/core/util/color/lists';
import type { ColorList } from '@molstar/core/util/color/color';
import {
  makeValuePaletteFunction,
  type ValuePaletteColor,
  type ValuePaletteFunction,
  type ValuePaletteSpec,
  type ValuePaletteSpecDefaults,
  valuePalettePropsFromSpec,
} from '@molstar/core/util/color/value-palette-eval';
import { fromRgb, resolveColor, type TokenTable, toRgb } from './colors.js';
import { isThemeNode, type PaletteKind, param, type ScaleParams, type ThemeNode } from './forms.js';

// Palettes (_spec/README.md "Themes", Palettes and scales): spelled in this language's words, evaluated by Mol*'s value
// palettes (`@molstar/core/util/color/value-palette-eval`, the evaluation without Mol*'s named lists), so a value has
// one color in MolViewSpec and here. A continuous palette is interpolated in Oklab (fixed) into a lookup table of 256
// colors, rebuilt only when the tokens' version or the domain changes; a categorical value's color is memoized per
// value.

/** The lookup table's size (fixed: indistinguishable steps at any width). */
export const LutSize = 256;

/** The palettes named by tokens: they follow the token table (light, dark, a host's `--seq-*`). */
export const TokenLists: Readonly<Record<string, readonly string[]>> = {
  ramp: ['seq-ramp-0', 'seq-ramp-1', 'seq-ramp-2', 'seq-ramp-3', 'seq-ramp-4'],
  diverging: ['seq-div-neg', 'seq-div-mid', 'seq-div-pos'],
  categories: ['seq-cat-0', 'seq-cat-1', 'seq-cat-2', 'seq-cat-3', 'seq-cat-4', 'seq-cat-5', 'seq-cat-6', 'seq-cat-7'],
};

/**
 * Mol*'s lists this language names (the ones its schemes use), by this language's name: no other list is offered, so a
 * spec reads the same wherever Seq* runs, and a list added to Mol* is not one here until it is named. Each imported by
 * name: the bundle ships these lists and none of Mol*'s others.
 */
export const NamedLists: Readonly<Record<string, ColorList>> = {
  viridis: ColorListViridis,
  cividis: ColorListCividis,
  magma: ColorListMagma,
  blues: ColorListBlues,
  greens: ColorListGreens,
  // the same stops, dark to light (a confident low end dark: AlphaFold DB's PAE plot)
  'greens-r': { ...ColorListGreens, list: [...ColorListGreens.list].reverse() },
  'rd-bu': ColorListRedBlue,
  'okabe-ito': ColorListOkabeIto,
  tableau10: ColorListTableau10,
  category10: ColorListCategory10,
};

/** A named list's colors (RGBA) in a token table; undefined for a name that is no list. */
export function namedList(name: string, tokens: TokenTable): number[] | undefined {
  const t = TokenLists[name];
  if (t) return t.map((n) => tokens.colors[n] ?? 0);
  if (!Object.hasOwn(NamedLists, name)) return undefined;
  return NamedLists[name].list.map((e) => fromRgb(Color.fromColorListEntry(e)));
}

/** Where a palette's colors failed to read: an unknown name or a color that is none. */
export type PaletteReport = (what: 'name' | 'color', value: string) => void;

/** Colors of a palette's `colors`: a named list or a list of colors, RGBA (an unreadable color is the missing token). */
function colorList(v: unknown, tokens: TokenTable, report: PaletteReport): number[] {
  if (typeof v === 'string') {
    const list = namedList(v, tokens);
    if (!list) report('name', v);
    return list ?? [];
  }
  if (!Array.isArray(v)) return [];
  return v.map((c) => color(c, tokens, report));
}

function color(c: unknown, tokens: TokenTable, report: PaletteReport): number {
  const out = resolveColor(c, tokens);
  if (out === undefined) report('color', String(c));
  return out ?? tokens.colors['seq-missing'] ?? 0;
}

/** Mol*'s spec of colors, RGB. */
const mol = (c: number): ValuePaletteColor => toRgb(c);

const Defaults: ValuePaletteSpecDefaults = {
  categorical: {
    kind: 'categorical',
    colors: [],
    repeat_color_list: false,
    sort: 'none',
    sort_direction: 'ascending',
    case_insensitive: false,
    missing_color: null,
  },
  discrete: { kind: 'discrete', colors: [], reverse: false, mode: 'absolute', value_domain: [null, null] },
  continuous: {
    kind: 'continuous',
    colors: [],
    reverse: false,
    mode: 'normalized',
    value_domain: [0, 1],
    underflow_color: null,
    overflow_color: null,
  },
};

/**
 * A palette ready to color: `color(v)` by a value (a string or a number; undefined: missing, the caller's missing
 * token), `num(x)` a number's color without a string. Numeric palettes carry their resolved domain.
 */
export interface Palette {
  readonly kind: PaletteKind;
  /** A value's color, RGBA; undefined: missing (NaN, absent, outside a log scale, a `direct` value that does not read). */
  color(v: string | number | undefined): number | undefined;
  /** A number's color (continuous and binned: the lookup table or the bins; others: as `color`). */
  num(x: number): number | undefined;
  /** Continuous: the stops (RGBA) at 0..1 of the domain; binned: one color a bin. */
  readonly colors: readonly number[];
  /** Continuous: the stops' positions, 0..1. */
  readonly at?: readonly number[];
  /** Binned: the n + 1 edges; continuous: the domain's ends. */
  readonly edges?: readonly number[];
  readonly domain?: readonly [number, number];
  /** Categorical: the listed values in order and their colors. */
  readonly values?: readonly string[];
  /** Categorical: unlisted values are hashed (into `others`). */
  readonly hashed?: boolean;
  readonly scale: 'linear' | 'log' | 'sqrt';
  readonly mid?: number;
  /** A log scale's non-positive values drawn at the low end (else missing). */
  readonly clamp?: boolean;
  /** The colors given below and above the domain (else its ends'). */
  readonly under?: number;
  readonly over?: number;
}

/** A value through the scale (log10, sqrt), NaN outside its domain. */
export function scaled(type: string | undefined, x: number): number {
  if (type === 'log') return x > 0 ? Math.log10(x) : Number.NaN;
  if (type === 'sqrt') return x >= 0 ? Math.sqrt(x) : Number.NaN;
  return x;
}

/** A value's place in the domain, 0..1 through the scale (and the center at 0.5 with `mid`); NaN for missing. */
export function normalize(x: number, scale: string, lo: number, hi: number, mid: number | undefined): number {
  const v = scaled(scale, x);
  if (Number.isNaN(v)) return Number.NaN;
  const a = scaled(scale, lo);
  const b = scaled(scale, hi);
  if (mid !== undefined) {
    const m = scaled(scale, mid);
    if (v < m) return m === a ? 0 : (0.5 * (v - a)) / (m - a);
    return b === m ? 1 : 0.5 + (0.5 * (v - m)) / (b - m);
  }
  return b === a ? 0.5 : (v - a) / (b - a);
}

export interface PaletteOptions {
  readonly tokens: TokenTable;
  /** The resolved domain (an explicit end wins over the data's). */
  readonly domain?: readonly [number, number];
  readonly scale?: ScaleParams;
  readonly report?: PaletteReport;
}

/**
 * The default palette of a field by its column's declared type: categories for text, the ramp for numbers, the
 * diverging tokens for numbers about a center (`mid`).
 */
export function defaultPalette(numeric: boolean, mid?: number): ThemeNode {
  if (!numeric) return { kind: 'categorical', params: { colors: 'categories' } };
  return { kind: 'continuous', params: { colors: mid !== undefined ? 'diverging' : 'ramp' } };
}

/** Compiles a palette node (one of `PaletteKinds`) to a `Palette`. */
export function compilePalette(node: ThemeNode, options: PaletteOptions): Palette {
  const { tokens } = options;
  const report = options.report ?? (() => {});
  const scale = options.scale ?? {};
  const type = scale.type ?? 'linear';
  const mid = (param(node, 'mid') as number | undefined) ?? scale.mid;
  // a value at or below zero on a log scale: missing, or the domain's low end with `non_positive: clamp`
  const clamp = type === 'log' && scale.non_positive === 'clamp';
  if (node.kind === 'direct') {
    return {
      kind: 'direct',
      colors: [],
      scale: 'linear',
      color: (v) => (typeof v === 'string' ? resolveColor(v, tokens) : undefined),
      num: () => undefined,
    };
  }
  if (node.kind === 'categorical') return categorical(node, tokens, report);
  const domain = options.domain ?? [0, 1];
  const under = param(node, 'under');
  const over = param(node, 'over');
  if (node.kind === 'binned') {
    let colors = colorList(param(node, 'colors') ?? 'ramp', tokens, report);
    let edges = param(node, 'edges') as number[] | undefined;
    if (!Array.isArray(edges)) {
      const n = (param(node, 'bins') as number | undefined) ?? colors.length;
      const a = scaled(type, domain[0]);
      const b = scaled(type, domain[1]);
      edges = [];
      for (let i = 0; i <= n; i++) {
        const s = a + ((b - a) * i) / n;
        edges.push(type === 'log' ? 10 ** s : type === 'sqrt' ? s * s : s);
      }
      if (colors.length !== n)
        colors = Array.from(
          { length: n },
          (_, i) => colors[Math.round((i * (colors.length - 1)) / Math.max(1, n - 1))] ?? 0,
        );
    }
    const ranges: [ValuePaletteColor, number, number][] = [];
    for (let i = 0; i + 1 < edges.length && i < colors.length; i++)
      ranges.push([mol(colors[i]), edges[i], edges[i + 1]]);
    const f = makeValuePaletteFunction(
      valuePalettePropsFromSpec({ kind: 'discrete', mode: 'absolute', colors: ranges, scale: type }, Defaults),
      [],
    );
    const underC = under !== undefined ? color(under, tokens, report) : colors[0];
    const overC = over !== undefined ? color(over, tokens, report) : colors[colors.length - 1];
    const lo = edges[0];
    const hi = edges[edges.length - 1];
    const o = { value: 0 };
    const num = (x: number) => {
      if (Number.isNaN(x)) return undefined;
      if (type === 'log' && x <= 0) {
        if (!clamp) return undefined;
        x = lo;
      }
      if (x < lo) return underC;
      if (x > hi) return overC;
      o.value = x;
      const c = f.at!(o);
      return c === undefined ? undefined : fromRgb(c);
    };
    return {
      kind: 'binned',
      colors,
      edges,
      domain: [lo, hi],
      scale: type,
      ...(clamp ? { clamp } : {}),
      ...(under !== undefined ? { under: underC } : {}),
      ...(over !== undefined ? { over: overC } : {}),
      color: (v) => num(typeof v === 'number' ? v : v === undefined ? Number.NaN : Number.parseFloat(v)),
      num,
    };
  }
  // continuous: the stops (a named list spaced evenly, or `[color, at]`; `:stops absolute` reads `at` as values)
  const raw = param(node, 'colors') ?? 'ramp';
  let colors: number[];
  let at: number[];
  let lo = domain[0];
  let hi = domain[1];
  if (Array.isArray(raw) && raw.length && raw.every((s) => Array.isArray(s) && s.length === 2)) {
    colors = raw.map((s) => color(s[0], tokens, report));
    at = raw.map((s) => Number(s[1]));
    if (param(node, 'stops') === 'absolute') {
      if (!options.domain) {
        lo = Math.min(...at);
        hi = Math.max(...at);
      }
      at = at.map((x) => normalize(x, type, lo, hi, undefined));
    }
  } else {
    colors = colorList(raw, tokens, report);
    at = colors.map((_, i) => (colors.length > 1 ? i / (colors.length - 1) : 0));
  }
  const lut = continuousLut(colors, at);
  const underC = under !== undefined ? color(under, tokens, report) : lut[0];
  const overC = over !== undefined ? color(over, tokens, report) : lut[LutSize - 1];
  const num = (x: number) => {
    if (type === 'log' && x <= 0) {
      if (!clamp) return undefined;
      x = lo;
    }
    // `normalize` inlined for the linear scale without a center: the common case on a frame
    const t =
      type === 'linear' && mid === undefined
        ? hi === lo
          ? 0.5
          : (x - lo) / (hi - lo)
        : normalize(x, type, lo, hi, mid);
    if (Number.isNaN(t)) return undefined;
    if (t < 0) return underC;
    if (t > 1) return overC;
    return lut[Math.round(t * (LutSize - 1))];
  };
  return {
    kind: 'continuous',
    colors,
    at,
    domain: [lo, hi],
    edges: [lo, hi],
    scale: type,
    ...(mid !== undefined ? { mid } : {}),
    ...(clamp ? { clamp } : {}),
    ...(under !== undefined ? { under: underC } : {}),
    ...(over !== undefined ? { over: overC } : {}),
    color: (v) => num(typeof v === 'number' ? v : v === undefined ? Number.NaN : Number.parseFloat(v)),
    num,
  };
}

/** 256 colors from stops: interpolated in Oklab (each segment densified), evaluated by Mol*'s continuous palette. */
export function continuousLut(colors: readonly number[], at: readonly number[]): Uint32Array {
  const lut = new Uint32Array(LutSize);
  if (!colors.length) return lut;
  const order = colors.map((_, i) => i).sort((a, b) => at[a] - at[b]);
  const stops: [ValuePaletteColor, number][] = [];
  const k = 16;
  for (let s = 0; s < order.length; s++) {
    const i = order[s];
    const j = order[s + 1];
    stops.push([mol(colors[i]), at[i]]);
    if (j === undefined || at[j] === at[i]) continue;
    for (let q = 1; q < k; q++)
      stops.push([
        Color.interpolateOklab(toRgb(colors[i]), toRgb(colors[j]), q / k),
        at[i] + ((at[j] - at[i]) * q) / k,
      ]);
  }
  const f = makeValuePaletteFunction(
    valuePalettePropsFromSpec(
      { kind: 'continuous', colors: stops, underflow_color: 'auto', overflow_color: 'auto' },
      Defaults,
    ),
    [],
  );
  for (let i = 0; i < LutSize; i++) lut[i] = fromRgb(f(i / (LutSize - 1)) ?? 0);
  return lut;
}

function categorical(node: ThemeNode, tokens: TokenTable, report: PaletteReport): Palette {
  const raw = param(node, 'colors') ?? 'categories';
  const values = (param(node, 'values') as string[] | undefined)?.map(String);
  const hashed = (param(node, 'unlisted') ?? 'hash') === 'hash';
  const dict = typeof raw === 'object' && raw !== null && !Array.isArray(raw) && !isThemeNode(raw);
  let spec: ValuePaletteSpec;
  let colors: number[];
  if (dict) {
    const entries = Object.entries(raw as Record<string, unknown>).map(
      ([k, c]) => [k, color(c, tokens, report)] as const,
    );
    colors = entries.map((e) => e[1]);
    const others = colorList(param(node, 'others') ?? 'categories', tokens, report);
    spec = {
      kind: 'categorical',
      colors: Object.fromEntries(entries.map(([k, c]) => [k, mol(c)])),
      unlisted: hashed ? 'hash' : 'missing',
      others: others.map(mol),
    };
  } else {
    colors = colorList(raw, tokens, report);
    const others = param(node, 'others');
    spec = {
      kind: 'categorical',
      colors: colors.map(mol),
      repeat_color_list: param(node, 'repeat') !== false,
      unlisted: hashed ? 'hash' : 'missing',
      ...(others !== undefined ? { others: colorList(others, tokens, report).map(mol) } : {}),
    };
  }
  const f: ValuePaletteFunction = makeValuePaletteFunction(valuePalettePropsFromSpec(spec, Defaults), values ?? []);
  const memo = new Map<string | number, number | undefined>();
  const colorOf = (v: string | number | undefined) => {
    if (v === undefined || (typeof v === 'number' && Number.isNaN(v))) return undefined;
    let c = memo.get(v);
    if (c === undefined && !memo.has(v)) {
      const m = f(v);
      c = m === undefined ? undefined : fromRgb(m);
      if (memo.size < 4096) memo.set(v, c);
    }
    return c;
  };
  return {
    kind: 'categorical',
    colors,
    values: dict ? Object.keys(raw as object) : (values ?? []),
    hashed,
    scale: 'linear',
    color: colorOf,
    num: colorOf,
  };
}

// --- MolViewSpec's palettes: one implementation, two spellings (a converter both ways).

/** A palette node as MolViewSpec's palette (colors resolved with the table; RGB hex strings). */
export function toMvsPalette(
  node: ThemeNode,
  tokens: TokenTable,
  scale: ScaleParams = {},
): ValuePaletteSpec | undefined {
  const hex = (c: number) => `#${(c >>> 8).toString(16).padStart(6, '0')}`;
  const report = () => {};
  const p = compilePalette(node, {
    tokens,
    scale,
    ...(scale.domain && scale.domain[0] !== null && scale.domain[1] !== null
      ? { domain: scale.domain as [number, number] }
      : {}),
  });
  const extra = {
    ...(p.scale !== 'linear' ? { scale: p.scale } : {}),
    ...(p.mid !== undefined ? { mid: p.mid } : {}),
    ...(p.clamp ? { non_positive: 'clamp' as const } : {}),
  };
  const flows = {
    ...(p.under !== undefined ? { underflow_color: hex(p.under) } : {}),
    ...(p.over !== undefined ? { overflow_color: hex(p.over) } : {}),
  };
  if (node.kind === 'categorical') {
    const raw = param(node, 'colors');
    if (raw && typeof raw === 'object' && !Array.isArray(raw)) {
      const dict: Record<string, string> = {};
      for (const [k, c] of Object.entries(raw)) dict[k] = hex(color(c, tokens, report));
      return { kind: 'categorical', colors: dict, unlisted: p.hashed ? 'hash' : 'missing' };
    }
    return { kind: 'categorical', colors: p.colors.map(hex), unlisted: p.hashed ? 'hash' : 'missing' };
  }
  if (node.kind === 'binned') {
    const e = p.edges ?? [];
    return {
      kind: 'discrete',
      mode: 'absolute',
      colors: p.colors.map((c, i) => [hex(c), e[i], e[i + 1]] as [string, number, number]),
      ...extra,
    };
  }
  if (node.kind === 'continuous') {
    const at = p.at ?? [];
    return {
      kind: 'continuous',
      colors: p.colors.map((c, i) => [hex(c), at[i]] as [string, number]),
      value_domain: [p.domain?.[0] ?? null, p.domain?.[1] ?? null],
      ...extra,
      ...flows,
    };
  }
  return undefined;
}

/**
 * The scale a MolViewSpec palette carries (its `scale`, `value_domain`, `mid`, `non_positive`), as a field form's
 * scale params: with `fromMvsPalette`'s node, the palette converted whole.
 */
export function fromMvsScale(spec: ValuePaletteSpec): ScaleParams {
  if (spec.kind === 'categorical') return {};
  const d = spec.kind === 'continuous' ? spec.value_domain : undefined;
  return {
    ...(spec.scale && spec.scale !== 'linear' ? { type: spec.scale } : {}),
    ...(d && (d[0] !== null || d[1] !== null) ? { domain: [d[0] ?? null, d[1] ?? null] as const } : {}),
    ...(spec.mid !== undefined && spec.mid !== null ? { mid: spec.mid } : {}),
    ...(spec.non_positive ? { non_positive: spec.non_positive } : {}),
  };
}

/** MolViewSpec's palette as a palette node of this language (its scale: `fromMvsScale`). */
export function fromMvsPalette(spec: ValuePaletteSpec): ThemeNode {
  const params: Record<string, unknown> = {};
  if (spec.kind === 'categorical') {
    if (spec.colors !== undefined) params.colors = spec.colors;
    if (spec.repeat_color_list !== undefined) params.repeat = spec.repeat_color_list;
    params.unlisted = spec.unlisted ?? 'missing';
    if (spec.others !== undefined) params.others = spec.others;
    return { kind: 'categorical', params };
  }
  if (spec.mid !== undefined && spec.mid !== null) params.mid = spec.mid;
  if (spec.kind === 'discrete') {
    const c = spec.colors;
    if (Array.isArray(c) && c.every((x) => Array.isArray(x) && x.length === 3)) {
      const ranges = c as readonly (readonly [ValuePaletteColor | null, number | null, number | null])[];
      params.edges = [ranges[0][1], ...ranges.map((r) => r[2])];
      params.colors = ranges.map((r) => r[0]);
    } else params.colors = c;
    return { kind: 'binned', params };
  }
  const c = spec.colors;
  params.colors = c;
  if (spec.mode === 'absolute') params.stops = 'absolute';
  if (spec.underflow_color && spec.underflow_color !== 'auto') params.under = spec.underflow_color;
  if (spec.overflow_color && spec.overflow_color !== 'auto') params.over = spec.overflow_color;
  return { kind: 'continuous', params };
}
