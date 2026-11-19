/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Oklab } from '@molstar/core/util/color/spaces/oklab';
import { contrastRgba, fromRgb, resolveColor, type TokenTable, toRgb } from './colors.js';
import {
  type ChannelName,
  ChannelWords,
  isThemeNode,
  param,
  type ScaleParams,
  type ThemeBody,
  ThemeChannels,
  type ThemeNode,
} from './forms.js';
import { compilePalette, defaultPalette, normalize, type Palette } from './palettes.js';
import { darkVariant, type SchemeEntry, type SchemeRegistry, type SchemeSubject } from './schemes.js';

// Compiling a form (_spec/README.md "Themes", Performance): once per form, token table and resolved domains, to a
// channel whose `at` reads a datum and returns the mark's value; a uniform form is a constant a layer reads once.

/** What a channel reads of one mark. A layer keeps one and sets its fields per mark (nothing allocated). */
export interface ThemeDatum {
  /** A field's value: text, a number (NaN missing), or undefined (no such field or no value). */
  field(name: string): string | number | undefined;
  /** A values layer's value (NaN missing). */
  value: number;
  /** A letters layer's letter (one character). */
  letter: string;
  /** The mark's extent on its data's system, `[start, end)`: what `in` and `within` read. */
  start: number;
  end: number;
}

/** The codes compiling and evaluating say (cataloged with the render codes). */
export type ThemeIssueCode =
  | 'seq.render.theme.channel'
  | 'seq.render.theme.scheme'
  | 'seq.render.theme.measure'
  | 'seq.render.theme.color'
  | 'seq.render.scale.log'
  | 'seq.render.line-style'
  | 'seq.render.color';

/** A measure a clause reads (a registry entry: the rows layer's column measures register here). */
export interface MeasureEntry {
  readonly name: string;
  read(d: ThemeDatum, params: Readonly<Record<string, unknown>>): number | string | undefined;
}

export interface ThemeEnv {
  readonly tokens: TokenTable;
  readonly schemes: SchemeRegistry;
  readonly measures?: ReadonlyMap<string, MeasureEntry>;
  /** A scale node's params by its address (`/views/conservation`). */
  scale?(address: string): ScaleParams | undefined;
  /** The held domain of an open scale by its key (unioned over the layers naming it at settle). */
  domain?(key: string): readonly [number, number] | undefined;
  /** Whether a field holds numbers (its column's declared type, never its values); undefined: text. */
  numeric?(field: string): boolean | undefined;
  /** A location's ranges on the layer's data system, `[start, end, ...]`, for `in` and `within`. */
  locate?(address: string): readonly number[] | undefined;
  /** The layer's alphabet (a letters layer: which schemes fit). */
  readonly alphabet?: 'protein' | 'dna' | 'rna';
  /** The reader's color-blind substitution (`theme.cvdSafe`). */
  readonly cvdSafe?: boolean;
  /** The key an inline open scale is held by (the layer's key). */
  readonly key?: string;
  /** A size channel's radius (a variant's head): its scale is `sqrt` unless it says `type`. */
  readonly radius?: boolean;
  report?(code: ThemeIssueCode, message: string): void;
}

/** What a legend reads of a compiled channel. */
export type LegendSource =
  | {
      readonly kind: 'palette';
      readonly field: string;
      readonly palette: Palette;
      readonly scale?: string;
      readonly open: boolean;
      readonly label?: string;
      readonly unit?: string;
    }
  | {
      readonly kind: 'swatches';
      /** Each swatch; a letters form's names the letters it colors (the legend lists those the layer drew). */
      readonly items: readonly { readonly label: string; readonly value: number | string; readonly letters?: string }[];
    }
  | {
      readonly kind: 'sizes';
      readonly field: string;
      readonly domain: readonly [number, number];
      readonly range: readonly [number, number];
      /** The size's scale (`sqrt` for a radius when the form names none). */
      readonly type: 'linear' | 'log' | 'sqrt';
      readonly scale?: string;
    };

export interface Channel {
  readonly name: ChannelName;
  /** A uniform channel's one value: a layer reads it once. */
  readonly constant?: number | string;
  /** The value for a datum: a color (RGBA) or a number or a word; undefined: missing. */
  at(d: ThemeDatum): number | string | undefined;
  /** The open domains it reads: each key's values in view union into the held domain at settle. */
  readonly open: readonly { readonly key: string; read(d: ThemeDatum): number }[];
  readonly legend?: LegendSource;
  /** The scheme drawn (after substitution), when the form named one. */
  readonly scheme?: { readonly name: string; readonly substituted?: string };
  /**
   * The one field its value depends on, when its form is a `field` form: a layer may evaluate it once per value of a
   * dictionary-coded column (a color per code, made once a tile) and read it per mark by the code.
   */
  readonly field?: string;
}

const MaxDepth = 4;

const subjects = (alphabet: ThemeEnv['alphabet']): SchemeSubject[] =>
  alphabet === 'protein' ? ['protein'] : alphabet ? ['nucleotide'] : ['numbers', 'categories', 'secondary', 'bands'];

/**
 * A channel from a slot's body; `fallback` is the layer type's default form for the channel (drawn where the body
 * names a scheme that does not fit, and the `:else` of rules without one).
 */
export function compileChannel(
  name: ChannelName,
  body: ThemeBody | undefined,
  fallback: ThemeBody,
  env: ThemeEnv,
): Channel {
  const report = env.report ?? (() => {});
  const open: { key: string; read(d: ThemeDatum): number }[] = [];
  let legend: LegendSource | undefined;
  let scheme: Channel['scheme'];
  const kind = ThemeChannels[name];

  const literal = (v: unknown): number | string | undefined => {
    if (kind === 'color') {
      const c = resolveColor(v, env.tokens);
      if (c === undefined)
        report(
          typeof v === 'string' && v.startsWith('token:') ? 'seq.render.color' : 'seq.render.theme.color',
          `The ${name} ${JSON.stringify(v)} is no color this table resolves`,
        );
      return c;
    }
    if (kind === 'number') {
      const n = typeof v === 'number' ? v : Number.NaN;
      if (Number.isNaN(n))
        report('seq.render.theme.channel', `The ${name} channel takes a number, not ${JSON.stringify(v)}`);
      return Number.isNaN(n) ? undefined : n;
    }
    const words = ChannelWords[name as keyof typeof ChannelWords];
    if (typeof v !== 'string' || !words.includes(v)) {
      report(
        name === 'line' ? 'seq.render.line-style' : 'seq.render.theme.channel',
        `The ${name} channel takes ${words.join(', ')}, not ${JSON.stringify(v)}`,
      );
      return undefined;
    }
    return v;
  };

  type Eval = { constant?: number | string; at(d: ThemeDatum): number | string | undefined };
  const constant = (v: number | string | undefined): Eval => ({ constant: v, at: () => v });

  const compile = (b: unknown, depth: number, top: boolean): Eval => {
    if (depth > MaxDepth) return constant(undefined);
    if (!isThemeNode(b)) return constant(literal(b));
    switch (b.kind) {
      case 'uniform':
        return constant(literal(param(b, 'value')));
      case 'field':
        return field(b, top);
      case 'rules':
        return rules(b, depth);
      case 'scheme':
        return schemeForm(b, depth, top);
      case 'letters':
        return letters(b, top);
      default:
        report('seq.render.theme.channel', `A ${b.kind} is no form of the ${name} channel`);
        return constant(undefined);
    }
  };

  const scaleOf = (b: ThemeNode): { params: ScaleParams; address?: string } => {
    const s = param(b, 'scale');
    if (typeof s === 'string') {
      const params = env.scale?.(s);
      if (!params) report('seq.render.theme.channel', `The scale ${s} names no scale node`);
      return { params: params ?? {}, address: s };
    }
    const inline: Record<string, unknown> = {};
    for (const k of ['type', 'domain', 'palette', 'range', 'mid', 'label', 'unit', 'non_positive'])
      if (param(b, k) !== undefined) inline[k] = param(b, k);
    return { params: inline as ScaleParams };
  };

  const field = (b: ThemeNode, top: boolean): Eval => {
    const fieldName = String(param(b, 'name') ?? 'value');
    const { params: sc, address } = scaleOf(b);
    const read = fieldReader(fieldName);
    const key = address ?? `${env.key ?? ''}|${fieldName}`;
    const domain = resolveDomain(sc, key, env);
    if (domain.open) open.push({ key, read: (d) => toNumber(read(d)) });
    if (kind === 'number') {
      const range = sc.range ?? [0, 1];
      const type = sc.type ?? (env.radius ? 'sqrt' : 'linear');
      const [lo, hi] = domain.value;
      if (top)
        legend = {
          kind: 'sizes',
          field: fieldName,
          domain: domain.value,
          range,
          type,
          ...(address ? { scale: address } : {}),
        };
      return {
        at: (d) => {
          const t = normalize(toNumber(read(d)), type, lo, hi, undefined);
          return Number.isNaN(t) ? undefined : range[0] + (range[1] - range[0]) * Math.max(0, Math.min(1, t));
        },
      };
    }
    if (kind !== 'color') {
      // a word channel by a field: the field holds the words
      return { at: (d) => (typeof read(d) === 'string' ? literal(read(d)) : undefined) };
    }
    const node =
      (param(b, 'palette') as ThemeNode | undefined) ??
      sc.palette ??
      defaultPalette(env.numeric?.(fieldName) === true, sc.mid);
    const palette = compilePalette(isThemeNode(node) ? node : defaultPalette(false), {
      tokens: env.tokens,
      domain: domain.value,
      scale: sc,
      report: (what, v) =>
        report(
          what === 'name' ? 'seq.render.theme.channel' : 'seq.render.theme.color',
          what === 'name' ? `No palette is named ${v}` : `The palette color ${v} is no color`,
        ),
    });
    const logMissing = sc.type === 'log' && (sc.non_positive ?? 'missing') === 'missing';
    if (top)
      legend = {
        kind: 'palette',
        field: fieldName,
        palette,
        open: domain.open,
        ...(address ? { scale: address } : {}),
        ...(sc.label !== undefined ? { label: sc.label } : {}),
        ...(sc.unit !== undefined ? { unit: sc.unit } : {}),
      };
    const direct = palette.kind === 'direct';
    let warned = false;
    let logWarned = false;
    return {
      at: (d) => {
        const v = read(d);
        if (typeof v === 'number') {
          if (logMissing && v <= 0 && !logWarned) {
            logWarned = true;
            report('seq.render.scale.log', `A ${fieldName} value at or below zero is missing on a log scale`);
          }
          return palette.num(v);
        }
        const c = palette.color(v);
        if (direct && c === undefined && v !== undefined && !warned) {
          warned = true;
          report('seq.render.theme.color', `A ${fieldName} value ${JSON.stringify(v)} is no color`);
        }
        return c;
      },
    };
  };

  const fieldReader = (n: string) => {
    if (n === 'value') return (d: ThemeDatum) => d.value as string | number | undefined;
    if (n === 'letter') return (d: ThemeDatum) => d.letter as string | number | undefined;
    return (d: ThemeDatum) => d.field(n);
  };

  const measure = (m: unknown): ((d: ThemeDatum) => number | string | undefined) | undefined => {
    if (!isThemeNode(m)) return undefined;
    if (m.kind === 'field') return fieldReader(String(param(m, 'name') ?? 'value'));
    if (m.kind === 'value') return (d) => d.value;
    if (m.kind === 'letter') return (d) => d.letter;
    const entry = env.measures?.get(m.kind);
    if (!entry) {
      report('seq.render.theme.measure', `The measure ${m.kind} is not available here`);
      return () => Number.NaN;
    }
    const p = m.params ?? {};
    return (d) => entry.read(d, p);
  };

  const clause = (c: unknown): ((d: ThemeDatum) => boolean) => {
    if (!isThemeNode(c)) return () => false;
    const args = c.args ?? [];
    if (c.kind === 'in' || c.kind === 'within') {
      const ranges = typeof args[0] === 'string' ? env.locate?.(args[0]) : undefined;
      if (!ranges) {
        report('seq.render.theme.measure', `The location of an ${c.kind} clause is not known here`);
        return () => false;
      }
      const within = c.kind === 'within';
      return (d) => {
        for (let i = 0; i + 1 < ranges.length; i += 2) {
          if (within ? d.start >= ranges[i] && d.end <= ranges[i + 1] : d.start < ranges[i + 1] && d.end > ranges[i])
            return true;
        }
        return false;
      };
    }
    const m = measure(args[0]);
    if (!m) return () => false;
    const rest = args.slice(1);
    const x = rest[0];
    // a comparison across types is false; NaN fails every clause (a measure that does not apply)
    const same = (v: unknown, w: unknown) =>
      typeof v === typeof w && v === w && !(typeof v === 'number' && Number.isNaN(v));
    const num = (v: unknown) => (typeof v === 'number' && !Number.isNaN(v) ? v : undefined);
    switch (c.kind) {
      case 'eq':
        return (d) => same(m(d), x);
      case 'ne': {
        return (d) => {
          const v = m(d);
          return typeof v === typeof x && !(typeof v === 'number' && Number.isNaN(v)) && v !== x;
        };
      }
      case 'one-of':
        return (d) => {
          const v = m(d);
          for (const r of rest) if (same(v, r)) return true;
          return false;
        };
      case 'lt':
      case 'le':
      case 'gt':
      case 'ge': {
        const k = c.kind;
        const t = num(x);
        if (t === undefined) return () => false;
        return (d) => {
          const v = num(m(d));
          if (v === undefined) return false;
          return k === 'lt' ? v < t : k === 'le' ? v <= t : k === 'gt' ? v > t : v >= t;
        };
      }
      default:
        report('seq.render.theme.channel', `A ${c.kind} is no clause`);
        return () => false;
    }
  };

  const rules = (b: ThemeNode, depth: number): Eval => {
    const list = (b.args ?? []).filter(isThemeNode).map((r) => ({
      when: ((param(r, 'when') as unknown[] | undefined) ?? []).map(clause),
      value: compile(param(r, 'value'), depth + 1, false),
      label: ruleLabel(r),
    }));
    const otherwise =
      param(b, 'else') !== undefined
        ? compile(param(b, 'else'), depth + 1, false)
        : compile(fallback, depth + 1, false);
    if (depth === 0 && kind === 'color') {
      const items: { label: string; value: number | string }[] = [];
      for (const r of list) if (r.value.constant !== undefined) items.push({ label: r.label, value: r.value.constant });
      if (otherwise.constant !== undefined) items.push({ label: 'otherwise', value: otherwise.constant });
      legend = { kind: 'swatches', items };
    }
    return {
      at: (d) => {
        outer: for (let i = 0; i < list.length; i++) {
          const r = list[i];
          for (let k = 0; k < r.when.length; k++) if (!r.when[k](d)) continue outer;
          return r.value.at(d);
        }
        return otherwise.at(d);
      },
    };
  };

  const schemeForm = (b: ThemeNode, depth: number, top: boolean): Eval => {
    const n = String(param(b, 'name') ?? '');
    let entry = env.schemes.get(n);
    const why = !entry ? `No scheme is named ${n}` : unfit(entry);
    if (!entry || why) {
      report('seq.render.theme.scheme', `${why}: the layer’s default is drawn`);
      return compile(fallback, depth + 1, top);
    }
    let substituted: string | undefined;
    if (env.cvdSafe && !entry.cvdSafe && entry.substitute) {
      const sub = env.schemes.get(entry.substitute);
      if (sub && !unfit(sub)) {
        substituted = entry.name;
        entry = sub;
      }
    }
    if (top) scheme = { name: entry.name, ...(substituted ? { substituted } : {}) };
    let form = entry.expand((b.params ?? {}) as Record<string, unknown>);
    if (env.tokens.theme === 'dark' && entry.dark === 'derive') form = darkForm(form);
    return compile(form, depth + 1, top);
  };

  const unfit = (entry: SchemeEntry): string | undefined => {
    const missing = (entry.needs ?? []).filter((m) => !env.measures?.has(m));
    if (missing.length) return `The scheme ${entry.name} needs ${missing.join(', ')}`;
    if (!subjects(env.alphabet).includes(entry.subject))
      return `The scheme ${entry.name} colors ${entry.subject}, not this layer`;
    return undefined;
  };

  const letters = (b: ThemeNode, top: boolean): Eval => {
    // a table by char code (ASCII), upper and lower case alike: the hot loop indexes it
    const table: (number | string | undefined)[] = new Array(128).fill(undefined);
    const set = (l: string, v: number | string | undefined) => {
      const u = l.toUpperCase().charCodeAt(0);
      const lo = l.toLowerCase().charCodeAt(0);
      if (u < 128 && table[u] === undefined) table[u] = v;
      if (lo < 128 && table[lo] === undefined) table[lo] = v;
    };
    const items: { label: string; value: number | string; letters: string }[] = [];
    const colors = param(b, 'colors');
    if (colors && typeof colors === 'object')
      for (const [l, c] of Object.entries(colors)) {
        const v = literal(c);
        set(l, v);
        if (v !== undefined) items.push({ label: l, value: v, letters: l });
      }
    for (const g of (param(b, 'groups') as unknown[] | undefined) ?? []) {
      if (!isThemeNode(g)) continue;
      const v = literal(param(g, 'value'));
      const ls = String(param(g, 'letters') ?? '');
      for (const l of ls) set(l, v);
      if (v !== undefined) items.push({ label: String(param(g, 'label') ?? ls), value: v, letters: ls });
    }
    const dflt = param(b, 'default') !== undefined ? literal(param(b, 'default')) : undefined;
    const rest =
      dflt === undefined && param(b, 'default') === undefined ? compile(fallback, MaxDepth, false) : undefined;
    // every letter and group with its letters: the legend keeps those drawn, then one a color (`legendOf`)
    if (top && kind === 'color') legend = { kind: 'swatches', items };
    return {
      at: (d) => {
        const code = d.letter.charCodeAt(0);
        const v = code < 128 ? table[code] : undefined;
        return v !== undefined ? v : rest ? rest.at(d) : dflt;
      },
    };
  };

  const top = body ?? fallback;
  const compiled = compile(top, 0, true);
  const only =
    isThemeNode(top) && top.kind === 'field' && param(top, 'data') === undefined
      ? String(param(top, 'name') ?? 'value')
      : undefined;
  return {
    name,
    ...(only !== undefined ? { field: only } : {}),
    ...(compiled.constant !== undefined ? { constant: compiled.constant } : {}),
    at: compiled.at,
    open,
    ...(legend ? { legend } : {}),
    ...(scheme ? { scheme } : {}),
  };
}

/** A domain: the scale's ends, an open end from the held union (the data in view), `[0, 1]` until known. */
function resolveDomain(sc: ScaleParams, key: string, env: ThemeEnv): { value: [number, number]; open: boolean } {
  const d = sc.domain;
  const lo = d?.[0] ?? null;
  const hi = d?.[1] ?? null;
  if (lo !== null && hi !== null) return { value: [lo, hi], open: false };
  const held = env.domain?.(key);
  const a = lo ?? held?.[0] ?? (hi !== null ? Math.min(0, hi) : 0);
  const b = hi ?? held?.[1] ?? Math.max(a + 1, 1);
  return { value: [a, b], open: true };
}

function toNumber(v: string | number | undefined): number {
  return typeof v === 'number' ? v : v === undefined ? Number.NaN : Number.parseFloat(v);
}

const ClauseWords: Readonly<Record<string, string>> = {
  eq: '=',
  ne: '≠',
  lt: '<',
  le: '≤',
  gt: '>',
  ge: '≥',
  'one-of': 'one of',
  in: 'in',
  within: 'within',
};

/** A rule's legend label: its clauses as people read them (`plddt ≥ 90`), or its own `label`. */
function ruleLabel(r: ThemeNode): string {
  if (typeof param(r, 'label') === 'string') return param(r, 'label') as string;
  const parts: string[] = [];
  for (const c of (param(r, 'when') as unknown[] | undefined) ?? []) {
    if (!isThemeNode(c)) continue;
    const [m, ...rest] = c.args ?? [];
    const what = isThemeNode(m) ? (m.kind === 'field' ? String(param(m, 'name')) : m.kind) : String(m);
    parts.push(
      `${what} ${ClauseWords[c.kind] ?? c.kind} ${rest.map((x) => (typeof x === 'string' ? x : JSON.stringify(x))).join(', ')}`,
    );
  }
  return parts.join(' and ');
}

/** A form's literal colors in the dark band (a scheme's derived dark variant); tokens and `light-dark` kept. */
function darkForm<T>(v: T): T {
  if (typeof v === 'string') return (v.startsWith('#') ? darkVariant(v) : v) as T;
  if (Array.isArray(v)) return v.map(darkForm) as T;
  if (v && typeof v === 'object') {
    const out: Record<string, unknown> = {};
    for (const [k, x] of Object.entries(v))
      out[k] = k === 'kind' || k === 'name' || k === 'letters' || k === 'label' ? x : darkForm(x);
    return out as T;
  }
  return v;
}

/**
 * An `ink` channel's color held at 3:1 against its fill (_spec/README.md "Themes", Colors): moved in Oklab lightness,
 * away from the fill's, until it reaches 3:1 (pure black or white at the end).
 */
export function inkAgainst(ink: number, fill: number): number {
  if (contrastRgba(ink, fill) >= 3) return ink;
  const fl = Oklab.lightness(toRgb(fill));
  const dir = fl > 0.6 ? -1 : 1;
  let l = Oklab.lightness(toRgb(ink));
  for (let i = 0; i < 40; i++) {
    l += dir * 0.025;
    if (l <= 0 || l >= 1) break;
    const c = fromRgb(Oklab.withLightness(toRgb(ink), l));
    if (contrastRgba(c, fill) >= 3) return c;
  }
  return dir < 0 ? 0x000000ff : 0xffffffff;
}

/** A color whose Oklab chroma is below this is a grey (a gap, `missing`, a neutral): `glyphAgainst` keeps it grey. */
export const GreyChroma = 0.02;

/**
 * A color made to read as a glyph on a surface (_spec/README.md "Render", The logo): its Oklab chroma raised to
 * `chroma`, hue kept, clipped to the sRGB gamut at its lightness, then its lightness moved away from the fill's until
 * 3:1 (`inkAgainst`'s rule; pure black or white at the end). A color above the floor keeps its chroma; a grey
 * (chroma below `GreyChroma`) stays grey, so only its lightness moves.
 */
export function glyphAgainst(ink: number, fill: number, chroma: number): number {
  const lab = Oklab.fromColor(Oklab.zero(), toRgb(ink));
  const c = Math.hypot(lab[1], lab[2]);
  if (c < GreyChroma || c >= chroma) return inkAgainst(ink, fill);
  const k = chroma / c;
  const fl = Oklab.lightness(toRgb(fill));
  const dir = fl > 0.6 ? -1 : 1;
  const l0 = lab[0];
  const out = Oklab.zero();
  for (let i = 0; i < 40; i++) {
    const l = l0 + dir * 0.025 * i;
    if (l <= 0 || l >= 1) break;
    out[0] = l;
    out[1] = lab[1] * k;
    out[2] = lab[2] * k;
    const col = fromRgb(Oklab.toColor(Oklab.clampChroma(out, out)));
    if (contrastRgba(col, fill) >= 3) return col;
  }
  return dir < 0 ? 0x000000ff : 0xffffffff;
}
