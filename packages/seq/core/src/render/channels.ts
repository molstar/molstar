/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column, issueOf } from '@molstar/core/repr/tree';
import type { RangesField } from '@molstar/seqx-core/data';
import {
  builtInSchemes,
  type Channel,
  type ChannelName,
  compileChannel,
  type Legend,
  legendOf,
  type MeasureEntry,
  type ScaleParams,
  type SchemeRegistry,
  type ThemeBody,
  type ThemeDatum,
  type ThemeEnv,
  type ThemeIssueCode,
} from '@molstar/seqx-core/theme/index';
import { RenderIssues } from './issues.js';
import type { LayoutContext } from './layer-type.js';

// The channels a layer type reads (_spec/README.md "Render", Layer types; "Themes"): each compiled from its slot's body
// once per token table, held domains and the reader's substitution, and kept in the layer's prepared form; a layer
// reads a uniform channel once and an evaluated one per mark through a datum it keeps.

/** What a layout reads of the themes: the registry, the measures, the scales, the held domains, the reader's setting. */
export interface ThemeContext {
  readonly schemes: SchemeRegistry;
  readonly measures?: ReadonlyMap<string, MeasureEntry>;
  /** The reader's color-blind substitution (`theme.cvdSafe`): never the spec's. */
  readonly cvdSafe: boolean;
  /** A scale node's params by its address. */
  scale(address: string): ScaleParams | undefined;
  /** An open scale's domain held since the view last settled (unioned over the layers naming it). */
  domain(key: string): readonly [number, number] | undefined;
  /** Folds a layer's extent of an open scale's values in view into this layout's union. */
  extent(key: string, lo: number, hi: number): void;
}

const Schemes = builtInSchemes();

/** A theme context without scales or held domains (a layout without a view's: tests, a lone layer). */
export const DefaultThemeContext: ThemeContext = {
  schemes: Schemes,
  cvdSafe: false,
  scale: () => undefined,
  domain: () => undefined,
  extent: () => {},
};

/** The render codes a theme says (cataloged in `RenderIssues`). */
const ThemeCodes: Readonly<Record<ThemeIssueCode, true>> = {
  'seq.render.theme.channel': true,
  'seq.render.theme.scheme': true,
  'seq.render.theme.measure': true,
  'seq.render.theme.color': true,
  'seq.render.scale.log': true,
  'seq.render.line-style': true,
  'seq.render.color': true,
};

/** A theme code reported for the layer at `ctx.key`, once a layout (the layout dedupes by code and address). */
export function reportTheme(ctx: LayoutContext, code: ThemeIssueCode, message: string): void {
  if (!ThemeCodes[code]) return;
  ctx.report(issueOf(RenderIssues, code, message, { address: ctx.key, hint: RenderIssues.codes[code].hint }));
}

interface Kept {
  readonly channel: Channel;
  readonly version: number;
  readonly cvd: boolean;
  readonly domains: string;
}

/** A layer's channels, compiled from its slots and kept while the tokens, the held domains and the reader agree. */
export class LayerChannels {
  private readonly kept = new Map<ChannelName, Kept>();
  constructor(readonly forms: Readonly<Partial<Record<ChannelName, unknown>>>) {}

  /** The channel `name` (its slot's body, else `fallback`, the layer type's default) for this layout. */
  get(name: ChannelName, fallback: ThemeBody, ctx: LayoutContext, extra: Partial<ThemeEnv> = {}): Channel {
    const theme = ctx.theme;
    const k = this.kept.get(name);
    if (k && k.version === ctx.tokens.version && k.cvd === theme.cvdSafe && k.domains === held(k.channel, theme))
      return k.channel;
    const channel = compileChannel(name, this.forms[name] as ThemeBody | undefined, fallback, {
      tokens: ctx.tokens,
      schemes: theme.schemes,
      ...(theme.measures ? { measures: theme.measures } : {}),
      scale: (a) => theme.scale(a),
      domain: (key) => theme.domain(key),
      cvdSafe: theme.cvdSafe,
      key: ctx.key,
      report: (code, message) => reportTheme(ctx, code, message),
      ...extra,
    });
    this.kept.set(name, { channel, version: ctx.tokens.version, cvd: theme.cvdSafe, domains: held(channel, theme) });
    return channel;
  }
}

function held(channel: Channel, theme: ThemeContext): string {
  let s = '';
  for (const o of channel.open) s += `${o.key}=${theme.domain(o.key)?.join(',') ?? ''};`;
  return s;
}

/** A layer's channel forms as its resolved layer holds them. */
export function channelForms(theme: Readonly<Record<string, unknown>>): LayerChannels {
  return new LayerChannels(theme as Partial<Record<ChannelName, unknown>>);
}

/**
 * A datum over a ranges table's rows: `row` set per mark, a field's column looked up once a table. A category reads
 * its text, a number column its number (NaN missing), a list its item.
 */
export class RowDatum implements ThemeDatum {
  value = Number.NaN;
  letter = '';
  start = 0;
  end = 0;
  row = 0;
  private fields: readonly RangesField[] | undefined;
  private readonly cols = new Map<string, Column | undefined>();

  bind(fields: readonly RangesField[] | undefined): void {
    if (fields === this.fields) return;
    this.fields = fields;
    this.cols.clear();
  }

  column(name: string): Column | undefined {
    let c = this.cols.get(name);
    if (c === undefined && !this.cols.has(name)) {
      c = this.fields?.find((f) => f.name === name)?.column;
      this.cols.set(name, c);
    }
    return c;
  }

  field(name: string): string | number | undefined {
    const c = this.column(name);
    if (!c) return undefined;
    if (c instanceof CategoryColumn) return c.get(this.row) ?? undefined;
    const v = (c as ArrayLike<unknown>)[this.row];
    return typeof v === 'number' ? v : v === null || v === undefined ? undefined : String(v);
  }

  /** Whether a field's column holds numbers (its declared type: a typed array), for the default palette. */
  numeric(name: string): boolean {
    const c = this.column(name);
    return c !== undefined && !(c instanceof CategoryColumn) && !Array.isArray(c);
  }
}

/** A datum with a value or a letter alone (values and letters layers). */
export class PlainDatum implements ThemeDatum {
  value = Number.NaN;
  letter = '';
  start = 0;
  end = 0;
  field(name: string): string | number | undefined {
    return name === 'value' ? this.value : name === 'letter' ? this.letter : undefined;
  }
}

/** The extents of a channel's open domains this layout: folded per mark, handed to the context at its end. */
export class Extents {
  private lo = new Float64Array(4);
  private hi = new Float64Array(4);
  private keys: readonly { readonly key: string; read(d: ThemeDatum): number }[] = [];

  /** Whether any channel reads an open domain (a layout folds nothing otherwise). */
  any = false;

  reset(channels: readonly (Channel | undefined)[]): void {
    const keys: { readonly key: string; read(d: ThemeDatum): number }[] = [];
    for (const c of channels) if (c) keys.push(...c.open);
    this.keys = keys;
    this.any = keys.length > 0;
    if (this.lo.length < this.keys.length) {
      this.lo = new Float64Array(this.keys.length);
      this.hi = new Float64Array(this.keys.length);
    }
    this.lo.fill(Number.POSITIVE_INFINITY);
    this.hi.fill(Number.NEGATIVE_INFINITY);
  }

  add(d: ThemeDatum): void {
    for (let i = 0; i < this.keys.length; i++) {
      const v = this.keys[i].read(d);
      if (Number.isNaN(v)) continue;
      if (v < this.lo[i]) this.lo[i] = v;
      if (v > this.hi[i]) this.hi[i] = v;
    }
  }

  flush(theme: ThemeContext): void {
    for (let i = 0; i < this.keys.length; i++)
      if (this.lo[i] <= this.hi[i]) theme.extent(this.keys[i].key, this.lo[i], this.hi[i]);
  }
}

/** A layer's legends from its channels (those with a key to say), titled by its track. */
export function channelLegends(
  channels: readonly (Channel | undefined)[],
  options: { readonly title: string; readonly layer: string; readonly painted?: Iterable<string> },
): Legend[] {
  const out: Legend[] = [];
  for (const c of channels) {
    if (!c) continue;
    const l = legendOf(c, options);
    if (l) out.push(l);
  }
  return out;
}
