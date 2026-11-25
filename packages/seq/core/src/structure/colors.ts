/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's colors as runs on units (_spec/README.md "Structures", Colors): its color channel as its track compiles it
// (`LayerColors`: the renderer's theme context, its held domains and the reader's substitution, the layer type's
// default and a colourless layer's categorical place), with the light and the dark tokens, evaluated per position
// (letters, values) or per item (ranges), grouped by color and carried onto the bound units. Unobserved positions
// cut, uncolored ones in no run, data not held pending, never an empty run (EC-169).

import { type Atom, type Issue, issueOf } from '@molstar/core/repr/tree';
import { type Location, location, rangesLocation } from '@molstar/seqx-core/query';
import { builtInTokens, type LayerColors, layerSource, type SeqTokens } from '@molstar/seqx-core/render';
import type { Channel, ThemeDatum } from '@molstar/seqx-core/theme';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { StructureSyncIssues } from './issues.js';
import { type StructureBinding, type StructureModel, structureReads } from './model.js';
import { annotationItems, type HeldOptions, readHeld, type StoreRead } from './rows.js';
import { type StructureTarget, toUnits } from './units.js';

/** One color and where it is, on units. */
export interface ColorRun {
  /** RGBA (`0xRRGGBBAA`) with the light tokens. */
  readonly color: number;
  /** The same with the dark tokens (a scheme's dark variant). */
  readonly dark: number;
  readonly target: StructureTarget;
}

export interface ColorRunsResult {
  /** Nothing pending: every position of the layer's data was read. */
  readonly complete: boolean;
  /** Positions of the layer's data the store does not hold yet. */
  readonly pending: Location;
  /** Whether the store refused the data's tiles as past its bound (their positions pending). */
  readonly tooMany: boolean;
  /** By the first position each color has. */
  readonly runs: readonly ColorRun[];
  readonly issues: readonly Issue[];
}

export interface ColorRunsOptions extends HeldOptions {
  /**
   * The view's layer colors as its renderer compiles them (`Renderer.layerColors()`); a headless host without a
   * renderer makes them with `layerColors(view, viewTheme(view, theme, domains))`.
   */
  readonly colors: LayerColors;
  /** The tokens of each theme; absent, the built-in light and dark ones. */
  readonly tokens?: { readonly light: SeqTokens; readonly dark: SeqTokens };
}

const BuiltIn = { light: builtInTokens('light'), dark: builtInTokens('dark') };

/** A ranges item's datum: its fields as the layer's rows would give them. */
class Datum implements ThemeDatum {
  value = Number.NaN;
  letter = '';
  start = 0;
  end = 0;
  fields: Readonly<Record<string, Atom>> = {};
  field(name: string): string | number | undefined {
    if (name === 'value') return this.value;
    if (name === 'letter') return this.letter;
    const v = this.fields[name];
    return v === null || v === undefined ? undefined : typeof v === 'number' || typeof v === 'string' ? v : String(v);
  }
}

interface Group {
  readonly color: number;
  readonly dark: number;
  readonly at: number;
  readonly ranges: number[];
  readonly parts: Location[];
}

/**
 * The colors of a letters, values or ranges layer (a resolved view's) on the units of the bound complexes. A snapshot:
 * a host rebuilds on the layer's changes, or the held domains', when it wants live colors.
 */
export function colorRuns(
  model: StructureModel,
  store: StoreRead,
  layer: ResolvedLayer,
  complexes: readonly StructureBinding[] | undefined,
  options: ColorRunsOptions,
): ColorRunsResult {
  const r = structureReads(model);
  const data = layer.data;
  const tokens = options.tokens ?? BuiltIn;
  const type = [layer.type, ...layer.fallbacks].find((t) => t === 'letters' || t === 'values' || t === 'ranges');
  const none = (why: string): ColorRunsResult => ({
    complete: true,
    pending: rangesLocation(data?.system ?? '', []),
    tooMany: false,
    runs: [],
    issues: [issueOf(StructureSyncIssues, 'seq.structure.layer', why, { address: layer.address ?? '' })],
  });
  if (!data || !type) return none(`The layer ${layer.address ?? ''} is no letters, values or ranges layer with data`);
  const datum = new Datum();
  const groups = new Map<string, Group>();
  const put = (light: unknown, dark: unknown, at: number, s: number, e: number, loc?: Location) => {
    if (typeof light !== 'number' || typeof dark !== 'number') return;
    const k = `${light}:${dark}`;
    let g = groups.get(k);
    if (!g) groups.set(k, (g = { color: light, dark, at, ranges: [], parts: [] }));
    if (loc) g.parts.push(loc);
    else if (g.ranges.length && g.ranges[g.ranges.length - 1] === s) g.ranges[g.ranges.length - 1] = e;
    else g.ranges.push(s, e);
  };
  const channels = (extra: Parameters<LayerColors['channel']>[2]): [Channel, Channel] | undefined => {
    const light = options.colors.channel(layer, tokens.light, extra);
    const dark = options.colors.channel(layer, tokens.dark, extra);
    return light && dark ? [light, dark] : undefined;
  };

  if (type === 'ranges') {
    const read = annotationItems(r, model, store, data.address, options.want);
    if (!read) return none(`The layer ${layer.address ?? ''} reads no ranges or pairs annotation`);
    const ch = channels({ numeric: (n: string) => read.numeric.has(n) });
    if (!ch) return none(`The layer ${layer.address ?? ''} is in no track of the colors' view`);
    for (const it of read.items) {
      datum.fields = it.fields;
      datum.start = it.start;
      datum.end = it.end;
      put(ch[0].at(datum), ch[1].at(datum), it.start, 0, 0, location({ parts: it.parts }));
    }
    return finish(model, groups, read.pending, read.tooMany, !read.unknown, complexes, options);
  }

  const entry = r.system(data.system);
  let system = data.system;
  let length = entry?.length ?? 0;
  let letters: ((p: number) => string) | undefined;
  let values: Float32Array | undefined;
  let pending: number[] = [];
  let tooMany = false;
  let known = true;
  if (type === 'letters' && entry?.letters !== undefined) {
    const l = entry.letters;
    letters = (p) => l[p] ?? '';
  } else {
    const src = layerSource(r.model, layer);
    if (!src) return none(`The layer ${layer.address ?? ''} reads no data the store holds`);
    system = src.dataSystem;
    length = src.length ?? 0;
    const held = readHeld(store, src, options.want);
    pending = held.pending;
    tooMany = held.tooMany;
    known = !held.unknown;
    if (type === 'letters') {
      const byStart = new Map<number, string>();
      for (const e of held.entries) {
        const t = e.tables.letters;
        if (!t) continue;
        const s = typeof t.letters === 'string' ? t.letters : codes(t.letters, t.alphabet);
        for (let i = 0; i < s.length; i++) byStart.set(t.start + i, s[i]);
      }
      letters = (p) => byStart.get(p) ?? '';
    } else {
      values = new Float32Array(length).fill(Number.NaN);
      for (const e of held.entries) {
        const t = e.tables.values;
        if (t?.form !== 'positions') continue;
        for (let i = 0; i < t.value.length; i++) if (t.start + i < length) values[t.start + i] = t.value[i];
      }
    }
  }
  const alphabet = r.system(system)?.alphabet ?? entry?.alphabet ?? 'protein';
  const ch = channels(type === 'letters' ? { alphabet } : {});
  if (!ch) return none(`The layer ${layer.address ?? ''} is in no track of the colors' view`);
  for (let p = 0; p < length; p++) {
    datum.start = p;
    datum.end = p + 1;
    if (letters) {
      datum.letter = letters(p);
      if (!datum.letter) continue;
    } else {
      datum.value = values![p];
      if (Number.isNaN(datum.value)) continue;
    }
    put(ch[0].at(datum), ch[1].at(datum), p, p, p + 1);
  }
  for (const g of groups.values()) g.parts.push(rangesLocation(system, g.ranges));
  return finish(model, groups, rangesLocation(system, pending), tooMany, known, complexes, options);
}

function codes(c: Uint8Array, alphabet: string | undefined): string {
  let s = '';
  for (let i = 0; i < c.length; i++) s += alphabet ? (alphabet[c[i]] ?? '') : String.fromCharCode(c[i]);
  return s;
}

function finish(
  model: StructureModel,
  groups: Map<string, Group>,
  pending: Location,
  tooMany: boolean,
  known: boolean,
  complexes: readonly StructureBinding[] | undefined,
  options: HeldOptions,
): ColorRunsResult {
  const runs: (ColorRun & { at: number })[] = [];
  const issues: Issue[] = [];
  for (const g of groups.values()) {
    const t = toUnits(model, location({ parts: g.parts.flatMap((l) => l.parts) }), complexes, {
      limits: null,
      ...(options.ligand ? { ligand: options.ligand } : {}),
    });
    for (const i of t.issues) if (!issues.some((j) => j.code === i.code)) issues.push(i);
    if (t.target.parts.length) runs.push({ color: g.color, dark: g.dark, target: t.target, at: g.at });
  }
  runs.sort((a, b) => a.at - b.at);
  return {
    complete: known && !pending.parts.length,
    pending,
    tooMany,
    runs: runs.map(({ at: _, ...run }) => run),
    issues,
  };
}
