/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's colors as runs on units (_spec/README.md "Structures", Colors): its color channel compiled as drawing
// compiles it (its slot's body, else the layer type's default, the built-in tokens light and dark, an open scale's
// domain the data's extent), evaluated per position (letters, values) or per item (ranges), grouped by color and
// carried onto the bound units. Unobserved positions cut, uncolored ones in no run, data not held pending (EC-169).

import { type Atom, type Issue, issueOf } from '@molstar/core/repr/tree';
import { type Location, location, rangesLocation } from '@molstar/seqx-core/query';
import { builtInTokens, DefaultThemeContext, layerSource, type SeqTokens } from '@molstar/seqx-core/render';
import {
  type Channel,
  compileChannel,
  lettersDefault,
  type ThemeBody,
  type ThemeDatum,
  type ThemeEnv,
} from '@molstar/seqx-core/theme';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { StructureSyncIssues } from './issues.js';
import { type StructureBinding, type StructureModel, structureReads } from './model.js';
import { annotationItems, type HeldOptions, readHeld, type StoreRead } from './rows.js';
import { type StructureTarget, toUnits } from './units.js';

/** One color and where it is, on units. */
export interface ColorRun {
  /** RGBA (`0xRRGGBBAA`) with the light tokens. */
  readonly color: number;
  /** The same with the dark tokens (a dark variant of a scheme). */
  readonly dark: number;
  readonly target: StructureTarget;
}

export interface ColorRunsResult {
  readonly complete: boolean;
  /** Positions of the layer's data the store does not hold yet. */
  readonly pending: Location;
  /** By the first position each color has. */
  readonly runs: readonly ColorRun[];
  readonly issues: readonly Issue[];
}

const Light = builtInTokens('light');
const Dark = builtInTokens('dark');

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

/**
 * The colors of a letters, values or ranges layer (a resolved view's) on the units of the bound complexes. A snapshot:
 * a host rebuilds on the layer's changes when it wants live colors.
 */
export function colorRuns(
  model: StructureModel,
  store: StoreRead,
  layer: ResolvedLayer,
  complexes?: readonly StructureBinding[],
  options: HeldOptions = {},
): ColorRunsResult {
  const r = structureReads(model);
  const data = layer.data;
  const type = [layer.type, ...layer.fallbacks].find((t) => t === 'letters' || t === 'values' || t === 'ranges');
  const none = (why: string): ColorRunsResult => ({
    complete: true,
    pending: rangesLocation(data?.system ?? '', []),
    runs: [],
    issues: [issueOf(StructureSyncIssues, 'seq.structure.layer', why, { address: layer.address ?? '' })],
  });
  if (!data || !type) return none(`The layer ${layer.address ?? ''} is no letters, values or ranges layer with data`);
  const body = layer.theme.color as ThemeBody | undefined;
  const datum = new Datum();
  // per color key: the light and dark colors and their ranges on the data's system, or their items' parts
  const groups = new Map<string, { color: number; dark: number; at: number; ranges: number[]; parts: Location[] }>();
  const put = (light: unknown, dark: unknown, at: number, s: number, e: number, loc?: Location) => {
    if (typeof light !== 'number' || typeof dark !== 'number') return;
    const k = `${light}:${dark}`;
    let g = groups.get(k);
    if (!g) groups.set(k, (g = { color: light, dark, at, ranges: [], parts: [] }));
    if (loc) g.parts.push(loc);
    else if (g.ranges.length && g.ranges[g.ranges.length - 1] === s) g.ranges[g.ranges.length - 1] = e;
    else g.ranges.push(s, e);
  };
  let pending: number[] = [];
  let system = data.system;

  if (type === 'ranges') {
    const read = annotationItems(r, model, store, data.address, options.want);
    if (!read) return none(`The layer ${layer.address ?? ''} reads no ranges or pairs annotation`);
    const env = { numeric: (n: string) => read.numeric.has(n) };
    const [light, dark] = channels(body, 'token:range', env, (d) => {
      for (const it of read.items) {
        setItem(datum, it);
        d(datum);
      }
    });
    for (const it of read.items) {
      setItem(datum, it);
      put(light.at(datum), dark.at(datum), it.start, 0, 0, location({ parts: it.parts }));
    }
    return finish(model, groups, read.pending, complexes, options);
  }

  const entry = r.system(data.system);
  let letters: ((p: number) => string) | undefined;
  let values: Float32Array | undefined;
  let length = entry?.length ?? 0;
  if (type === 'letters' && entry?.letters !== undefined) {
    const l = entry.letters;
    letters = (p) => l[p] ?? '';
  } else {
    const src = layerSource(r.model, layer);
    if (!src) return none(`The layer ${layer.address ?? ''} reads no data the store holds`);
    system = src.dataSystem;
    length = src.length ?? 0;
    const held = readHeld(store, src.source, src.system, src.length, options.want);
    pending = held.pending;
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
  const [light, dark] =
    type === 'letters'
      ? channels(body, { kind: 'scheme', params: { name: lettersDefault(alphabet) } }, { alphabet }, () => {})
      : channels(body, 'token:value', { numeric: (n) => n === 'value' }, (d) => {
          for (let p = 0; p < length; p++) {
            datum.value = values![p];
            d(datum);
          }
        });
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
    put(light.at(datum), dark.at(datum), p, p, p + 1);
  }
  for (const g of groups.values()) g.parts.push(rangesLocation(system, g.ranges));
  return finish(
    model,
    groups,
    pending.length ? rangesLocation(system, pending) : rangesLocation(system, []),
    complexes,
    options,
  );
}

function setItem(d: Datum, it: { fields: Readonly<Record<string, Atom>>; start: number; end: number }) {
  d.fields = it.fields;
  d.start = it.start;
  d.end = it.end;
  const v = it.fields.value;
  d.value = typeof v === 'number' ? v : Number.NaN;
}

function codes(c: Uint8Array, alphabet: string | undefined): string {
  let s = '';
  for (let i = 0; i < c.length; i++) s += alphabet ? (alphabet[c[i]] ?? '') : String.fromCharCode(c[i]);
  return s;
}

/**
 * The color channel with the light and the dark tokens; an open scale's domain the extent `each` folds (a two-pass
 * compile: the first finds the open domains, the second reads them, as a view's settle holds them).
 */
function channels(
  body: ThemeBody | undefined,
  fallback: ThemeBody,
  extra: Partial<ThemeEnv>,
  each: (fold: (d: ThemeDatum) => void) => void,
): [Channel, Channel] {
  const compile = (tokens: SeqTokens, domain?: (key: string) => readonly [number, number] | undefined) =>
    compileChannel('color', body, fallback, {
      tokens,
      schemes: DefaultThemeContext.schemes,
      cvdSafe: false,
      ...(domain ? { domain } : {}),
      ...extra,
    });
  const first = compile(Light);
  if (!first.open.length) return [first, compile(Dark)];
  const lo = new Map<string, number>();
  const hi = new Map<string, number>();
  each((d) => {
    for (const o of first.open) {
      const v = o.read(d);
      if (Number.isNaN(v)) continue;
      lo.set(o.key, Math.min(lo.get(o.key) ?? v, v));
      hi.set(o.key, Math.max(hi.get(o.key) ?? v, v));
    }
  });
  const domain = (key: string): readonly [number, number] | undefined =>
    lo.has(key) ? [lo.get(key)!, hi.get(key)!] : undefined;
  return [compile(Light, domain), compile(Dark, domain)];
}

function finish(
  model: StructureModel,
  groups: Map<string, { color: number; dark: number; at: number; parts: Location[] }>,
  pending: Location,
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
    complete: !pending.parts.length,
    pending,
    runs: runs.map(({ at: _, ...run }) => run),
    issues,
  };
}
