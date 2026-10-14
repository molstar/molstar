/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `rules` theme (seqspec:viewer#extension-points; seqspec:format#themes): thresholds and composition as data, for
// every channel kind. Its rules compile at `create` into flat arrays (per rule its clauses, value, nested theme and
// alpha measure; per clause its measure, bounds and letter mask), and its measures into typed arrays: a field per data
// row, a column's value, the reference row's letters, and one composition per column (letter counts over the
// alignment's held rows,
// or its `profile`'s frequencies) that `share` and `rank` read. So `value` reads typed arrays and the location's
// numbers and allocates nothing. On a sequence (no alignment) a position is a column of its own letter: Clustal X
// colors a lone protein as the old theme did. `share` or `rank` on a lazy alignment with no profile is not applicable:
// the reference's fallback, else `letter`, is drawn. The named tables (`clustal_x`, `identity_to_row`,
// `conservation_shading`) are `rules` expanded at resolve.

import {
  type ChannelEntry,
  type ChannelKind,
  Color,
  canonical,
  numberColumn,
  paletteValue,
  type ThemeRef,
} from "@molstar/seq-core/format";
import { type Annotation, Annotations, CoordinateSystemKind, KeyKind, type KeyTarget } from "@molstar/seq-core/model";
import { SchemeTheme } from "./color-schemes.ts";
import { NO_COLOR } from "./draw/backend.ts";
import {
  coverGranularity,
  dataIndex,
  Granularity as G,
  KernelThemes,
  NO_GLYPH,
  NO_SIZE,
  OPAQUE,
  SizeUnit,
  sizeOf,
  type Theme,
  type ThemeContext,
  type ThemeLocation,
  type ThemeProvider,
  themeProvider,
} from "./themes.ts";

interface Clause {
  readonly measure?: unknown;
  readonly min?: number | null;
  readonly max?: number | null;
  readonly in?: readonly string[] | null;
}
interface Rule {
  readonly when?: readonly Clause[];
  readonly theme?: ThemeRef | null;
  readonly alpha?: unknown;
  readonly label?: string | null;
}
type Params = Readonly<Record<string, unknown>>;

/** Measure codes: the element's value, a row field (or a string field's membership in `in`), a column's value, the
 * share of its letter or of a set, its letter's rank, whether it equals a row's letter. */
const M = { value: 0, field: 1, fieldIn: 2, column: 3, share: 4, shareOf: 5, rank: 6, equalsRow: 7 } as const;
const MEASURE_GRANULARITY = [G.value, G.item, G.item, G.column, G.column, G.column, G.column, G.column];

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const isObject = (m: unknown): m is Record<string, unknown> => m !== null && typeof m === "object";
const usesColumns = (m: unknown) => isObject(m) && ("share" in m || "rank" in m);
/** The section's alignment (an alignment axis has one block, the columns); -1: none. */
const alignmentOf = (tc: ThemeContext) => tc.section.blocks[0]?.alignment ?? -1;
const annotationOf = (tc: ThemeContext, ref: unknown): Annotation | null =>
  typeof ref === "string" && tc.set.byName(ref, t) && t.kind === KeyKind.annotation
    ? (tc.set.annotations[t.id] ?? null)
    : null;
/** The letter's index 0..25 (case-blind); -1: not a letter. */
const letterOf = (c: number) => ((c | 32) >= 97 && (c | 32) <= 122 ? (c & ~32) - 65 : -1);

/** Letter counts per column (26 a column) and their totals: the alignment's held rows, or a profile's frequencies. */
function composition(tc: ThemeContext, profile: unknown): { counts: Float64Array; totals: Float64Array } | null {
  const aln = tc.set.alignments[alignmentOf(tc)];
  const p = annotationOf(tc, profile);
  if (!aln && !p) return null;
  const n = aln?.columnCount ?? 0;
  const counts = new Float64Array(n * 26);
  const totals = new Float64Array(n);
  const add = (c: number, k: number, w: number) => {
    if (k < 0 || c < 0 || c >= n) return;
    counts[c * 26 + k] += w;
    totals[c] += w;
  };
  if (p?.matrix) {
    const S = p.matrix.states.length;
    const letters = Array.from(p.matrix.states, (id) => letterOf(tc.strings.get(id).charCodeAt(0)));
    for (let r = 0; r < p.count; r++)
      for (let s = 0; s < S; s++) add(p.start[r], letters[s], p.matrix.values[r * S + s]);
  } else if (aln && !aln.lazy) {
    const codes = aln.rows.codes;
    for (let r = 0; r < aln.rows.count; r++) for (let c = 0; c < n; c++) add(c, letterOf(codes[r * n + c]), 1);
  }
  return { counts, totals };
}

/** The measures of one `rules` theme: deduplicated, each a code and its data, read through one function. */
class Measures {
  readonly code: number[] = [];
  /** field: per data annotation, per row; column: per column; equalsRow: the row's codes; shareOf: letter indices. */
  readonly rows: (Float64Array | null)[][] = [];
  readonly data: (Float64Array | Uint8Array | null)[] = [];
  readonly keys = new Map<string, number>();
  comp: { counts: Float64Array; totals: Float64Array } | null | undefined;
  readonly tc: ThemeContext;
  readonly profile: unknown;
  constructor(tc: ThemeContext, profile: unknown) {
    this.tc = tc;
    this.profile = profile;
  }

  /** A string field (no number column in the layer's data): `in` tests its value, not the letter. */
  stringField(m: unknown): boolean {
    if (!isObject(m) || typeof m.field !== "string") return false;
    return !this.tc.annotations.some((a) => a && numberColumn(a, this.tc.set, m.field as string, this.tc.svs));
  }
  add(m: unknown, inSet: readonly string[] | null): number {
    const key = canonical([m, inSet]);
    const held = this.keys.get(key);
    if (held !== undefined) return held;
    const { tc } = this;
    const k = this.code.length;
    let code: number = M.value;
    let data: Float64Array | Uint8Array | null = null;
    let rows: (Float64Array | null)[] = [];
    const aln = tc.set.alignments[alignmentOf(tc)];
    if (isObject(m) && typeof m.field === "string") {
      const name = m.field;
      code = inSet ? M.fieldIn : M.field;
      const of = new Set(inSet);
      rows = tc.annotations.map((a) => {
        const column = a && !inSet ? numberColumn(a, tc.set, name, tc.svs) : null;
        return Float64Array.from({ length: a?.count ?? 0 }, (_, r) => {
          if (column) return column[r];
          const v = paletteValue(a, tc.set, r, name, tc.svs);
          if (inSet) return v == null ? Number.NaN : +of.has(String(v));
          return typeof v === "number" ? v : typeof v === "string" ? Number.parseFloat(v) : Number.NaN;
        });
      });
    } else if (isObject(m) && "column" in m) {
      code = M.column;
      const a = annotationOf(tc, m.column);
      const n = aln?.columnCount ?? 0;
      const slot = a ? Annotations.slotOf(a, { kind: CoordinateSystemKind.columns, id: aln?.id ?? -1, row: -1 }) : -1;
      // No slot on this alignment's columns: no data, so the measure is NaN (never another target's values).
      data = Float64Array.from({ length: a && slot >= 0 ? n : 0 }, (_, c) =>
        Annotations.valueAt(a as Annotation, slot, c),
      );
    } else if (isObject(m) && "equals_row" in m) {
      code = M.equalsRow;
      const row = typeof m.equals_row === "number" ? m.equals_row : Math.max(0, aln?.referenceRow ?? 0);
      // A lazy alignment's row where it is held (the reference row is fetched); else none (`isApplicable` says no).
      const slot = aln ? aln.rows.slotOf(row as never) : -1;
      data = aln && slot >= 0 ? aln.rows.codes.subarray(slot * aln.columnCount, (slot + 1) * aln.columnCount) : null;
    } else if (isObject(m) && ("share" in m || "rank" in m)) {
      const of = typeof m.share === "string" && m.share !== "letter" ? m.share : null;
      code = "rank" in m ? M.rank : of !== null ? M.shareOf : M.share;
      if (of !== null)
        data = Uint8Array.from(new Set([...of].map((c) => letterOf(c.charCodeAt(0))).filter((x) => x >= 0)));
      if (this.comp === undefined) this.comp = composition(tc, this.profile);
    }
    this.code.push(code);
    this.data.push(data);
    this.rows.push(rows);
    this.keys.set(key, k);
    return k;
  }
}

/** One `rules` theme: compiled; `value` evaluates the first rule that holds. */
function compile(tc: ThemeContext, p: Params, channel: ChannelEntry) {
  const kind = channel.kind;
  const rules = (p.rules ?? []) as readonly Rule[];
  const R = rules.length;
  const ms = new Measures(tc, p.profile ?? null);
  const C = rules.reduce((n, r) => n + (r.when?.length ?? 0), 0);
  const cs = new Int32Array(R + 1);
  const cm = new Int32Array(C).fill(-1);
  const lo = new Float64Array(C).fill(Number.NEGATIVE_INFINITY);
  const hi = new Float64Array(C).fill(Number.POSITIVE_INFINITY);
  const hasIn = new Uint8Array(C);
  const mask = new Uint8Array(C * 128);
  const alphaM = new Int32Array(R).fill(-1);
  const vals = new Int32Array(R).fill(-1);
  /** Sizes boxed once, so returning one allocates nothing. */
  const sizes: number[] = [];
  /** The same as a size theme gives them (`sizeOf`): integers. */
  const sized = new Int32Array(R).fill(NO_SIZE);
  const themes: (Theme | null)[] = [];
  let g: number = G.uniform;
  let c = 0;
  rules.forEach((rule, r) => {
    cs[r] = c;
    for (const cl of rule.when ?? []) {
      const fieldIn = cl.in != null && ms.stringField(cl.measure);
      if (cl.measure != null) {
        cm[c] = ms.add(cl.measure, fieldIn ? (cl.in ?? null) : null);
        g = coverGranularity(g, MEASURE_GRANULARITY[ms.code[cm[c]]]);
      }
      if (cl.min != null) lo[c] = cl.min;
      if (cl.max != null) hi[c] = cl.max;
      if (fieldIn) lo[c] = 1;
      else if (cl.in != null) {
        hasIn[c] = 1;
        g = coverGranularity(g, G.letter);
        for (const s of cl.in) {
          const code = s.length === 1 ? s.charCodeAt(0) : -1;
          if (code >= 0 && code < 128)
            mask[c * 128 + code] = mask[c * 128 + (letterOf(code) >= 0 ? code ^ 32 : code)] = 1;
        }
      }
      c++;
    }
    const v = (rule as Params)[kind];
    if (kind === "color") vals[r] = v == null ? NO_COLOR : (Color.decodeFirst(v) ?? NO_COLOR);
    else if (kind === "glyph") vals[r] = v == null ? NO_GLYPH : tc.glyphs.intern(v);
    sizes.push(typeof v === "number" ? v : Number.NaN);
    if (kind === "size") sized[r] = typeof v === "number" ? sizeOf(v) : NO_SIZE;
    const nested = rule.theme ? tc.resolve(rule.theme, channel) : null;
    themes.push(nested);
    if (nested) g = coverGranularity(g, nested.granularity);
    if (kind === "color" && rule.alpha != null) {
      alphaM[r] = ms.add(rule.alpha, null);
      g = coverGranularity(g, MEASURE_GRANULARITY[ms.code[alphaM[r]]]);
    }
  });
  cs[R] = c;
  const dflt = p.default ? tc.resolve(p.default as ThemeRef, channel) : null;
  if (dflt) g = coverGranularity(g, dflt.granularity);

  const code = Int32Array.from(ms.code);
  const { rows, data } = ms;
  const comp = ms.comp ?? null;
  const counts = comp?.counts ?? new Float64Array(0);
  const totals = comp?.totals ?? new Float64Array(0);
  const x = new Float64Array(1);
  /** A measure at the element into x[0]; NaN where it does not apply. */
  const read = (k: number, loc: ThemeLocation) => {
    const letter = loc.letter >= 0 && loc.letter < 128 ? letterOf(loc.letter) : -1;
    const col = loc.column;
    const op = code[k];
    const d = data[k];
    if (op === M.value) x[0] = loc.value;
    else if (op === M.field || op === M.fieldIn) {
      const at = loc.annotation >= 0 && loc.row >= 0 ? dataIndex(tc, loc.annotation) : -1;
      const v = at >= 0 ? rows[k][at] : null;
      x[0] = v && loc.row < v.length ? v[loc.row] : Number.NaN;
    } else if (op === M.column) x[0] = d && col >= 0 && col < d.length ? d[col] : Number.NaN;
    else if (op === M.equalsRow)
      x[0] = d && letter >= 0 && col >= 0 && col < d.length ? +(letterOf(d[col]) === letter) : Number.NaN;
    else if (letter < 0) x[0] = Number.NaN;
    else if (!comp) x[0] = op === M.rank ? 0 : op === M.share || (d as Uint8Array).indexOf(letter) >= 0 ? 1 : 0;
    else if (col < 0 || col >= totals.length || totals[col] === 0) x[0] = Number.NaN;
    else if (op === M.rank) {
      let n = 0;
      for (let j = 0; j < 26; j++) if (counts[col * 26 + j] > counts[col * 26 + letter]) n++;
      x[0] = n;
    } else {
      let n = 0;
      if (op === M.share) n = counts[col * 26 + letter];
      else for (let j = 0; j < (d as Uint8Array).length; j++) n += counts[col * 26 + (d as Uint8Array)[j]];
      x[0] = n / totals[col];
    }
  };
  const holds = (k: number, loc: ThemeLocation): boolean => {
    if (hasIn[k] === 1 && (loc.letter < 0 || loc.letter >= 128 || mask[k * 128 + loc.letter] === 0)) return false;
    const m = cm[k];
    if (m < 0) return true;
    read(m, loc);
    return x[0] >= lo[k] && x[0] <= hi[k];
  };
  /** The first rule that holds; -1: none. */
  const first = (loc: ThemeLocation): number => {
    for (let r = 0; r < R; r++) {
      let k = cs[r];
      while (k < cs[r + 1] && holds(k, loc)) k++;
      if (k === cs[r + 1]) return r;
    }
    return -1;
  };
  const fades = alphaM.some((k) => k >= 0) || themes.some((n) => n?.alpha) || !!dflt?.alpha;
  // A draw asks `value` and then `alpha` of one location: the rule found for it is kept, by what the measures read
  // (letter, column, row, annotation and value), so the rules are searched once (typed slots: no double is boxed).
  const memo = new Int32Array(5).fill(-2);
  const memoValue = new Float64Array(1);
  const search = !fades
    ? first
    : (loc: ThemeLocation): number => {
        if (
          memo[0] === loc.letter &&
          memo[1] === loc.column &&
          memo[2] === loc.row &&
          memo[3] === loc.annotation &&
          Object.is(memoValue[0], loc.value)
        )
          return memo[4];
        const r = first(loc);
        memo[0] = loc.letter;
        memo[1] = loc.column;
        memo[2] = loc.row;
        memo[3] = loc.annotation;
        memo[4] = r;
        memoValue[0] = loc.value;
        return r;
      };
  const none = kind === "color" ? NO_COLOR : kind === "glyph" ? NO_GLYPH : NO_SIZE;
  const value =
    kind === "size"
      ? (loc: ThemeLocation) => {
          const r = search(loc);
          const n = r >= 0 ? themes[r] : dflt;
          return n ? n.value(loc) : r >= 0 ? sized[r] : none;
        }
      : (loc: ThemeLocation) => {
          const r = search(loc);
          const n = r >= 0 ? themes[r] : dflt;
          return n ? n.value(loc) : r >= 0 ? vals[r] : none;
        };
  const alpha = (loc: ThemeLocation): number => {
    const r = search(loc);
    const n = r >= 0 ? themes[r] : dflt;
    if (r < 0 || alphaM[r] < 0) return n?.alpha ? n.alpha(loc) : n || r >= 0 ? OPAQUE : 0;
    read(alphaM[r], loc);
    const a = x[0];
    return a >= 1 ? OPAQUE : a > 0 ? Math.round(a * 255) : 0;
  };
  const labeled = rules.some((r) => r.label);
  return {
    value,
    granularity: g,
    alpha: kind === "color" && fades ? alpha : undefined,
    legend: labeled
      ? (out: Parameters<NonNullable<Theme["legend"]>>[0]) => {
          rules.forEach((rule, r) => {
            if (!rule.label) return;
            if (kind === "size")
              out.size(rule.label, sizes[r], channel.unit === "band" ? SizeUnit.band : SizeUnit.px, {
                channel: channel.name,
              });
            else if (kind === "glyph") out.item(rule.label, NO_COLOR, vals[r]);
            else if (vals[r] !== NO_COLOR) out.item(rule.label, vals[r], NO_GLYPH);
          });
        }
      : undefined,
  };
}

/** `rules` for the given kinds (the format's entry serves color, glyph and size). */
export function rulesTheme(name: string, kinds: readonly ChannelKind[]): ThemeProvider {
  return themeProvider(name, kinds, G.position, compile, {
    content: true,
    // `share` and `rank` need the columns held, or a profile.
    isApplicable: (tc, _kind, p) => {
      const aln = tc.set.alignments[alignmentOf(tc)];
      if (!aln?.lazy) return true;
      const rules = (p.rules ?? []) as readonly Rule[];
      // `equals_row` reads its row: a lazy alignment's only where it is held.
      const unheld = (m: unknown) => {
        if (!isObject(m) || !("equals_row" in m)) return false;
        const row = typeof m.equals_row === "number" ? m.equals_row : Math.max(0, aln.referenceRow ?? 0);
        return aln.rows.slotOf(row as never) < 0;
      };
      const reads = (m: unknown) => (p.profile == null && usesColumns(m)) || unheld(m);
      return !rules.some((r) => reads(r.alpha) || (r.when ?? []).some((cl) => reads(cl.measure)));
    },
  }) as ThemeProvider;
}
const RulesTheme = rulesTheme("rules", ["color", "glyph", "size"]);

/** Every built-in theme; `DefaultRegistry` registers them. */
export const DefaultThemes: readonly ThemeProvider[] = [...KernelThemes, RulesTheme, SchemeTheme];
