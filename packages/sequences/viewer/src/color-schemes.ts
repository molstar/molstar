/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Color schemes in the viewer (seqspec:interaction#color-schemes-at-runtime): the `scheme` theme (the theme its scheme
// stands for, the skin's variant, made through `tc.resolve`: nothing of its own per letter), the runtime override per
// layer (`ui.colorSchemes`, put in place in the graph's theme stage), and the legend's picker both UIs draw.

import {
  Color,
  type ColorListName,
  ColorLists,
  type ColorSchemeEntry,
  type DiscretePalette,
  PaletteDefaults,
  type SchemeParams,
  type SchemeSubject,
  type SvsRegistries,
  schemeFits,
  schemeLegend,
  schemeOf,
  schemeSubject,
  schemesFor,
  schemeTheme,
  type ThemeRef,
} from "@molstar/seq-core/format";
import type { Alphabet } from "@molstar/seq-core/model";
import type { SeqViewerContext } from "./context.js";
import { NO_COLOR } from "./draw/backend.js";
import { Granularity, NO_GLYPH, ofAlphabet, type Theme, type ThemeContext, type ThemeProvider } from "./themes.js";

/** The layer's color subject (letters, values, categories) in its theme context; null: no scheme fits it. */
export function subjectOf(tc: ThemeContext): SchemeSubject | null {
  const type = tc.layer.type >= 0 ? tc.svs?.get("layerTypes", tc.strings.get(tc.layer.type)) : undefined;
  const of = ofAlphabet(tc);
  return schemeSubject(
    tc.set,
    tc.section,
    type,
    tc.annotations,
    of === undefined ? undefined : (of as Alphabet | null),
  );
}

/**
 * `scheme` (seqspec:format#themes): its scheme's theme for the skin, created as a base is, with the scheme's labeled
 * groups listed in the legend before the theme's own. Applies where the scheme fits the layer's subject; else its
 * fallback.
 */
export const SchemeTheme: ThemeProvider = {
  name: "scheme",
  label: "scheme",
  category: "kernel",
  kinds: ["color"],
  granularity: Granularity.letter,
  content: true,
  isApplicable(tc, _kind, p) {
    const def = schemeOf(tc.svs, p as SchemeParams);
    const s = subjectOf(tc);
    return !!def && !!s && schemeFits(def, s);
  },
  create(tc, p, channel): Theme {
    const ref = schemeTheme(tc.svs, p as SchemeParams, { dark: tc.dark });
    const inner = ref && tc.resolve(ref, channel);
    if (!inner)
      return { kind: "color", granularity: 0, reduction: "classes", value: () => NO_COLOR, contextHash: tc.stamp() };
    const items = schemeLegend(tc.svs, p as SchemeParams, tc.dark).map(
      (x) => [x.label, Color.decodePrimary(x.color)] as const,
    );
    if (items.length === 0) return inner;
    const own = inner.legend;
    return {
      ...inner,
      legend: (out) => {
        for (const [label, color] of items) out.item(label, color, NO_GLYPH);
        own?.(out);
      },
    };
  },
};

/** The color an override puts in place of the layer's own: the scheme, with the own theme's field and domain. */
export function overrideTheme(own: ThemeRef | null, name: string): ThemeRef {
  const p = (own?.params ?? {}) as Record<string, unknown>;
  const given =
    own?.theme === "scheme" ? p.domain : (p.palette as { value_domain?: unknown } | undefined)?.value_domain;
  const domain = Array.isArray(given) && given.some((x) => x != null) ? given : null;
  return { theme: "scheme", params: { name, field: typeof p.field === "string" ? p.field : null, domain } };
}
/**
 * A scheme picked at runtime over the layer's own color: the override, falling back to the own theme, so a pick kept
 * across a patch after which it no longer fits draws the layer's own color.
 */
export const pickedTheme = (own: ThemeRef | null, name: string): ThemeRef => ({
  ...overrideTheme(own, name),
  fallback: own,
});

/**
 * A color theme over a range set at runtime (the legend's editor): a `scheme`'s `domain`, or a normalized continuous
 * or discrete `field_palette`'s `value_domain`, with the ends' colors beyond it (the range may be wider or narrower
 * than the data: a continuous palette's flow colors, a discrete one's first and last bins open-ended); any other theme,
 * an absolute palette's included (its checkpoints are its range), is as it was.
 */
export function withRange(ref: ThemeRef, range: readonly [number, number]): ThemeRef {
  const p = (ref.params ?? {}) as Record<string, unknown>;
  const palette = p.palette as { kind?: string; mode?: string } | undefined;
  if (ref.theme === "scheme") return { ...ref, params: { ...p, domain: range } };
  if (ref.theme !== "field_palette" || !palette || palette.kind === "categorical" || palette.mode === "absolute")
    return ref;
  if (palette.kind === "discrete")
    return { ...ref, params: { ...p, palette: { ...openEnded(palette as DiscretePalette), value_domain: range } } };
  const flow = { underflow_color: "auto", overflow_color: "auto" };
  return { ...ref, params: { ...p, palette: { ...palette, ...flow, value_domain: range } } };
}

/** A discrete palette's bins as `[color, from, to]` ranges, its first open below and its last above. */
function openEnded(palette: DiscretePalette): DiscretePalette {
  const given = palette.colors ?? PaletteDefaults.discrete.colors;
  const list = (typeof given === "string" ? (ColorLists[given as ColorListName] ?? []) : given) as readonly unknown[];
  if (list.length === 0) return palette;
  let bins: [unknown, number | null, number | null][];
  if (list.every((c) => !Array.isArray(c))) {
    const colors = palette.reverse ? list.slice().reverse() : list;
    bins = colors.map((c, i) => [c, i / colors.length, (i + 1) / colors.length]);
  } else if (list.every((c) => Array.isArray(c) && c.length === 2)) {
    const pairs = list as readonly (readonly [unknown, number])[];
    bins = pairs.map(([c, from], i) => [c, from, pairs[i + 1]?.[1] ?? null]);
  } else bins = (list as readonly [unknown, number | null, number | null][]).map((b) => [...b]);
  bins[0][1] = null;
  bins[bins.length - 1][2] = null;
  return { ...palette, colors: bins as never, reverse: false };
}

/** What a theme node's key holds for its schemes: the entries named (nested in `rules` too), and the override's. */
export function schemeEntries(svs: SvsRegistries, refs: readonly (ThemeRef | null)[], override?: string): unknown[] {
  const out: unknown[] = [];
  const walk = (r: ThemeRef | null | undefined) => {
    if (!r) return;
    const p = r.params as
      | { name?: unknown; rules?: { theme?: ThemeRef | null }[]; default?: ThemeRef | null }
      | undefined;
    if (r.theme === "scheme" && typeof p?.name === "string") out.push(svs.get("colorSchemes", p.name));
    for (const x of p?.rules ?? []) walk(x.theme);
    walk(p?.default);
  };
  for (const r of refs) walk(r);
  if (override !== undefined) out.push(override, svs.get("colorSchemes", override));
  return out;
}

/** The schemes that fit a layer's color, by its subject (`PreparedLayer.subject`); empty: none, or no subject. */
export function layerSchemes(ctx: SeqViewerContext, layer: string): readonly ColorSchemeEntry[] {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return [];
  const i = Array.from(view.resolved.layers.ref).indexOf(view.set.strings.find(layer));
  const subject = i >= 0 ? view.prepared[i]?.subject : null;
  return subject ? schemesFor(ctx.registry.svs, subject) : [];
}

export interface ColorSchemePicker {
  /** The layer it switches. */
  readonly layer: string;
  readonly name: string;
  /** The scheme shown: the override, else the document's scheme name; null: the document's own color. */
  readonly value: string | null;
  readonly options: readonly { readonly value: string | null; readonly label: string }[];
}
/**
 * A picker for each of the track's layers whose color more than one scheme fits (seqspec:interaction#the-legend):
 * "document" (back to the document's own), then the schemes, colorblind-safe ones marked. Edge-rate: once per assembled
 * view.
 */
export function colorSchemePickers(ctx: SeqViewerContext, track: number): ColorSchemePicker[] {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return [];
  const { resolved: rv, set, prepared } = view;
  const out: ColorSchemePicker[] = [];
  for (let i = 0; i < rv.layers.count; i++) {
    const subject = prepared[i]?.subject;
    if (rv.layers.track[i] !== track || !subject) continue;
    const schemes = schemesFor(ctx.registry.svs, subject);
    if (schemes.length < 2) continue;
    const layer = set.strings.get(rv.layers.ref[i]);
    const own = rv.layers.themes[i][0];
    const named = own?.theme === "scheme" ? (own.params as { name?: unknown } | undefined)?.name : null;
    const value = ctx.ui.colorSchemes.value.get(layer) ?? (typeof named === "string" ? named : null);
    const label = (e: ColorSchemeEntry) => (e.colorblind ? `${e.label} (colorblind-safe)` : e.label);
    out.push({
      layer,
      name: "Colors",
      value,
      options: [{ value: null, label: "document" }, ...schemes.map((e) => ({ value: e.name, label: label(e) }))],
    });
  }
  return out;
}
