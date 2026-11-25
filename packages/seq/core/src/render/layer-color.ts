/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's color channel as its track compiles it (_spec/README.md "Themes"; "Render", Layer colors): the defaults the
// letters, values and ranges types fall back to, a colourless ranges layer's categorical place in its view, the view's
// theme context (its scales, held domains, the reader's substitution), and the one compile through the layer's kept
// channels. The layer types and `layerColors` read these, so a structure colored by a track agrees with the track.

import type { SequenceAlphabet } from '@molstar/core/util/sequence-constants';
import {
  type Channel,
  lettersDefault,
  type ScaleParams,
  type ThemeBody,
  type ThemeEnv,
} from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer, ResolvedTrack, ResolvedView } from '@molstar/seqx-core/view';
import {
  type ChannelContext,
  channelForms,
  DefaultThemeContext,
  type LayerChannels,
  type ThemeContext,
} from './channels.js';
import { layerKey } from './tracks-columns.js';
import { CategoryCount, type SeqTokens } from './tokens.js';

/** A letters layer's color default and environment by alphabet, made once (a layout makes none to read a kept channel). */
export const LettersColorFallback: Record<SequenceAlphabet, ThemeBody> = {
  protein: { kind: 'scheme', params: { name: lettersDefault('protein') } },
  dna: { kind: 'scheme', params: { name: lettersDefault('dna') } },
  rna: { kind: 'scheme', params: { name: lettersDefault('rna') } },
};
export const LettersEnv: Record<SequenceAlphabet, { alphabet: SequenceAlphabet }> = {
  protein: { alphabet: 'protein' },
  dna: { alphabet: 'dna' },
  rna: { alphabet: 'rna' },
};

/** A values layer's color default, and its environment: its one field is the value (made once). */
export const ValueColorFallback: ThemeBody = 'token:value';
export const ValueEnv = { numeric: (n: string) => n === 'value' };

/** The categorical palette's tokens, a colourless ranges layer's default by its place (`LayoutContext.colorIndex`). */
export const CategoryTokens = Array.from({ length: CategoryCount }, (_, k) => `token:cat-${k}`);

/** A ranges layer's color default: its categorical place, else `token:range`. */
export function rangesColorFallback(colorIndex: number | undefined): ThemeBody {
  return colorIndex === undefined ? 'token:range' : CategoryTokens[colorIndex % CategoryTokens.length];
}

/**
 * What a layer with no color is colored by (_spec/README.md "Themes", A layer with no color): its key, the data it
 * shows (a ring and its detail showing one node: one color) or, for a template's layer, the template's (its instances:
 * one color); undefined for a layer the spec colors or a type that has its own default.
 */
export function colorKey(track: ResolvedTrack, layer: ResolvedLayer): string | undefined {
  if (layer.theme.color !== undefined) return undefined;
  return placeKey(track, layer);
}

/** A ranges layer's key in the palette's places, coloured or not (its data, or a template's instances). */
function placeKey(track: ResolvedTrack, layer: ResolvedLayer): string | undefined {
  if (layer.type !== 'ranges') return undefined;
  return track.template ? `${track.address}:${layer.index}` : layer.data?.address;
}

/**
 * The colourless layers' places in the categorical palette, by their keys' first appearance in the view (a layer's
 * order), counted over every ranges layer, coloured or not: a layer's default never depends on another's scheme
 * (coloring one moved the next one's color, and its preview too). Made once a view.
 */
const DefaultColors = new WeakMap<ResolvedView, Map<string, number>>();
export function defaultColors(view: ResolvedView): Map<string, number> {
  let m = DefaultColors.get(view);
  if (m) return m;
  m = new Map();
  for (const s of view.sections)
    for (const t of s.tracks)
      for (const l of t.layers) {
        const k = placeKey(t, l);
        if (k !== undefined && !m.has(k)) m.set(k, m.size % CategoryCount);
      }
  DefaultColors.set(view, m);
  return m;
}

/** What a view's theme reads beside its scales: the registry's schemes and measures, and the reader's setting. */
export type ThemeInput = Partial<Pick<ThemeContext, 'schemes' | 'measures' | 'cvdSafe'>>;

/** A view's theme context: its scales, the held domains given, the registry's schemes, the reader's substitution. */
export function viewTheme(
  view: ResolvedView,
  input: ThemeInput | undefined,
  domain: (key: string) => readonly [number, number] | undefined,
  extent: ThemeContext['extent'] = () => {},
): ThemeContext {
  return {
    schemes: input?.schemes ?? DefaultThemeContext.schemes,
    ...(input?.measures ? { measures: input.measures } : {}),
    cvdSafe: input?.cvdSafe ?? false,
    scale: (a) => view.scales?.[a] as ScaleParams | undefined,
    domain,
    extent,
  };
}

/** The color channels of a view's layers, compiled as their tracks compile them. */
export interface LayerColors {
  /**
   * A letters, values or ranges layer's color channel with these tokens (its kept channels: a second call with the
   * same tokens reads them); `extra` is the type's environment beside the default (a ranges layer's `numeric`, a letters
   * layer's `alphabet`). Undefined for a layer of the view's no track holds, or of another type.
   */
  channel(layer: ResolvedLayer, tokens: SeqTokens, extra?: Partial<ThemeEnv>): Channel | undefined;
}

const NoReport = () => {};

/** The layer colors of a view under a theme context (a renderer's: `Renderer.layerColors()`, with its held domains). */
export function layerColors(view: ResolvedView, theme: ThemeContext): LayerColors {
  const kept = new WeakMap<ResolvedLayer, LayerChannels>();
  const trackOf = (layer: ResolvedLayer) => {
    for (const s of view.sections) for (const t of s.tracks) if (t.layers.includes(layer)) return t;
    return undefined;
  };
  return {
    channel(layer, tokens, extra = {}) {
      const track = trackOf(layer);
      const type = [layer.type, ...layer.fallbacks].find((t) => t === 'letters' || t === 'values' || t === 'ranges');
      if (!track || !type) return undefined;
      let channels = kept.get(layer);
      if (!channels) kept.set(layer, (channels = channelForms(layer.theme)));
      const ctx: ChannelContext = { tokens, theme, key: layerKey(track, layer), report: NoReport };
      if (type === 'letters') {
        const alphabet = extra.alphabet ?? 'protein';
        return channels.get('color', LettersColorFallback[alphabet], ctx, { ...LettersEnv[alphabet], ...extra });
      }
      if (type === 'values') return channels.get('color', ValueColorFallback, ctx, { ...ValueEnv, ...extra });
      const k = colorKey(track, layer);
      return channels.get(
        'color',
        rangesColorFallback(k === undefined ? undefined : defaultColors(view).get(k)),
        ctx,
        extra,
      );
    },
  };
}
