/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';
import {
  builtInHeaderStates,
  drawsSectionHeader,
  type HeaderFacts,
  type HeaderLine,
  hiddenLayersOf,
  hiddenTracksOf,
  type SectionHeaderLine,
  sectionHeaderLine,
  trackHeaderLine,
} from './headers.js';
import type { Box } from './layer-type.js';
import { layerKey } from './tracks-columns.js';
import type { LayerLayout, LayoutInput } from './tracks-types.js';

// The tracks layout's headers (_spec/README.md "Render", Track headers): each header's facts read from the view and
// this layout (the rows a filter left out, the hidden layers and tracks), its parts from them (`LayoutInput.headers`,
// the built-in states when absent), and its line laid out.

const NoFiltered: HeaderFacts['filtered'] = [];
const NoKeys: readonly string[] = [];

/** A section's header line when it draws one (`section.header`), from `top`; its height is `header.sectionPx`. */
export function sectionHeaderOf(
  input: LayoutInput,
  section: ResolvedSection,
  id: string,
  top: number,
  config: RenderConfig,
): SectionHeaderLine | undefined {
  if (input.picture) return undefined;
  const facts: HeaderFacts = {
    section,
    filtered: NoFiltered,
    hiddenLayers: NoKeys,
    hiddenTracks: hiddenTracksOf(input.model, section),
  };
  const parts = (input.headers ?? builtInHeaderStates)(facts);
  if (!drawsSectionHeader(parts, config)) return undefined;
  return sectionHeaderLine(parts, section, id, top, config, input.text, input.figure === true);
}

/**
 * A track's header line over its well: its facts (its layers' filtered rows, then its `section` layers', the first
 * layer with some first; its hidden layers), its parts, its line; and its label's baseline.
 */
export function trackHeaderOf(
  input: LayoutInput,
  section: ResolvedSection,
  id: string,
  track: ResolvedTrack,
  layers: readonly LayerLayout[],
  sectionLays: readonly LayerLayout[],
  well: Box,
  config: RenderConfig,
): { readonly header: HeaderLine; readonly labelY: number } {
  let filtered: { layer: string; count: number }[] | undefined;
  for (const own of [layers, sectionLays])
    for (const l of own) {
      if (own === sectionLays && !track.layers.includes(l.layer)) continue;
      const count = l.type.filtered?.(l.lay) ?? 0;
      if (count) (filtered ??= []).push({ layer: own === sectionLays ? layerKey(track, l.layer) : l.key, count });
    }
  const facts: HeaderFacts = {
    section,
    track,
    filtered: filtered ?? NoFiltered,
    hiddenLayers: hiddenLayersOf(track, (l) => layerKey(track, l)),
    hiddenTracks: hiddenTracksOf(input.model, section),
  };
  let parts = (input.headers ?? builtInHeaderStates)(facts);
  // a rows track's label column is its tree and its rows' labels: its slots wait for its title row (the package's TODO)
  if (parts && (parts.slots || parts.menu) && track.layers.some((l) => l.type === 'rows'))
    parts = parts.states ? { states: parts.states } : undefined;
  return trackHeaderLine(parts, track.address, id, well, config, input.text, input.figure === true);
}
