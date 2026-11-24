/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ResolvedSection, ResolvedTrack, ResolvedTrackGroup } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';
import {
  builtInHeaderStates,
  drawsSectionHeader,
  type HeaderInputFacts,
  type HeaderLine,
  hiddenTracksOf,
  type SectionHeaderLine,
  sectionHeaderLine,
  trackHeaderLine,
} from './headers.js';
import type { Box } from './layer-type.js';
import { layerKey } from './tracks-columns.js';
import type { RangesLay } from './ranges-lay.js';
import type { LayerLayout, LayoutInput } from './tracks-types.js';

// The tracks layout's headers (_spec/README.md "Render", Track headers): each header's facts read from the view and
// this layout (the rows a filter left out, the hidden layers and tracks), its parts from them (`LayoutInput.headers`,
// the built-in states when absent), and its line laid out.

const NoFiltered: HeaderInputFacts['filtered'] = [];
const NoKeys: readonly string[] = [];
const HiddenLayers = new WeakMap<ResolvedTrack, readonly string[]>();

/** A track's layers the spec hides (`hidden: true`), by key; kept a resolved track. */
function hiddenLayersOf(track: ResolvedTrack): readonly string[] {
  let out = HiddenLayers.get(track);
  if (out) return out;
  let found: string[] | undefined;
  for (const l of track.layers) if (l.params.hidden === true) (found ??= []).push(layerKey(track, l));
  out = found ?? NoKeys;
  HiddenLayers.set(track, out);
  return out;
}

/** A section's header line when it draws one (`section.header`), from `top`; its height is `header.sectionPx`. */
export function sectionHeaderOf(
  input: LayoutInput,
  section: ResolvedSection,
  id: string,
  top: number,
  config: RenderConfig,
): SectionHeaderLine | undefined {
  if (input.picture) return undefined;
  const facts: HeaderInputFacts = {
    section,
    filtered: NoFiltered,
    hiddenLayers: NoKeys,
    hiddenTracks: hiddenTracksOf(input.model, section),
  };
  const parts = (input.headers ?? builtInHeaderStates)(facts);
  if (!drawsSectionHeader(parts, config)) return undefined;
  return sectionHeaderLine(parts, section, id, top, config, input.text, input.figure === true);
}

/** What a facts record spreads for a track outside a folded group: one object, none allocated a header. */
const NoGroup = Object.freeze({});

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
  group?: ResolvedTrackGroup,
): { readonly header: HeaderLine; readonly labelY: number } | undefined {
  // a navigator's picture has no label lane
  if (input.picture) return undefined;
  // its layers', then its `section` layers' (pairs: laid out once a section)
  let filtered: { layer: string; count: number }[] | undefined;
  for (let i = 0; i < layers.length; i++) {
    const l = layers[i];
    const count = l.type.filtered?.(l.lay) ?? 0;
    if (count) (filtered ??= []).push({ layer: l.key, count });
  }
  let scrolls = false;
  for (let i = 0; i < layers.length && !scrolls; i++) {
    const l = layers[i];
    scrolls = l.type.name === 'ranges' && (l.lay as RangesLay).scrolls && (l.lay as RangesLay).scrollMax > 0;
  }
  for (let i = 0; i < sectionLays.length; i++) {
    const l = sectionLays[i];
    if (!track.layers.includes(l.layer)) continue;
    const count = l.type.filtered?.(l.lay) ?? 0;
    if (count) (filtered ??= []).push({ layer: layerKey(track, l.layer), count });
  }
  const facts: HeaderInputFacts = {
    section,
    track,
    filtered: filtered ?? NoFiltered,
    hiddenLayers: hiddenLayersOf(track),
    hiddenTracks: hiddenTracksOf(input.model, section),
    scrolls: scrolls || undefined,
    ...(group ? { group } : NoGroup),
  };
  // a rows track's label column is its tree and its rows' labels: its header waits for its title row (the package's
  // TODO): no parts and no label hit over them
  const rows = track.layers.some((l) => l.type === 'rows');
  const parts = rows ? undefined : (input.headers ?? builtInHeaderStates)(facts);
  // a folded group's row: its line the group's (its slots and its menu act on the group, Decision 4 of the groups)
  const out = trackHeaderLine(
    parts,
    group?.address ?? track.address,
    id,
    well,
    config,
    input.text,
    input.figure === true,
    !rows,
    // the label the states beside it leave whole (a folded group's row: the group's own, drawn elsewhere)
    group || rows
      ? 0
      : input.text.measure(config.labelFont, track.label) +
          (track.group && !track.summaryOf ? config.header.groupIndentPx : 0),
  );
  return group ? { header: { ...out.header, group: group.address }, labelY: out.labelY } : out;
}
