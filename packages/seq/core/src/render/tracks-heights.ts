/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { ResolvedLayer, ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';
import { RenderIssues } from './issues.js';
import type { LayerType } from './layer-type.js';
import {
  BuiltInSizings,
  ContentSizing,
  heightsFit,
  knownVersion,
  lanesHeight,
  type LayerContent,
  type SizedTrack,
} from './sizing.js';
import { TracksLayout } from './layout/tracks.js';
import type { SectionLayoutType } from './layout/types.js';
import { RowHeight } from './rows-model.js';
import type { LayoutInput, LayoutState, SizedHeights } from './tracks-types.js';
import { ownPieces } from './tracks-axes.js';
import { titleRowOf } from './tracks-columns.js';

// A section's tracks' heights (_spec/README.md "Render", Track heights): the sizing over the tracks' contents, a
// track's own height and the host's port, kept while their key is; a rows track's height and the ruler a rail sits on.

/** Two lists of the same tracks (by identity). */
function sameTracks(a: readonly ResolvedTrack[], b: readonly ResolvedTrack[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/**
 * A section's tracks' heights (_spec/README.md "Render", Track heights): the section's `heights` sizing over its
 * tracks' content heights (a layout with heights of its own, text, keeps them), a track's own height over it, the
 * host's port over all. Kept in the state by the section (and `lead`, its leader's reading) while its key is: its
 * shown tracks, the px a position of its whole extent, the config, the layout, the metrics' version, the registry and
 * the port; so a pan runs neither the sizing nor the port. Each track's content (one lane block a track) is kept with
 * it.
 */
export function bandHeights(
  input: LayoutInput,
  state: LayoutState,
  config: RenderConfig,
  types: Readonly<Record<string, LayerType<any, any>>>,
  dataW: number,
  section: ResolvedSection,
  tracks: readonly ResolvedTrack[],
  layout: SectionLayoutType<any>,
  lead: boolean,
  say: (i: Issue) => void,
): SizedHeights {
  // the px a position with the section's own whole extent across the data box, whatever its window shows (a pan over
  // a cyclic origin, a follower's range, a zoom): a label fitting there fits at every zoom in
  const own = layout.heights ? undefined : ownPieces(section, input.model, input.data);
  let total = 0;
  if (own) for (let k = 0; k + 1 < own.length; k += 2) total += Math.abs(own[k + 1] - own[k]);
  const fit = total > 0 ? dataW / total : 0;
  const memo = lead ? state.leadSizes : state.sizes;
  const kept = memo.get(section.address);
  if (
    kept &&
    kept.fit === fit &&
    kept.config === config &&
    kept.layout === layout &&
    kept.version === input.text.version &&
    kept.sizings === input.sizings &&
    kept.hook === input.trackSizing &&
    kept.known === knownVersion(types) &&
    sameTracks(kept.tracks, tracks)
  ) {
    if (kept.issue) say(kept.issue);
    return kept;
  }
  const made = sizeTracks(input, state, config, types, section, tracks, layout, fit);
  memo.set(section.address, made);
  if (made.issue) say(made.issue);
  return made;
}

/** The sizing pass itself (`bandHeights` keeps what it gives). */
function sizeTracks(
  input: LayoutInput,
  state: LayoutState,
  config: RenderConfig,
  types: Readonly<Record<string, LayerType<any, any>>>,
  section: ResolvedSection,
  tracks: readonly ResolvedTrack[],
  layout: SectionLayoutType<any>,
  fit: number,
): SizedHeights {
  const byType = layout.heights?.(config) ?? config.heights;
  const fixed = tracks.map((t) => trackHeight(t, config, byType));
  const base = {
    tracks: [...tracks],
    fit,
    config,
    layout,
    version: input.text.version,
    sizings: input.sizings,
    hook: input.trackSizing,
    known: knownVersion(types),
  };
  if (layout.heights || !tracks.length) return { ...base, heights: fixed, contents: [], empties: [] };
  const cap = layout.name === 'rings' ? config.rings.trackMaxPx : config.track.maxPx;
  const cctx = { config, text: input.text, fit, labelLane: layout === TracksLayout && !input.picture };
  const sized: SizedTrack[] = [];
  const contents: (LayerContent | undefined)[] = [];
  // the tracks with no items (the tracks layout's: a ring's or a text's track is its own shape), by their `empty`
  const empties: ('note' | 'collapse' | undefined)[] = [];
  for (let j = 0; j < tracks.length; j++) {
    const t = tracks[j];
    let lanes: number | undefined;
    let labels = false;
    let partial = false;
    let other = 0;
    let kind = '';
    for (const l of t.layers) {
      // a hidden layer is no part of the track's kind or height
      if (l.params.hidden === true) continue;
      kind = kind ? `${kind}+${l.type}` : l.type;
      const type = l.placed ? firstType(l, types) : undefined;
      let c: LayerContent | undefined;
      if (type?.content) {
        let prep = state.preps.get(l);
        if (prep === undefined) state.preps.set(l, (prep = type.prepare(l, input.model)));
        c = type.content(prep, cctx);
      }
      if (c) {
        lanes = Math.max(lanes ?? 0, c.lanes);
        labels ||= c.labels;
        partial ||= c.partial === true;
      } else other = Math.max(other, byType[l.type as keyof RenderConfig['heights']] ?? 0);
    }
    // one lane block a track: the most lanes its layers need, one label lane for all of them
    // no items in what it covers: its content known whole (not a source's partial lanes: pending is not empty), no lane
    // of any layer, and no layer of a fixed height beside
    const mode = t.empty ?? config.track.empty;
    empties.push(
      layout === TracksLayout && lanes === 0 && !partial && other === 0 && mode !== 'keep' && t.height === undefined
        ? mode
        : undefined,
    );
    contents.push(lanes === undefined ? undefined : partial ? { lanes, labels, partial } : { lanes, labels });
    // its lanes' px to the cap (`track.maxPx`; a ring's radial width `rings.trackMaxPx`), past which they scroll inside
    // the track (or compact, on a ring) with the label lane above them kept
    const lanesPx =
      lanes !== undefined ? Math.min(Math.max(other, lanesHeight(lanes, config)), Math.max(cap, other)) : 0;
    const content =
      lanes !== undefined && rowsFill(t) === undefined ? lanesPx + (labels ? config.laneHeight : 0) : undefined;
    sized.push({
      track: t,
      kind,
      fixed: fixed[j],
      // its label lane its own: what it adds over its lanes (a kind shares lanes, never a label lane)
      ...(content !== undefined ? { content, label: labels ? content - lanesPx : 0 } : {}),
      ...(t.height !== undefined ? { height: t.height } : {}),
    });
  }
  const name = section.params.heights;
  let sizing = typeof name === 'string' ? (BuiltInSizings[name] ?? input.sizings?.[name]) : ContentSizing;
  let issue: Issue | undefined;
  if (!sizing) {
    issue = issueOf(
      RenderIssues,
      'seq.render.heights',
      `The track sizing ${name} is not registered here: tracks fit their content`,
      { address: section.address, hint: RenderIssues.codes['seq.render.heights'].hint },
    );
    sizing = ContentSizing;
  }
  let named = sizing.heights(sized, config);
  if (!heightsFit(named, sized.length)) named = ContentSizing.heights(sized, config);
  const ownHeights = named.map((h, j) => sized[j].height ?? (empties[j] ? config.track.emptyPx : h));
  const hook = input.trackSizing?.(section, sized, ownHeights);
  const heights = heightsFit(hook, sized.length) ? [...hook] : ownHeights;
  return { ...base, heights, contents, empties, ...(issue ? { issue } : {}) };
}

/** The first shown track drawing a ruler (with height: a hidden ruler is no rail's), or −1. */
export function rulerTrack(shown: readonly ResolvedTrack[], heights: readonly number[]): number {
  for (let j = 0; j < shown.length; j++)
    if (heights[j] > 0 && shown[j].layers.some((l) => l.type === 'ruler' && l.placed && l.params.hidden !== true))
      return j;
  return -1;
}

/** How a track's rows layer takes its height: `fill`, a number of rows, or none (no rows layer). */
export function rowsFill(track: ResolvedTrack): 'fill' | number | undefined {
  for (const l of track.layers) {
    if (l.type !== 'rows') continue;
    const h = l.params.height;
    return typeof h === 'number' ? h : 'fill';
  }
  return undefined;
}

function trackHeight(track: ResolvedTrack, config: RenderConfig, heights: RenderConfig['heights']): number {
  if (track.height !== undefined) return track.height;
  const rows = rowsFill(track);
  if (typeof rows === 'number') return rows * RowHeight + titleRowOf(track, config);
  let h = 0;
  for (const l of track.layers) h = Math.max(h, heights[l.type as keyof RenderConfig['heights']] ?? 0);
  return h || heights.ranges;
}

/** A layer's type: the first of its fallbacks drawn here (a loop: no array made a layer a frame). */
export function firstType(
  layer: ResolvedLayer,
  types: Readonly<Record<string, LayerType<any, any>>>,
): LayerType<any, any> | undefined {
  const f = layer.fallbacks;
  for (let i = 0; i < f.length; i++) {
    const t = types[f[i]];
    if (t !== undefined) return t;
  }
  return undefined;
}
