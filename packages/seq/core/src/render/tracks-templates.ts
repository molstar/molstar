/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { type Coverage, nodeKey, type Want } from '@molstar/seqx-core/data';
import { projectRanges } from '@molstar/seqx-core/resolve';
import type { ResolvedSection, ResolvedTrack, TrackTemplate } from '@molstar/seqx-core/view';
import { AxisTransform } from './axis.js';
import type { Overlay, StoreReader } from './layer-type.js';
import { boundedLevel } from './layer-data.js';
import type { LayoutInput, LayoutState } from './tracks-types.js';

// A section's track templates (_spec/README.md "Render", Shown nodes): what a template's source is asked for over the
// view, and whether an instance is in view (laid out) or not.

/**
 * A template source's region in view, on the section's system as the source names it, at the level that fits; and for
 * a template below its source, each nested source's systems' regions in view (through the overlay's path from the
 * section), so a nested source that sends its nodes by region is asked for them.
 */
export function wantTemplate(
  input: LayoutInput,
  data: StoreReader,
  section: ResolvedSection,
  template: TrackTemplate,
  axis: AxisTransform,
  overlay: Overlay,
  out: (w: Want) => void,
) {
  const pieces = axis.visible([]);
  const at = templateSource(input, data, section, template, axis);
  if (at)
    for (let i = 0; i + 1 < pieces.length; i += 2)
      out({
        source: template.source,
        system: at.system,
        start: Math.floor(pieces[i]),
        end: Math.ceil(pieces[i + 1]),
        level: at.level,
      });
  if (template.depth < 1) return;
  for (const nested of overlay.nested(new Map([[template.source, template.depth]]))) {
    const cov = data.coverage(nested.address);
    for (const s of cov?.systems ?? []) {
      const target = overlay.systemNamed(s.name, nested.parent);
      if (!target) continue;
      const path =
        target.address === section.axis.system
          ? undefined
          : overlay.path(section.axis.system, target.address, section.via).path;
      if (!path && target.address !== section.axis.system) continue;
      const ranges = HalfOpenRanges.normalize(path ? projectRanges(path, pieces) : pieces);
      const ratio = path?.step === 1 ? 3 : path?.step === -1 ? 1 / 3 : 1;
      // bounded as a layer's level is: a long nested system at the fit is one coarse tile, never 81,307 at level 1
      const level = boundedLevel(cov!.levels, ratio / axis.devicePx, (axis.u1 - axis.u0) * ratio, data);
      for (let i = 0; i + 1 < ranges.length; i += 2)
        out({
          source: nested.address,
          system: s.name,
          start: Math.floor(ranges[i]),
          end: Math.ceil(ranges[i + 1]),
          level,
        });
    }
  }
}

/**
 * A template source's system as it names the section's, its coverage, and the level the template wants of it: the one
 * that fits a device pixel, coarser while the view would want more than the store's tiles a want (`boundedLevel`).
 */
function templateSource(
  input: LayoutInput,
  data: StoreReader,
  section: ResolvedSection,
  template: TrackTemplate,
  axis: AxisTransform,
): { readonly coverage: Coverage; readonly system: string; readonly level: number } | undefined {
  const coverage = data.coverage(template.source);
  const name = input.model.systems.get(section.axis.system)?.name;
  const info = coverage?.systems.find((s) => s.name === name || (name !== undefined && s.aliases?.includes(name)));
  if (!coverage || !info) return undefined;
  return {
    coverage,
    system: info.name,
    level: boundedLevel(coverage.levels, 1 / axis.devicePx, axis.u1 - axis.u0, data),
  };
}

/**
 * Whether a template's instance is in view (_spec/README.md "Render", Shown nodes): its node's top node (the one of
 * the template's source it came under) is listed up front in that source's coverage, or carried by a loaded tile of it
 * over the view at any level whose tiles there are within the store's bound; while the template's tiles over the view
 * are pending (a zoom's new level), an instance shown last time stays. Out of view it is not laid out: its track goes
 * and its hold with it, and its node is evicted as any other once the store needs room. Where the template's source
 * does not name the section's system it cannot be told: shown.
 */
export function instanceInView(
  input: LayoutInput,
  data: StoreReader,
  section: ResolvedSection,
  template: TrackTemplate,
  axis: AxisTransform,
  overlay: Overlay,
  instance: ResolvedTrack,
  state: LayoutState,
): boolean {
  const at = templateSource(input, data, section, template, axis);
  let address: string | undefined;
  for (const l of instance.layers) if (l.data && !address) address = l.data.address;
  let root = address ? overlay.node(address) : undefined;
  while (root && root.depth > 0) root = overlay.node(root.source);
  if (!at || !root) return true;
  const top = root;
  if (at.coverage.nodes?.some((x) => x.key === top.key)) return true;
  const key = nodeKey(top.source, top.key);
  const pieces = axis.visible([]);
  let pending = false;
  for (const level of at.coverage.levels)
    for (let i = 0; i + 1 < pieces.length; i += 2) {
      const want = {
        source: template.source,
        system: at.system,
        start: Math.floor(pieces[i]),
        end: Math.ceil(pieces[i + 1]),
        level,
      };
      // bounded by what it reads: a level whose tiles over the view pass the store's bound is not read
      if (data.tooMany?.(want)) continue;
      for (const k of data.keysOf(want)) {
        const e = data.get(k);
        if (e?.status === 'loaded') {
          if (e.nodes.includes(key)) return remember(state, instance.address, true);
        } else if (level === at.level && e?.status !== 'empty' && e?.status !== 'failed') pending = true;
      }
    }
  return remember(state, instance.address, pending && state.instancesShown.has(instance.address));
}

/** Records whether an instance was shown, for the next frame's pending rule; returns it. */
function remember(state: LayoutState, address: string, shown: boolean): boolean {
  if (shown) state.instancesShown.add(address);
  else state.instancesShown.delete(address);
  return shown;
}
