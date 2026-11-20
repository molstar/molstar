/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, parseAddress } from '@molstar/core/repr/tree';
import type { Viewport as RenderViewport } from '@molstar/seqx-core/render';
import type { Axis, ResolvedSection, ResolvedView } from '@molstar/seqx-core/view';
import { axesEqual, axisKey, axisTotal, breakHolding, nearestOffsetOf, offsetOf, positionAt } from './axis.js';

// The view's scroll groups (_spec/README.md "ViewModel"): a section's state and the group whose range it shows, sections
// of one layout on equal axes grouped (EC-113), a detached follower its own group; groups made again for a new
// resolved view, each range carried across it (EC-111); a group's range as its layout reads it.

export interface Group {
  readonly key: string;
  readonly axis: Axis;
  /** Its sections' layout: a group holds one (scroll groups across layouts would read one range two ways). */
  readonly layout: string;
  start: number;
  end: number;
  /** px a position, where the layout reads one (text). */
  scale?: number;
  /** It shows the whole axis and stays so across resizes and updates. */
  fitted: boolean;
  /** The last width (px) it was placed at: a resize keeps its zoom against it. */
  placedWidth: number;
  readonly members: string[];
}

export interface SectionState {
  readonly section: ResolvedSection;
  /** Its group's key; undefined for an attached follower. */
  group?: string;
  width: number;
  top: number;
  /** A grid's own y window (`follow_y: none`): kept apart from x's group (_spec/README.md "ViewModel", Grids). */
  y?: { start: number; end: number };
}

/** The key prefix of a detached follower's own group (its address after it). */
export const DetachedPrefix = '#';

/**
 * The sections and groups of a new resolved view, each viewport kept by the rules of EC-111: a section keeps its
 * width and top; a detached follower stays detached while its axis is the same; the others are grouped by axis and
 * layout (`groupEqualAxes`, else each its own), each group's range carried from the group of its key or the
 * section's last.
 */
export function groupSections(
  view: ResolvedView | undefined,
  root: Node,
  oldSections: ReadonlyMap<string, SectionState>,
  oldGroups: ReadonlyMap<string, Group>,
  groupEqualAxes: boolean,
): { sections: Map<string, SectionState>; groups: Map<string, Group> } {
  const sections = new Map<string, SectionState>();
  const groups = new Map<string, Group>();
  for (const section of view?.sections ?? []) {
    const old = oldSections.get(section.address);
    const state: SectionState = { section, width: old?.width ?? 0, top: old?.top ?? 0 };
    sections.set(section.address, state);
    if (section.follow) {
      // a detached follower stays detached while its axis is the same
      const kept = old?.group?.startsWith(DetachedPrefix) ? oldGroups.get(old.group) : undefined;
      if (kept && axesEqual(kept.axis, section.axis)) {
        groups.set(kept.key, { ...kept, axis: section.axis });
        state.group = kept.key;
      }
      continue;
    }
    // a pending axis has no pieces: its key holds the section's own range, so only equal ranges group
    // groups only among sections of one layout: a range is read its layout's way (_spec/README.md "ViewModel")
    const axisPart = section.axis.ranges
      ? axisKey(section.axis)
      : `${axisKey(section.axis)}|${rangeParam(root, section.address)}`;
    const key = groupEqualAxes ? `${axisPart}|${section.layout}` : section.address;
    state.group = key;
    const existing = groups.get(key);
    if (existing) {
      existing.members.push(section.address);
      continue;
    }
    const same = oldGroups.get(key);
    const before = old?.group !== undefined ? oldGroups.get(old.group) : undefined;
    groups.set(key, carry(key, section, same && axesEqual(same.axis, section.axis) ? same : before));
  }
  return { sections, groups };
}

/**
 * A group's viewport across an update (EC-111): an equal axis keeps it; a changed one keeps its center and zoom when
 * the center is still shown, else fits; an axis that was pending (no length) keeps its intent where offsets are
 * positions from 0 on the same system. A new group fits.
 */
function carry(key: string, section: ResolvedSection, prev: Group | undefined): Group {
  const axis = section.axis;
  const total = axisTotal(axis) ?? 0;
  const group: Group = {
    key,
    axis,
    layout: section.layout,
    start: 0,
    end: total,
    fitted: true,
    placedWidth: 0,
    members: [section.address],
  };
  if (!prev) return group;
  group.placedWidth = prev.placedWidth;
  // the letter size is kept within one layout (a text across a patch of its params)
  if (prev.scale !== undefined && prev.layout === section.layout) group.scale = prev.scale;
  if (prev.fitted) return group;
  if (axesEqual(prev.axis, axis)) return { ...group, start: prev.start, end: prev.end, fitted: false };
  if (prev.axis.system !== axis.system) return group;
  if (!prev.axis.ranges) {
    // a pending axis's offsets are positions from 0: kept where the new axis starts there
    return axis.ranges?.[0] === 0 ? { ...group, start: prev.start, end: prev.end, fitted: false } : group;
  }
  const center = positionAt(prev.axis, (prev.start + prev.end) / 2);
  const c = center === undefined ? undefined : centerOffset(axis, center);
  if (c === undefined) return group;
  const half = (prev.end - prev.start) / 2;
  return { ...group, start: c - half, end: c + half, fitted: false };
}

/**
 * A centre position's offset on an axis: where it is shown; in a stretch a break hides (a collapse), the break's (the
 * nearest piece edge, zero wide); outside the axis, none.
 */
export function centerOffset(axis: Axis, center: number): number | undefined {
  const o = offsetOf(axis, center);
  if (o !== undefined) return o;
  return breakHolding(axis, Math.floor(center)) >= 0 ? nearestOffsetOf(axis, center) : undefined;
}

/** A group's range as its layout reads it. */
export function rangeOf(g: Group): RenderViewport {
  return g.scale !== undefined ? { start: g.start, end: g.end, scale: g.scale } : { start: g.start, end: g.end };
}

/** A section's own `range` param as written ('all' without one): part of a pending axis's group key. */
function rangeParam(root: Node, address: string): string {
  const p = parseAddress(address).pattern;
  let n: Node | undefined = root;
  for (const seg of p?.segments ?? []) n = typeof seg === 'string' ? n?.children?.find((c) => c.id === seg) : undefined;
  const r = n?.params?.range;
  return Array.isArray(r) ? r.join(',') : 'all';
}
