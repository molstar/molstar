/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's tracks in the spec through its groups (_spec/README.md "Resolved view", Groups): the one walker every
// reader of a section's tracks takes, so a track in a group is read as the section's in spec order.

import type { Node } from '@molstar/core/repr/tree';
import { SeqLanguage } from '@molstar/seqx-core/lang';

/** A track of a section in the spec, with the group holding it (its node and address) when it sits in one. */
export interface SpecTrack {
  readonly node: Node;
  /** Its address (`<section>/<track>`, or `<section>/<group>/<track>`). */
  readonly address: string;
  readonly group?: { readonly node: Node; readonly address: string };
}

/** A group of a section's tracks as the resolved view says it (`ResolvedSection.trackGroups`). */
export interface ResolvedTrackGroup {
  readonly address: string;
  readonly label: string;
  /** The opening state (`collapsed`, false when absent): the viewer's folds are session state after it. */
  readonly collapsed: boolean;
  /** Its authored summary track's address (`summary` naming one of its tracks); absent for `auto` and `none`. */
  readonly summary?: string;
  /** `none`: a folded group draws its header line alone. */
  readonly none?: boolean;
  /** Its tracks shown open (its summary track not counted). */
  readonly count: number;
}

const Kinds = SeqLanguage.kinds as Readonly<Record<string, unknown>>;

/**
 * A section's tracks in spec order, a group's in place (depth first), hidden ones included. A child of a kind the
 * language lacks (a newer minor's container, as `group` is to an older reader: `FlattenedKinds`) is read through.
 */
export function tracksOf(section: Node | undefined, address: string): SpecTrack[] {
  const out: SpecTrack[] = [];
  for (const c of section?.children ?? []) {
    if (c.id === undefined) continue;
    if (c.kind === 'track') out.push({ node: c, address: `${address}/${c.id}` });
    else if (c.kind === 'group' || !Object.hasOwn(Kinds, c.kind)) {
      const g = { node: c, address: `${address}/${c.id}` };
      for (const t of c.children ?? [])
        if (t.kind === 'track' && t.id !== undefined)
          out.push({ node: t, address: `${g.address}/${t.id}`, ...(c.kind === 'group' ? { group: g } : {}) });
    }
  }
  return out;
}

/** The `summary` param of a group: `auto`, `none`, or the id of one of its tracks. */
export function groupSummary(group: Node): string {
  const s = group.params?.summary;
  return typeof s === 'string' ? s : 'auto';
}

/** Whether a track is drawn at all: neither it nor its group hidden. */
export function specTrackShown(t: SpecTrack): boolean {
  return t.node.params?.hidden !== true && t.group?.node.params?.hidden !== true;
}

/** The resolved groups from what the resolve loop counted (a group with no track shown left out: it draws nothing). */
export function resolvedTrackGroups(
  groups: ReadonlyMap<string, { readonly node: Node; readonly count: number; readonly summary?: string }>,
): ResolvedTrackGroup[] {
  const out: ResolvedTrackGroup[] = [];
  for (const [address, g] of groups) {
    if (g.count === 0) continue;
    const label = g.node.params?.label;
    out.push({
      address,
      label: typeof label === 'string' ? label : (g.node.id ?? ''),
      collapsed: g.node.params?.collapsed === true,
      ...(g.summary ? { summary: g.summary } : {}),
      ...(groupSummary(g.node) === 'none' ? { none: true } : {}),
      count: g.count,
    });
  }
  return out;
}

/** Whether a section holds a group of tracks (what the groups' lazy parts load for, `withViewKind`). */
export function holdsGroup(section: Node): boolean {
  const c = section.children;
  if (c) for (let i = 0; i < c.length; i++) if (c[i].kind === 'group') return true;
  return false;
}
