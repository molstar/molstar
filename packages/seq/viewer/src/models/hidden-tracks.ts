/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import type { Node } from '@molstar/core/repr/tree';
import { type HiddenTrack, hiddenTrackName } from '@molstar/seqx-core/render';
import { idRunOf } from '@molstar/seqx-core/lang/check';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { tracksOf } from '@molstar/seqx-core/view';

/** A track the spec hides, with the section it is in (`SeqViewer.hiddenTracks`). */
export interface HiddenTrackEntry extends HiddenTrack {
  /** The section's address. */
  readonly section: string;
  /** A hidden group of tracks is one entry ("Variants (3 tracks)"): its `track.show` shows the group. */
  readonly kind?: 'group';
  /** The group a hidden track sits in (shown), its address and label. */
  readonly group?: { readonly address: string; readonly label: string };
}

const None: readonly HiddenTrackEntry[] = [];

function sameEntries(a: readonly HiddenTrackEntry[], b: readonly HiddenTrackEntry[]): boolean {
  return (
    a === b ||
    (a.length === b.length &&
      a.every((x, i) => x.address === b[i].address && x.label === b[i].label && x.section === b[i].section))
  );
}

/** Every hidden track (`hidden: true`) of every section of a spec, in its order. */
export function hiddenTracksOf(root: Node): readonly HiddenTrackEntry[] {
  let out: HiddenTrackEntry[] | undefined;
  const walk = (node: Node, path: readonly string[]) => {
    for (const c of node.children ?? []) {
      if (c.id === undefined) continue;
      const at = [...path, c.id];
      if (c.kind === 'section') {
        const section = `/${at.join('/')}`;
        const groups = new Set<string>();
        for (const { node: t, address, group } of tracksOf(c, section)) {
          if (group && group.node.params?.hidden === true) {
            // a hidden group: one entry, its tracks counted
            if (groups.has(group.address)) continue;
            groups.add(group.address);
            const n = (group.node.children ?? []).filter((x) => x.kind === 'track').length;
            const label = typeof group.node.params?.label === 'string' ? group.node.params.label : group.node.id!;
            (out ??= []).push({
              address: group.address,
              label: `${label} (${n} track${n === 1 ? '' : 's'})`,
              section,
              kind: 'group',
            });
            continue;
          }
          if (t.params?.hidden !== true) continue;
          const g = group && {
            address: group.address,
            label: typeof group.node.params?.label === 'string' ? group.node.params.label : group.node.id!,
          };
          (out ??= []).push({ address, label: hiddenTrackName(t), section, ...(g ? { group: g } : {}) });
        }
        continue;
      }
      if (c.kind !== 'track') walk(c, at);
    }
  };
  walk(root, []);
  return out ?? None;
}

/**
 * The hidden tracks of the spec shown (_spec/README.md "Hidden tracks"): the application's to list and restore, never
 * drawn by the view. It reads the spec, not the layout: a section whose tracks are all hidden has none to read there.
 */
export function hiddenTracksState(model: ReadonlyState<StaticModel>): ReadonlyState<readonly HiddenTrackEntry[]> {
  return derive([model], (m) => hiddenTracksOf(m.spec), sameEntries);
}

/**
 * Whether what data owners make (an item's key) is drawn only on hidden tracks: every layer of the spec whose data is
 * (a form over) one of the owners lies on a hidden track (or in a hidden group), and there is at least one. A key also on a shown track is
 * not (a jump goes there).
 */
export function onlyOnHiddenTracks(root: Node, owners: ReadonlySet<string>): boolean {
  let hidden = 0;
  let shown = 0;
  const walk = (node: Node, inHidden: boolean) => {
    for (const c of node.children ?? []) {
      if (c.kind === 'layer') {
        const data = c.params?.data;
        if (typeof data === 'string' && owners.has(`/${idRunOf(data).join('/')}`)) {
          if (inHidden) hidden++;
          else shown++;
        }
      } else walk(c, inHidden || ((c.kind === 'track' || c.kind === 'group') && c.params?.hidden === true));
    }
  };
  walk(root, false);
  return hidden > 0 && shown === 0;
}
