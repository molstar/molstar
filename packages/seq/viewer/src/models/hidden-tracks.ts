/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import type { Node } from '@molstar/core/repr/tree';
import type { HiddenTrack } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';

/** A track the spec hides, with the section it is in (`SeqViewer.hiddenTracks`). */
export interface HiddenTrackEntry extends HiddenTrack {
  /** The section's address. */
  readonly section: string;
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
        for (const t of c.children ?? []) {
          if (t.kind !== 'track' || t.id === undefined || t.params?.hidden !== true) continue;
          const label = t.params?.label;
          (out ??= []).push({
            address: `${section}/${t.id}`,
            label: typeof label === 'string' ? label : t.id,
            section,
          });
        }
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
