/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { StructureSource, type StructureFormat } from '@molstar/seq-molstar-poc';
import { merge } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import type { SequenceView3D } from '../sequence/link.js';
import type { HostActions, HostStructure, StructureColoring } from '../sequence/model.js';
import { structureNames } from './naming.js';
import type { SourceScope } from './state.js';

// The sequence's side of the viewer (sequence.md §2–§5, plan 2.4): which structure it follows, the parse it shares
// with the 3D view, and the link's 3D side. Reads come from the plugin here; every write (highlight, selection, focus,
// the camera, coloring from a track's menu) is an operation with origin `link` (sequence.md §4, review F21), so it is
// recorded and, where it changes content, a history step.

export interface SequenceSource {
  /** The structure the sequence follows (sequence.md §3: the focused one, else the latest selected, else the first). */
  readonly structure: ReadonlyState<HostStructure | undefined>;
  /** Follows this structure (by ref) whatever the rule says; undefined: the rule. */
  setFollowed(ref: string | undefined): void;
  readonly followed: ReadonlyState<string | undefined>;
  /** What the track menus may do to the 3D view: through the API. */
  readonly actions: HostActions;
  /** The link's 3D side over a structure. */
  view3d(structure: Structure, resetCamera?: () => void): SequenceView3D;
}

/** The format the plugin parsed a trajectory from, as the sequence's structure source names it. */
function formatOf(plugin: PluginContext, s: StructureRef): StructureFormat | undefined {
  const trajectory = s.model?.trajectory?.cell;
  if (!trajectory) return undefined;
  const name = trajectory.transform.transformer.definition.name;
  if (name === 'trajectory-from-pdb') return 'pdb';
  if (name !== 'trajectory-from-mmcif') return undefined;
  // mmCIF from text or from BinaryCIF: the data cell under the parse says which.
  const parent = plugin.state.data.cells.get(trajectory.transform.parent);
  const data = parent && plugin.state.data.cells.get(parent.transform.parent);
  return data?.obj?.data instanceof Uint8Array ? 'bcif' : 'mmcif';
}

/** Runs writes one after another; a write still waiting when a newer one comes is dropped. */
function channel(): (write: () => Promise<unknown>) => void {
  let tail: Promise<unknown> = Promise.resolve();
  let latest = 0;
  return (write) => {
    const id = ++latest;
    tail = tail.then(() => (id === latest ? write() : undefined)).catch(() => undefined);
  };
}

export function createSequenceSource(scope: SourceScope, plugin: PluginContext, link: ViewerApi): SequenceSource {
  const structure = observableState<HostStructure | undefined>(undefined, { compare: (a, b) => a?.key === b?.key });
  const followed = observableState<string | undefined>(undefined);
  let latestSelected: string | undefined;
  const keys = new WeakMap<Structure, number>();
  let serial = 0;
  const keyOf = (ref: string, s: Structure) => {
    let k = keys.get(s);
    if (k === undefined) keys.set(s, (k = ++serial));
    return `${ref}:${k}`;
  };

  const pick = (): StructureRef | undefined => {
    const list = plugin.managers.structure.hierarchy.current.structures.filter((s) => !!s.cell.obj?.data);
    if (!list.length) return undefined;
    const byRef = (ref: string | undefined) => (ref ? list.find((s) => s.cell.transform.ref === ref) : undefined);
    const focus = plugin.managers.structure.focus.current?.loci.structure.root;
    return (
      byRef(followed.value) ??
      (focus ? list.find((s) => s.cell.obj!.data.root === focus) : undefined) ??
      byRef(latestSelected) ??
      list[0]
    );
  };

  const update = () => {
    const s = pick();
    const data = s?.cell.obj?.data;
    const trajectory = s?.model?.trajectory?.cell.obj?.data;
    const format = s && formatOf(plugin, s);
    if (!s || !data || !trajectory || !format) {
      structure.set(undefined);
      return;
    }
    const ref = s.cell.transform.ref;
    const key = keyOf(ref, data);
    if (structure.value?.key === key) return;
    structure.set({
      key,
      title: structureNames(plugin.managers.structure.hierarchy.current.structures).get(ref) ?? ref,
      source: StructureSource.of(trajectory, format),
      structure: data,
    });
  };

  scope.add(
    plugin.managers.structure.selection.events.changed.subscribe(() => {
      // The structure whose selection changed last (an addition, not a clear).
      const list = plugin.managers.structure.hierarchy.current.structures;
      for (const s of list) {
        const data = s.cell.obj?.data;
        if (!data) continue;
        const loci = plugin.managers.structure.selection.getLoci(data);
        if (StructureElement.Loci.is(loci) && loci.elements.length) latestSelected = s.cell.transform.ref;
      }
      update();
    }),
  );
  scope.add(
    merge(
      plugin.managers.structure.hierarchy.behaviors.selection,
      plugin.managers.structure.focus.behaviors.current,
    ).subscribe(() => update()),
  );
  scope.add(followed.asReadonly.subscribe(() => update()));

  const actions: HostActions = {
    coloring() {
      // The polymer's representations (not the ligands' ball and stick): one color theme among them, or none of the three.
      let name: string | undefined;
      for (const st of plugin.managers.structure.hierarchy.current.structures) {
        for (const c of st.components) {
          for (const r of c.representations) {
            const params = r.cell.transform.params as { type?: { name: string }; colorTheme?: { name: string } };
            if (params?.type?.name === 'ball-and-stick') continue;
            const theme = params?.colorTheme?.name;
            if (name !== undefined && theme !== name) return undefined;
            name = theme;
          }
        }
      }
      // The default coloring under any of its names; any other theme by its own.
      return name === 'default' || name === 'chain-id' || name === 'uniform' ? 'default' : name;
    },
    async colorStructure(coloring: StructureColoring) {
      // sequence.md §5: a track's "Color by this track" is `representation.color` on the polymer, as the link.
      const theme = coloring === 'default' ? 'chain-id' : coloring;
      await link.representation.color({
        target: 'polymer',
        theme,
        structure: structure.value?.key.split(':')[0] ?? '',
      });
    },
  };

  return {
    structure: structure.asReadonly,
    followed: followed.asReadonly,
    setFollowed: (ref) => followed.set(ref),
    actions,
    view3d(target: Structure, resetCamera?: () => void): SequenceView3D {
      const onto = (l: unknown): StructureElement.Loci | null =>
        StructureElement.Loci.is(l) && l.structure.root === target.root ? StructureElement.Loci.remap(l, target) : null;
      // The link calls without waiting; its writes apply in the order it sent them, one channel at a time, and a write
      // superseded while it waited is skipped (only the latest hover or selection matters).
      const highlight = channel();
      const select = channel();
      const focus = channel();
      return {
        highlight: (loci) =>
          highlight(() => (loci?.elements.length ? link.highlight.set({ target: { loci } }) : link.highlight.clear())),
        select: (loci) =>
          select(() =>
            loci?.elements.length ? link.select.set({ target: { loci }, mode: 'set' }) : link.select.clear(),
          ),
        focus: (loci) => focus(() => link.camera.frame({ target: { loci } })),
        onHover: (fn) => {
          const sub = plugin.behaviors.interaction.hover.subscribe((e) => fn(onto(e.current.loci)));
          return () => sub.unsubscribe();
        },
        onSelect: (fn) => {
          const sub = plugin.managers.structure.selection.events.changed.subscribe(() => {
            const loci = plugin.managers.structure.selection.getLoci(target);
            fn(
              StructureElement.Loci.is(loci) && loci.elements.length ? StructureElement.Loci.remap(loci, target) : null,
            );
          });
          return () => sub.unsubscribe();
        },
        ...(resetCamera ? { resetCamera } : {}),
      };
    },
  };
}
