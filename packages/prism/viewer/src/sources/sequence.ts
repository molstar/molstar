/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Representation } from '@molstar/graphics/repr/representation';
import { Model, Structure, StructureElement, type Trajectory } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { LociLabelProvider } from '@molstar/plugin/state/manager/loci-label';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { StructureFormat, StructureSource } from '../sequence/molstar/source.js';
import { merge } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import type { SequenceView3D } from '../sequence/link.js';
import type { DataByUrl, HostActions, HostStructure, StructureColoring } from '../sequence/model.js';
import { secondaryAssignmentOf } from '../sequence/tracks.js';
import { cifFrameOf, createDataByUrl } from './data-by-url.js';
import { structureNames } from './naming.js';
import { onceATask, type SourceScope } from './state.js';

// The sequence's side of the viewer (sequence.md §2–§5): which structure it follows, the parse it shares
// with the 3D view, and the link's 3D side. Reads come from the plugin here; every write (highlight, selection, focus,
// the camera, coloring from a track's menu) is an operation with origin `link` (sequence.md §4), so it is
// recorded and, where it changes content, a history step.

export interface SequenceSource {
  /** The structure the sequence follows (sequence.md §3: the focused one, else the latest selected, else the first). */
  readonly structure: ReadonlyState<HostStructure | undefined>;
  /** Follows this structure (by ref) whatever the rule says; undefined: the rule. */
  setFollowed(ref: string | undefined): void;
  readonly followed: ReadonlyState<string | undefined>;
  /**
   * Shows this structure (by ref) until the person focuses or loads another, follows another (`setFollowed`), or it goes:
   * a loaded SeqViewSpec's structure (sequence.md §6). Undefined: the rule again.
   */
  show(ref: string | undefined): void;
  /** The structure `show` shows; undefined once that ended. */
  readonly shows: ReadonlyState<string | undefined>;
  /** What the track menus may do to the 3D view: through the API. */
  readonly actions: HostActions;
  /** The plugin's data cells and their parses by address: the sequence's parse cache (sequence.md §2). */
  readonly dataByUrl: DataByUrl;
  /** The link's 3D side over a structure. */
  view3d(structure: Structure, resetCamera?: () => void): SequenceView3D;
  /**
   * Whether the sequence's hover maps into the structure now (the link lit residues in 3D): then the hover label is
   * Mol*'s, of those residues (its label providers label every highlight), and Seq*'s own label is for what does not
   * map (sequence.md §10).
   */
  readonly hoverMapped: ReadonlyState<boolean>;
  /** The band in the page's shell (its dock shown or not), once the shell has mounted it; headless: none. */
  readonly band: ReadonlyState<SequenceBand | undefined>;
  /** The shell's band (`sequence/shell-ui.tsx`); returns its removal. */
  setBand(band: SequenceBand): () => void;
}

/** Whether the page shows the band, and showing or hiding it (the Sequence menu's No sequence). */
export interface SequenceBand {
  readonly shown: ReadonlyState<boolean>;
  show(on: boolean): void;
  /** Opens the command menu at Add track (the commands that add a track to the structure shown). */
  openAddTrack?(): void;
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

/** Writes, one after another; a write still waiting when a newer one comes is dropped. */
interface Channel {
  (write: () => Promise<unknown>): void;
  /** A write is queued or running: what it writes is not in the plugin yet. */
  readonly pending: boolean;
}

function channel(): Channel {
  let tail: Promise<unknown> = Promise.resolve();
  let latest = 0;
  let open = 0;
  const send = (write: () => Promise<unknown>) => {
    const id = ++latest;
    open++;
    tail = tail
      .then(() => (id === latest ? write() : undefined))
      .catch(() => undefined)
      .finally(() => {
        open--;
      });
  };
  return Object.defineProperty(send, 'pending', { get: () => open > 0 }) as Channel;
}

export function createSequenceSource(scope: SourceScope, plugin: PluginContext, link: ViewerApi): SequenceSource {
  // A new key or a new frame of it (the same key, another Structure: sequence.md §2).
  const structure = observableState<HostStructure | undefined>(undefined, {
    compare: (a, b) => a?.key === b?.key && a?.structure === b?.structure && a?.assignment === b?.assignment,
  });
  const followed = observableState<string | undefined>(undefined);
  const shows = observableState<string | undefined>(undefined);
  const hoverMapped = observableState(false);
  /** What was there when `show` was called: a structure added since, or another focus, ends it. */
  let showsSince: { refs: ReadonlySet<string>; focus: unknown } | undefined;
  let latestSelected: string | undefined;
  /** What each structure had selected at the last change, to see which selection grew. */
  const selectedSizes = new Map<string, number>();
  const keys = new WeakMap<Structure, number>();
  let serial = 0;
  const keyOf = (ref: string, s: Structure) => {
    let k = keys.get(s);
    if (k === undefined) keys.set(s, (k = ++serial));
    return `${ref}:${k}`;
  };

  /** Two frames of one structure: the same units, atoms and hierarchy (only coordinates may differ). */
  const sameAtoms = (a: Structure, b: Structure) =>
    Structure.areUnitIdsAndIndicesEqual(a, b) &&
    a.models.length === b.models.length &&
    a.models.every((m, i) => Model.areHierarchiesEqual(m, b.models[i]));

  // The followed structure's parse, one per trajectory (its models and structures cached across keys), with the CIF
  // block it came from: SIFTS numberings, a predicted model's origin and the scheme's auth labels read from it.
  const sources = new WeakMap<Trajectory, StructureSource>();
  const sourceOf = (s: StructureRef, trajectory: Trajectory, format: StructureFormat): StructureSource => {
    const hit = sources.get(trajectory);
    if (hit?.format === format) return hit;
    const source: StructureSource = {
      kind: 'molstar:structure-source',
      format,
      trajectory,
      frame: cifFrameOf(plugin, s.model!.trajectory!.cell),
      structures: new Map(),
    };
    sources.set(trajectory, source);
    return source;
  };

  const pick = (): StructureRef | undefined => {
    const list = plugin.managers.structure.hierarchy.current.structures.filter((s) => !!s.cell.obj?.data);
    if (!list.length) return undefined;
    const byRef = (ref: string | undefined) => (ref ? list.find((s) => s.cell.transform.ref === ref) : undefined);
    const focus = plugin.managers.structure.focus.current?.loci.structure.root;
    return (
      byRef(shows.value) ??
      byRef(followed.value) ??
      (focus ? list.find((s) => s.cell.obj!.data.root === focus) : undefined) ??
      byRef(latestSelected) ??
      list[0]
    );
  };

  /** Ends `show` when the person moved on: another structure loaded or focused, or its own gone (sequence.md §6). */
  const endShow = () => {
    const ref = shows.value;
    if (!ref || !showsSince) return;
    const list = plugin.managers.structure.hierarchy.current.structures;
    const own = list.find((s) => s.cell.transform.ref === ref);
    const added = list.some((s) => !showsSince!.refs.has(s.cell.transform.ref));
    const focus = plugin.managers.structure.focus.current;
    const focusedElsewhere =
      !!focus && focus !== showsSince.focus && focus.loci.structure.root !== own?.cell.obj?.data.root;
    if (!own || added || focusedElsewhere) {
      showsSince = undefined;
      shows.set(undefined);
    }
  };

  const update = () => {
    endShow();
    const s = pick();
    const data = s?.cell.obj?.data;
    const trajectory = s?.model?.trajectory?.cell.obj?.data;
    const format = s && formatOf(plugin, s);
    if (!s || !data || !trajectory || !format) {
      structure.set(undefined);
      return;
    }
    const ref = s.cell.transform.ref;
    const now = structure.value;
    // A frame of the structure shown (a trajectory step, the model-index animation): the same atoms, so the same key:
    // the sequence keeps its document and the link follows the new Structure (sequence.md §2).
    if (now?.ref === ref && now.structure !== data && now.source.trajectory === trajectory) {
      if (sameAtoms(now.structure, data)) {
        now.source.drawn = data;
        keys.set(data, keys.get(now.structure)!);
        structure.set({ ...now, structure: data, assignment: secondaryAssignmentOf(data) });
        return;
      }
    }
    // The structure the 3D view draws, for what the sequence reads of its computed properties (its secondary structure).
    sourceOf(s, trajectory, format).drawn = data;
    const key = keyOf(ref, data);
    if (now?.key === key) return;
    structure.set({
      assignment: secondaryAssignmentOf(data),
      key,
      ref,
      title: structureNames(plugin.managers.structure.hierarchy.current.structures).get(ref) ?? ref,
      // Seq*'s source over the 3D view's parse (`StructureSource.of`, made here: Seq* is the sequence's own chunk).
      source: sourceOf(s, trajectory, format),
      structure: data,
    });
  };

  scope.add(
    // Once a task: a select over many structures tells `changed` once a structure, and this reads every structure.
    onceATask(plugin.managers.structure.selection.events.changed).subscribe(() => {
      // The structure whose selection grew in this change (an addition, not a clear); the last of them if several.
      const list = plugin.managers.structure.hierarchy.current.structures;
      const seen = new Set<string>();
      for (const s of list) {
        const data = s.cell.obj?.data;
        if (!data) continue;
        const ref = s.cell.transform.ref;
        const loci = plugin.managers.structure.selection.getLoci(data);
        const size = StructureElement.Loci.is(loci) ? StructureElement.Loci.size(loci) : 0;
        if (size > (selectedSizes.get(ref) ?? 0)) latestSelected = ref;
        else if (size === 0 && latestSelected === ref) latestSelected = undefined;
        selectedSizes.set(ref, size);
        seen.add(ref);
      }
      for (const ref of [...selectedSizes.keys()]) if (!seen.has(ref)) selectedSizes.delete(ref);
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
  // The shown structure's secondary structure set otherwise in the 3D view (its custom properties updated, the Structure
  // the same): told again, so the sequence's track follows (`secondaryAssignmentOf` names the settings in its data).
  scope.add(
    plugin.state.data.events.object.updated.subscribe((e) => {
      const now = structure.value;
      const data = e.obj?.data;
      if (!now || !(data instanceof Structure) || data.root !== now.structure.root) return;
      structure.set({ ...now, assignment: secondaryAssignmentOf(now.structure) });
    }),
  );
  scope.add(shows.asReadonly.subscribe(() => update()));

  const actions: HostActions = {
    coloring() {
      // The followed structure's polymer representations (not the ligands' ball and stick): one color theme among them,
      // or none of the three.
      const ref = structure.value?.ref;
      if (!ref) return undefined;
      let name: string | undefined;
      for (const st of plugin.managers.structure.hierarchy.current.structures) {
        if (st.cell.transform.ref !== ref) continue;
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
      // Of the structure the sequence shows; with none, nothing (an empty structure would mean every one).
      const ref = structure.value?.ref;
      if (!ref) return;
      const theme = coloring === 'default' ? 'chain-id' : coloring;
      await link.representation.color({ target: 'polymer', theme, structure: ref });
    },
  };

  const band = observableState<SequenceBand | undefined>(undefined);
  return {
    band: band.asReadonly,
    setBand: (b) => {
      band.set(b);
      return () => {
        if (band.value === b) band.set(undefined);
      };
    },
    hoverMapped: hoverMapped.asReadonly,
    structure: structure.asReadonly,
    dataByUrl: createDataByUrl(scope, plugin),
    followed: followed.asReadonly,
    setFollowed: (ref) => {
      // Following a structure ends what `show` showed.
      showsSince = undefined;
      shows.set(undefined);
      followed.set(ref);
    },
    shows: shows.asReadonly,
    show: (ref) => {
      showsSince = ref
        ? {
            refs: new Set(plugin.managers.structure.hierarchy.current.structures.map((s) => s.cell.transform.ref)),
            focus: plugin.managers.structure.focus.current,
          }
        : undefined;
      shows.set(ref);
    },
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
        highlight: (loci) => {
          hoverMapped.set(!!loci?.elements.length);
          return highlight(() =>
            loci?.elements.length ? link.highlight.set({ target: { loci } }) : link.highlight.clear(),
          );
        },
        select: (loci) =>
          select(() =>
            loci?.elements.length ? link.select.set({ target: { loci }, mode: 'set' }) : link.select.clear(),
          ),
        focus: (loci) => focus(() => link.camera.frame({ target: { loci } })),
        onHover: (fn) => {
          // Mol* emits the hover after every pick while the pointer moves: the same element again is not passed on.
          let last: Representation.Loci | undefined;
          const sub = plugin.behaviors.interaction.hover.subscribe((e) => {
            if (last && Representation.Loci.areEqual(last, e.current)) return;
            last = e.current;
            fn(onto(e.current.loci));
          });
          return () => sub.unsubscribe();
        },
        onSelect: (fn) => {
          // Once a task, as the source's own listener: a select over S structures says `changed` once per structure,
          // and this remaps the whole target selection each time.
          const sub = onceATask(plugin.managers.structure.selection.events.changed).subscribe(() => {
            const loci = plugin.managers.structure.selection.getLoci(target);
            fn(
              StructureElement.Loci.is(loci) && loci.elements.length ? StructureElement.Loci.remap(loci, target) : null,
            );
          });
          return () => sub.unsubscribe();
        },
        annotateHover: (read) => {
          // After every other provider: the residue first, then the track's lines (sequence.md §10). Only while the
          // sequence's hover maps: what the link lit is then the 3D highlight of this structure.
          const provider: LociLabelProvider = {
            label: (loci) =>
              hoverMapped.value && StructureElement.Loci.is(loci) && loci.structure.root === target.root
                ? read()
                : undefined,
            priority: -1000,
          };
          const labels = plugin.managers.lociLabels;
          labels.addProvider(provider);
          return {
            // A highlight on its way labels itself with the latest lines when it lands: labelling the one before it
            // with them would pair the last residue with the new track lines, and label twice.
            refresh: () => {
              if (hoverMapped.value && !highlight.pending) labels.refresh();
            },
            dispose: () => labels.removeProvider(provider),
          };
        },
        ...(resetCamera ? { resetCamera } : {}),
      };
    },
  };
}
