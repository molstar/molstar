/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-viewer/link-structure (_spec/README.md "Structure link"): the wiring between a viewer and a structure
// viewer a host translates for (its `side`, plain data): the viewer's hover, selection and jumps out as targets on
// units, the side's hover, selection and focus in, through `@molstar/seqx-core/structure`. What it adds over the
// viewer's events is per link: echoes never sent back, one selection message a frame, a stale `null` ignored, the
// clearing at either end. No behavior, no port: a host that never calls it carries none of it.

import { EmptyKeys, isEmptyLocation, type Location, locationsEqual } from '@molstar/seqx-core/query';
import {
  DefaultSyncLimits,
  fromUnits,
  type StructureBinding,
  type StructureModel,
  type StructureTarget,
  type SyncLimits,
  toUnits,
} from '@molstar/seqx-core/structure';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import type { SeqViewer } from './viewer.js';

export type { StructureBinding, StructureTarget, SyncLimits };

/** A structure viewer as its host translates it: what Seq\* tells it, and its own events. */
export interface StructureSide {
  /** Seq\*'s hover, selection and jumps, as targets on units (`null`: cleared). */
  highlight(target: StructureTarget | null): void;
  select(target: StructureTarget | null): void;
  focus(target: StructureTarget, options: { readonly durationMs: number }): void;
  /** The structure viewer's own events, translated by the host; each returns its unsubscribe. */
  onHover?(fn: (target: StructureTarget | null) => void): () => void;
  onSelect?(fn: (target: StructureTarget | null) => void): () => void;
  onFocus?(fn: (target: StructureTarget) => void): () => void;
  /** Aborted when the structure viewer goes: the link clears and disposes itself. */
  readonly signal?: AbortSignal;
}

/** Which way a path is linked: both ways, Seq\* to the side only, the side to Seq\* only, or not. */
export type LinkDirection = 'both' | 'out' | 'in' | 'none';

export interface LinkStructureOptions {
  /** The complexes the link is bound to and the maps their paths take; absent: every complex, no `via`. */
  readonly complexes?: readonly StructureBinding[];
  readonly hover?: LinkDirection;
  readonly select?: LinkDirection;
  readonly focus?: LinkDirection;
  /** What a selection from the side does to the view: revealed (the zoom kept), focused, or nothing. */
  readonly revealOnSelect?: 'reveal' | 'focus' | 'none';
  /** Whether back and forward bring the camera (a step that names a Location). */
  readonly followJumps?: boolean;
  /** The size gate a direction (`DefaultSyncLimits`). */
  readonly limits?: { readonly hover?: SyncLimits; readonly select?: SyncLimits };
  /** A molecule unit lit whole (the default) or on its atoms. */
  readonly ligand?: 'whole' | 'atoms';
}

const out = (d: LinkDirection | undefined) => d === undefined || d === 'both' || d === 'out';
const inward = (d: LinkDirection | undefined) => d === undefined || d === 'both' || d === 'in';
const AtomsSuffix = '/(atoms)';

/**
 * Links a viewer to a structure viewer's side both ways (_spec/README.md "Structure link"); returns the disposer. The
 * disposer, the side's `signal` or the viewer retiring clears both ends once (`highlight(null)`, `select(null)` for the
 * directions sent, the link's hover dropped) and ends the link.
 */
export function linkStructure(viewer: SeqViewer, side: StructureSide, options: LinkStructureOptions = {}): () => void {
  const complexes = options.complexes;
  const hoverLimits = options.limits?.hover ?? DefaultSyncLimits.hover;
  const selectLimits = options.limits?.select ?? DefaultSyncLimits.select;
  const ligand = options.ligand ?? 'whole';
  const reveal = options.revealOnSelect ?? 'reveal';
  const followJumps = options.followJumps ?? true;
  /** The owner of the hover the side gives (EC-141): the pointer taking the hover takes it from the link. */
  const owner = { link: 'structure' };
  const offs: (() => void)[] = [];
  let disposed = false;
  /** The last hover sent (`null`: cleared or never). */
  let sentHover: StructureTarget | null = null;
  /** The last selection sent, and the last set from the side (echoes by value). */
  let sentSelection: StructureTarget | null = null;
  let fromSide: Location | undefined;
  /** The frame callback a selection message waits for (one a frame). */
  let selectionDue: number | undefined;
  /** While a message from the side is applied: what it moves is not sent back. */
  let inbound = false;

  const model = (): StructureModel => {
    const m = viewer.spec.model.value;
    const overlay = viewer.store.currentOverlay(m);
    // no source's nodes held (none yet, before the first frame: an empty overlay that knows no path): the model alone
    return overlay.nodes.length ? { model: m, overlay } : m;
  };
  const onUnits = (loc: Location, limits: SyncLimits): StructureTarget | null => {
    if (isEmptyLocation(loc)) return null;
    const t = toUnits(model(), loc, complexes, { limits, ligand }).target;
    return isEmptyLocation(t) ? null : t;
  };
  const same = (a: StructureTarget | null, b: StructureTarget | null) =>
    a === b || (a !== null && b !== null && locationsEqual(a, b));
  /** Atoms kept as atoms where a section or a block shows a unit's atoms form; else collapsed onto residues. */
  const showsAtoms = () => {
    for (const s of viewer.view.resolved.value?.sections ?? []) {
      if (s.axis.system.endsWith(AtomsSuffix)) return true;
      for (const b of s.blocks ?? []) if (b.axis.system.endsWith(AtomsSuffix)) return true;
    }
    return false;
  };
  const onSequences = (target: StructureTarget) =>
    fromUnits(model(), target, complexes, { atoms: showsAtoms() }).location;
  const apply = (fn: () => void) => {
    if (disposed) return;
    inbound = true;
    try {
      fn();
    } finally {
      inbound = false;
    }
  };

  // ---- out: the viewer's hover, selection and jumps ----

  if (out(options.hover))
    offs.push(
      viewer.hover.version.listen(() => {
        // the side's own hover, carried in, is never sent back: it only ends the pointer's sent before it
        const loc = viewer.hover.owner === owner ? undefined : viewer.hover.location();
        const t = loc ? onUnits(loc, hoverLimits) : null;
        if (same(t, sentHover)) return;
        sentHover = t;
        side.highlight(t);
      }),
    );

  const sendSelection = () => {
    selectionDue = undefined;
    if (disposed) return;
    const shown = viewer.selection.value.value.shown;
    if (fromSide && locationsEqual(shown, fromSide)) return;
    fromSide = undefined;
    const t = onUnits(shown, selectLimits);
    if (same(t, sentSelection)) return;
    sentSelection = t;
    side.select(t);
  };
  if (out(options.select))
    offs.push(
      viewer.selection.value.listen(() => {
        // a drag moves it every pointer event: one message a frame
        if (selectionDue !== undefined) return;
        selectionDue = requestFrameCallback(sendSelection);
      }),
    );

  if (out(options.focus))
    offs.push(
      viewer.view.jumps.listen(() => {
        const j = viewer.view.lastJump;
        if (inbound || !j?.location) return;
        if ((j.kind === 'back' || j.kind === 'forward') && !followJumps) return;
        const t = onUnits(j.location, selectLimits);
        if (t) side.focus(t, { durationMs: j.durationMs });
      }),
    );

  // ---- in: the side's hover, selection and focus ----

  if (inward(options.hover) && side.onHover)
    offs.push(
      side.onHover((target) =>
        apply(() => {
          const loc = target ? onSequences(target) : undefined;
          // a stale null (the pointer took the hover since) clears nothing: only the link's own
          if (loc && !isEmptyLocation(loc)) viewer.hover.setLocation(loc, owner);
          else viewer.hover.clear(owner);
        }),
      ),
    );

  if (inward(options.select) && side.onSelect)
    offs.push(
      side.onSelect((target) =>
        apply(() => {
          const loc = target ? onSequences(target) : undefined;
          if (!loc || isEmptyLocation(loc)) {
            viewer.selection.clear();
            fromSide = viewer.selection.value.value.shown;
            sentSelection = null;
            return;
          }
          viewer.selection.set({ elements: EmptyKeys, location: loc, shown: loc });
          fromSide = loc;
          sentSelection = target;
          if (reveal === 'reveal') viewer.addresses.revealLocation(loc);
          else if (reveal === 'focus') viewer.addresses.focusLocation(loc);
        }),
      ),
    );

  if (inward(options.focus) && side.onFocus)
    offs.push(
      side.onFocus((target) =>
        apply(() => {
          const loc = onSequences(target);
          if (!isEmptyLocation(loc)) viewer.addresses.focusLocation(loc);
        }),
      ),
    );

  // ---- the end: the disposer, the side going, the viewer retiring ----

  const dispose = () => {
    if (disposed) return;
    disposed = true;
    for (const off of offs.splice(0)) off();
    if (selectionDue !== undefined) cancelFrameCallback(selectionDue);
    if (out(options.hover)) side.highlight(null);
    if (out(options.select)) side.select(null);
    viewer.hover.clear(owner);
  };
  const ends = [side.signal, viewer.retiredSignal];
  for (const s of ends) {
    if (!s) continue;
    if (s.aborted) {
      dispose();
      break;
    }
    s.addEventListener('abort', dispose, { once: true });
    offs.push(() => s.removeEventListener('abort', dispose));
  }
  return dispose;
}
