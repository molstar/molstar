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

import {
  EmptyKeys,
  isEmptyLocation,
  type Location,
  location,
  locationsEqual,
  unionLocations,
} from '@molstar/seqx-core/query';
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
  /**
   * `at` is when the structure viewer's hover happened (`performance.now()`); absent, it is taken as now. A hover older
   * than the pointer's last own hover is dropped (a late message never takes the hover from the pointer).
   */
  onHover?(fn: (target: StructureTarget | null, at?: number) => void): () => void;
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
    return { model: m, overlay: viewer.store.currentOverlay(m) };
  };
  /** Whether the last hover and selection sent wait for a described map (EC-169): sent again once it resolves. */
  let waitingHover = false;
  /** When the hover last changed other than by the side's (the pointer's own, or its leaving). */
  let lastOwnHover = Number.NEGATIVE_INFINITY;
  let waitingSelection = false;
  /** A Location on units, and whether a described map it could take waits for its description. */
  const unitsOf = (loc: Location, limits: SyncLimits): { t: StructureTarget | null; pending: boolean } => {
    if (isEmptyLocation(loc)) return { t: null, pending: false };
    const r = toUnits(model(), loc, complexes, { limits, ligand });
    return { t: isEmptyLocation(r.target) ? null : r.target, pending: !!r.pending };
  };
  const onUnits = (loc: Location, limits: SyncLimits): StructureTarget | null => unitsOf(loc, limits).t;
  const same = (a: StructureTarget | null, b: StructureTarget | null) =>
    a === b || (a !== null && b !== null && locationsEqual(a, b));
  /** The systems the view's sections and blocks show: a unit's atoms form among them where it is shown. */
  const shown = () => {
    const out = new Set<string>();
    for (const s of viewer.view.resolved.value?.sections ?? []) {
      out.add(s.axis.system);
      for (const b of s.blocks ?? []) out.add(b.axis.system);
    }
    return out;
  };
  /**
   * A target on Seq*'s systems: a unit's atoms kept as atoms where its atoms form is shown, else collapsed onto its
   * residues, unit by unit (a part on a shown system read with `atoms`, the rest without).
   */
  const onSequences = (target: StructureTarget): Location => {
    const systems = shown();
    const kept = target.parts.filter((p) => systems.has(p.system));
    const m = model();
    if (!kept.length) return fromUnits(m, target, complexes).location;
    const rest = target.parts.filter((p) => !systems.has(p.system));
    const atoms = fromUnits(m, location({ parts: kept, rects: target.rects }), complexes, { atoms: true }).location;
    return rest.length ? unionLocations(atoms, fromUnits(m, location({ parts: rest }), complexes).location) : atoms;
  };
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

  const sendHover = () => {
    // the side's own hover, carried in, is never sent back: it only ends the pointer's sent before it
    const loc = viewer.hover.owner === owner ? undefined : viewer.hover.location();
    const { t, pending } = loc ? unitsOf(loc, hoverLimits) : { t: null, pending: false };
    // waiting for a described map: held (the last sent stays), never sent as nothing; sent once it resolves
    waitingHover = pending;
    if (pending && t === null) return;
    if (same(t, sentHover)) return;
    sentHover = t;
    side.highlight(t);
  };
  if (out(options.hover))
    offs.push(
      viewer.hover.version.listen(() => {
        if (viewer.hover.owner !== owner) lastOwnHover = performance.now();
        sendHover();
      }),
    );

  const sendSelection = () => {
    selectionDue = undefined;
    if (disposed) return;
    const shown = viewer.selection.value.value.shown;
    if (fromSide && locationsEqual(shown, fromSide)) return;
    fromSide = undefined;
    const { t, pending } = unitsOf(shown, selectLimits);
    waitingSelection = pending;
    if (pending && t === null) return;
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

  // a description arriving (the store's nodes): what waited for it is sent again
  if (out(options.hover) || out(options.select))
    offs.push(
      viewer.store.store.nodeVersion.listen(() => {
        if (disposed) return;
        if (waitingHover) sendHover();
        if (waitingSelection && selectionDue === undefined) sendSelection();
      }),
    );

  // ---- in: the side's hover, selection and focus ----

  if (inward(options.hover) && side.onHover)
    offs.push(
      side.onHover((target, at = performance.now()) =>
        apply(() => {
          // a hover the pointer has since overtaken (it happened before the pointer's last own hover) is stale
          if (at < lastOwnHover) return;
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
