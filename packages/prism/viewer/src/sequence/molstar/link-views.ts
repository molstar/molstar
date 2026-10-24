/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The link between a viewer and a 3D view over `ctx.systems` (seqspec:libraries#the-mol-adapter): hover and selection
// both ways, sent by size, focus. Prism's sequence link (../link.ts) and the Seq* Storybook's Mol* demo run on it.

import { Extents, type Loci } from '@molstar/seq-core/model';
import { type Coordinates, type LociSize, lociSize, type SeqViewerContext } from '@molstar/seq-viewer';
import { OrderedSet } from '@molstar/core/data/int';
import { StructureElement, type Unit } from '@molstar/model/model/structure';
import type { ChainInstance } from './chains.js';
import type { StructureLink } from './link.js';

const sameLoci = (a: StructureElement.Loci | null, b: StructureElement.Loci | null) =>
  a === null || b === null ? a === b : StructureElement.Loci.areEqual(a, b);

/** Seq* loci (a selection's parts, a hover) as one set of runs. Edge rate. */
export function lociRuns(list: readonly (Loci | null)[]): Extents {
  const e = Extents.create(8);
  for (const l of list) if (l) Extents.addLoci(e, l);
  return e;
}

/** A 3D view as the link needs it: Prism's 3D view implements it, as does the Seq* Storybook's bare Mol* plugin. */
export interface StructureView3D {
  highlight(loci: StructureElement.Loci | null): void;
  select(loci: StructureElement.Loci | null): void;
  focus(loci: StructureElement.Loci): void;
  /** Mol* hover and selection, as loci of the linked structure (null: off it). Returns the unsubscribe. */
  onHover(fn: (loci: StructureElement.Loci | null) => void): () => void;
  onSelect(fn: (loci: StructureElement.Loci | null) => void): () => void;
}

/** Runs in a link's chain systems as `Coordinates` by key, what the viewer's channels take. Edge rate. */
export function coordinatesOf(link: StructureLink, runs: Extents): Coordinates {
  const systems: string[] = [];
  // Each system's index by key (an `indexOf` a run was runs × systems).
  const index = new Map<string, number>();
  const out: number[] = [];
  for (let i = 0; i < runs.count; i++) {
    const m = link.memberOf(runs, i);
    if (!m) continue;
    let k = index.get(m.key);
    if (k === undefined) {
      k = systems.push(m.key) - 1;
      index.set(m.key, k);
    }
    out.push(k, runs.start[i], runs.stop[i]);
  }
  return { systems, runs: Int32Array.from(out) };
}

/**
 * How the link reacts to a Seq* set by its size (seqspec:libraries#the-mol-adapter), read before it resolves: within
 * `maxItems` items and `maxLength` positions, residue by residue; past either, `over`: the whole chains it touches, or
 * nothing.
 */
export interface LinkReaction {
  readonly maxItems: number;
  readonly maxLength: number;
  readonly over: 'chains' | 'skip';
  /**
   * A non-polymer chain's position (one heavy atom) lights and selects its whole residue in 3D, as the unit holds it,
   * hydrogens and alternate locations included: a ligand as a whole, a glycan's monosaccharide (a branched chain's
   * residue, not the whole glycan); waters are not linked. False: the atoms the positions name. A polymer's position
   * is its residue either way. Default true (seqspec:libraries#the-mol-adapter).
   */
  readonly wholeResidue: boolean;
}
/** The link's reactions (`linkViews`' options): its composition, as the viewer spec is the viewer's. */
export interface LinkViewsOptions {
  /** Seq*'s hover to 3D and back; false: neither. */
  readonly hover?: false | Partial<LinkReaction>;
  /** Selections both ways; false: neither. */
  readonly select?: false | Partial<LinkReaction>;
  /** A selection focuses the other view. */
  readonly focus?: boolean;
}
export const DefaultLinkViewsOptions: {
  readonly hover: LinkReaction;
  readonly select: LinkReaction;
  readonly focus: boolean;
} = {
  hover: { maxItems: 2_000, maxLength: 20_000, over: 'chains', wholeResidue: true },
  select: { maxItems: 20_000, maxLength: 200_000, over: 'chains', wholeResidue: true },
  focus: true,
};

/**
 * The whole chains a set's runs touch, as one Mol* loci (a large set's summary). Edge rate. Each touched chain is its
 * own interval: Mol* merges consecutive one-atom chains (ions) into one unit, and the chains between two touched ones
 * are not touched. The units are the structure the link maps now (a frame of the one it was attached to).
 */
export function wholeChains(link: StructureLink, runs: Extents): StructureElement.Loci {
  const members = link.toMembers(runs, Extents.create(16));
  const byUnit = new Map<number, Map<ChainInstance, readonly [number, number]>>();
  for (let i = 0; i < members.count; i++) {
    const c = link.memberOf(members, i)?.chain;
    if (!c || c.count === 0) continue;
    let chains = byUnit.get(c.unit.id);
    if (!chains) byUnit.set(c.unit.id, (chains = new Map()));
    // the chain's atoms, not first to last element: a ligand's elements are its heavy atoms, its hydrogens follow
    // them (perf-audit-2 10/B5)
    chains.set(c, c.span);
  }
  const elements: { unit: Unit; indices: OrderedSet<StructureElement.UnitIndex> }[] = [];
  for (const unit of link.structure.units) {
    const chains = byUnit.get(unit.id);
    if (!chains) continue;
    const spans = [...chains.values()].sort((a, b) => a[0] - b[0]);
    let indices = OrderedSet.ofBounds(spans[0][0], spans[0][1]);
    for (let k = 1; k < spans.length; k++)
      indices = OrderedSet.union(indices, OrderedSet.ofBounds(spans[k][0], spans[k][1]));
    elements.push({ unit, indices: indices as OrderedSet<StructureElement.UnitIndex> });
  }
  return StructureElement.Loci(link.structure, elements);
}

/**
 * A set's loci in 3D under a reaction: residue by residue within its limits (a non-polymer position its whole residue
 * under `wholeResidue`), else its chains or nothing.
 */
function react(link: StructureLink, r: LinkReaction, size: LociSize | null, read: () => readonly (Loci | null)[]) {
  if (!size) return null;
  const within = size.items <= r.maxItems && size.length <= r.maxLength;
  if (!within && r.over === 'skip') return null;
  const runs = lociRuns(read());
  if (runs.count === 0) return null;
  return within ? link.toLoci(runs, { wholeResidue: r.wholeResidue }) : wholeChains(link, runs);
}

/**
 * Hover and selection both ways between a viewer and a 3D view over one link: Seq* hover lights the residues in 3D and
 * the 3D hover lights Seq*'s (owner `3d`); a selection in either selects in the other, and a 3D selection also focuses
 * the Seq* view on it, a Seq* one the 3D view. A Seq* set is sent by its size (`LinkReaction`, read before it resolves):
 * `DefaultLinkViewsOptions` unless `o` says otherwise. Replays and carries are no one's move
 * (seqspec:interaction#replays-and-carries): linking clears nothing in 3D, and a carry (origin `graph`) is sent only
 * where it changed the selection's 3D image last seen (a 3D selection's own included), never focused. Returns the unlink.
 */
export function linkViews(
  ctx: SeqViewerContext,
  view3d: StructureView3D,
  link: StructureLink,
  o: LinkViewsOptions = {},
  /** False on a frame step (sequence.md §2): the 3D view holds the selection already, so it is not sent again. */
  replay = true,
): () => void {
  const stops: (() => void)[] = [];
  const scratch = Extents.create(16);
  if (o.hover !== false) {
    const r = { ...DefaultLinkViewsOptions.hover, ...o.hover };
    // What was last sent (undefined: not known, the 3D view hovered its own): the same again (a relayout's refresh,
    // the replay at link time) is not sent.
    let lit: StructureElement.Loci | null | undefined = null;
    const sub = ctx.hover.live.version.subscribe(() => {
      if (ctx.hover.live.current.owner === '3d') {
        lit = undefined;
        return;
      }
      const next = react(link, r, ctx.hover.size(), () => [ctx.hover.toLoci() as Loci | null]);
      if (lit !== undefined && sameLoci(next, lit)) return;
      lit = next;
      view3d.highlight(next);
    });
    stops.push(() => sub.unsubscribe());
    stops.push(
      view3d.onHover((loci) => {
        scratch.count = 0;
        ctx.hover.setLoci(loci ? coordinatesOf(link, link.fromLoci(loci, scratch)) : null, '3d');
      }),
    );
  }
  if (o.select !== false) {
    const r = { ...DefaultLinkViewsOptions.select, ...o.select };
    // What was last sent to 3D: Mol*'s event echoing it is not sent back (Seq*'s echo is told by its origin). `seen`:
    // the 3D image of the Seq* selection last seen, the one 3D made too: a carry is compared with it, so a 3D selection
    // Seq* holds only in part (an atom, a unit no system covers) is not replaced by its image on a data edit.
    let sent3d: StructureElement.Loci | null = null;
    let seen: StructureElement.Loci | null = null;
    let replaying = true;
    const sub = ctx.selection.current.subscribe((sel) => {
      const first = replaying;
      replaying = false;
      const parts = [sel.seq, sel.aligned, sel.rows, sel.items, sel.boundaries] as (Loci | null)[];
      const size = { items: 0, length: 0 };
      for (const l of parts) {
        const z = lociSize(l);
        size.items += z?.items ?? 0;
        size.length += z?.length ?? 0;
      }
      const loci = react(link, r, size, () => parts);
      const next = loci?.elements.length ? loci : null;
      // The replay at link time sends a selection Seq* holds, never a clear; a carry only what it changed.
      const same = sameLoci(next, seen);
      seen = next;
      if (first && !replay) {
        // a frame step: the selection's image on the new frame is what 3D holds (perf-audit-2 10/B6)
        sent3d = next;
        return;
      }
      if (sel.origin === 'host:3d' || ((first || sel.origin === 'graph') && same)) return;
      sent3d = next;
      view3d.select(sent3d);
      if (sent3d && o.focus !== false && !first && sel.origin !== 'graph') view3d.focus(sent3d);
    });
    stops.push(() => sub.unsubscribe());
    stops.push(
      view3d.onSelect((loci) => {
        const same = loci === null ? sent3d === null : sent3d !== null && StructureElement.Loci.areEqual(loci, sent3d);
        if (same) return;
        sent3d = loci;
        scratch.count = 0;
        const coords = loci ? coordinatesOf(link, link.fromLoci(loci, scratch)) : null;
        const seq = coords?.runs.length ? coords : null;
        ctx.selection.set(seq, 'host:3d');
        if (seq && o.focus !== false) ctx.focus.fromLoci(seq, { origin: 'host:3d' });
      }),
    );
  }
  return () => {
    for (const s of stops) s();
  };
}
