/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EmptyKeys, location } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { Emphasized } from '@molstar/seqx-viewer/models';
import type { BehaviorContext, BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';

// The `links.hover` behavior (_spec/README.md "Behaviors", Links hover): hovering a
// molecule (a chip, a residue, an atom) lights the polymer positions its links reach; hovering a linked site lights
// the molecule (its whole unit) or, on a link between two polymers (a disulfide), the other end. Written to a
// highlight channel of its own (`linked`), cleared when the hover goes. Off by default: how events from different
// parts compose is the application's; a host lists it (the examples do). A link is never a map: nothing is projected
// across it, the lit ends are its rows' own (a merged copy's drawn on the block it is merged into).

/** The channel the behavior writes: its color is `highlight.colors.linked`, else `highlight.color`. */
export const LinksHoverChannel = 'linked';

export interface LinksHoverBehavior {
  readonly channel: string;
  dispose(): void;
}

/** A system's residues system: an atoms system's, else itself. */
function residuesOf(model: StaticModel, system: string): string {
  return model.molecular.atoms.get(system)?.of ?? system;
}

/** Each merged copy's system, to the system of the block it is merged into (the view's blocks). */
function mergedSystems(view: ResolvedView | undefined): Map<string, string> {
  const out = new Map<string, string>();
  for (const s of view?.sections ?? [])
    for (const b of s.blocks ?? []) for (const sys of b.systems) if (sys !== b.axis.system) out.set(sys, b.axis.system);
  return out;
}

/**
 * The places a hover on `systems` (a block's merged copies too) at `position` lights through the model's links: what
 * the hovered molecule's links reach, or the molecule (or other end) a hovered site's links reach. `[system, start,
 * end]` triples.
 */
export function linkedBy(
  model: StaticModel,
  systems: readonly string[],
  position: number,
): { system: string; start: number; end: number }[] {
  const out: { system: string; start: number; end: number }[] = [];
  const units = model.molecular.units;
  const isMolecule = (s: string) => {
    const u = units.get(s);
    return u !== undefined && u.class !== 'polymer';
  };
  const whole = (s: string) => {
    const u = units.get(s);
    return u?.end === undefined ? 1 : u.end - u.start;
  };
  for (const a of model.annotations.values()) {
    const d = a.pairs;
    if (!d || !a.units || !d.aUnit || !d.bUnit) continue;
    for (let k = 0; k < d.n; k++) {
      const as = a.units[d.aUnit[k]];
      const bs = a.units[d.bUnit[k]];
      for (let side = 0; side < 2; side++) {
        const here = side ? bs : as;
        const there = side ? as : bs;
        if (!systems.includes(here)) continue;
        const s0 = side ? d.bStart[k] : d.aStart[k];
        const e0 = side ? d.bEnd[k] : d.aEnd[k];
        const t0 = side ? d.aStart[k] : d.bStart[k];
        const t1 = side ? d.aEnd[k] : d.bEnd[k];
        // a molecule hovered anywhere lights its links' sites; a site hovered at its end lights the other side
        if (isMolecule(here)) out.push({ system: there, start: t0, end: t1 });
        else if (position >= s0 && position < e0)
          out.push(
            isMolecule(there) ? { system: there, start: 0, end: whole(there) } : { system: there, start: t0, end: t1 },
          );
      }
    }
  }
  return out;
}

export function linksHoverBehavior(ctx: BehaviorContext): LinksHoverBehavior {
  const { hover, highlight, spec, view } = ctx;
  let merged = new Map<string, string>();
  let mergedFor: unknown;
  /** The last places lit (their key) and the value written for them: a hover within them writes nothing new. */
  let lastKey = '';

  const read = () => {
    const r = hover.current;
    const model = spec.model.value;
    if (r.owner === undefined || r.system === undefined || Number.isNaN(r.position)) {
      lastKey = '';
      highlight.clear(LinksHoverChannel);
      return;
    }
    const v = view.resolved.value;
    if (mergedFor !== v) {
      mergedFor = v;
      merged = mergedSystems(v);
    }
    // a merged block stands for every copy: its hover is each copy's
    const section = v?.sections.find((s) => s.address === r.section);
    const block = r.block === undefined ? undefined : section?.blocks?.[r.block];
    const own = residuesOf(model, r.system);
    const systems = block && block.axis.system === r.system ? block.systems.map((s) => residuesOf(model, s)) : [own];
    const lit = linkedBy(model, systems, Math.floor(r.position));
    const parts = new Map<string, number[]>();
    for (const l of lit) {
      const at = merged.get(l.system) ?? l.system;
      let ranges = parts.get(at);
      if (!ranges) parts.set(at, (ranges = []));
      ranges.push(l.start, l.end);
    }
    const key = [...parts].map(([s, r]) => `${s}:${r.join(',')}`).join(';');
    if (key === lastKey) return;
    lastKey = key;
    if (!parts.size) {
      highlight.clear(LinksHoverChannel);
      return;
    }
    const loc = location({ parts: [...parts].map(([system, ranges]) => ({ system, ranges })) });
    const value: Emphasized = { elements: EmptyKeys, location: loc, shown: loc };
    highlight.set(LinksHoverChannel, value);
  };

  const stops = [
    hover.version.listen(read),
    view.resolved.listen(() => {
      lastKey = '';
      read();
    }),
    spec.model.listen(() => {
      lastKey = '';
      read();
    }),
  ];
  read();
  return {
    channel: LinksHoverChannel,
    dispose() {
      for (const s of stops) s();
      stops.length = 0;
      highlight.clear(LinksHoverChannel);
    },
  };
}

/** The `links.hover` behavior's provider: not in `DefaultSeqViewerSpec()`; a host (an example) lists it. */
export const LinksHoverProvider: BehaviorProvider<undefined, LinksHoverBehavior> = {
  name: 'links.hover',
  create: (ctx) => linksHoverBehavior(ctx),
};
