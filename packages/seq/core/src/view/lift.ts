/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's data carried through a map read by region (_spec/README.md "Resolved view", A liftover; EC-212): data over
// one assembly drawn on a block of another, through a `map` whose runs are a pairs source's rows (`map.source`). No
// static path joins them (the runs arrive by region, EC-211): the layer is placed with the map and the end of its runs
// on the block's system, and the render reads the runs by that end before the data along them. Loaded with the genome
// chunk (`view.ts`): a liftover is over assemblies.

import type { Node } from '@molstar/core/repr/tree';
import type { AccessPlan, StaticModel } from '@molstar/seqx-core/resolve';
import type { LayerData } from './view.js';

/** What a lifted layer reads: the map, the end of its runs on the section's system, and the map's plan there. */
export interface LiftData {
  /** The map read by region (`map.source`) the data is carried through. */
  readonly map: string;
  /** The end of a run on the section's system: `b` on the map's `to` (a liftover's new assembly), `a` on its `from`. */
  readonly side: 'a' | 'b';
  /** The map's plan on the section's system: its source, the name it lists the system by, its side when sided. */
  readonly near: AccessPlan;
}

/** The maps read by region of a spec, by address (`map.source`), kept per spec. */
const RegionMaps = new WeakMap<Node, readonly string[]>();
export function regionMaps(spec: Node): readonly string[] {
  let maps = RegionMaps.get(spec);
  if (maps) return maps;
  const out: string[] = [];
  for (const d of spec.children ?? [])
    if (d.kind === 'data' && d.id !== undefined)
      for (const n of d.children ?? [])
        if (n.kind === 'map' && n.id !== undefined && n.params?.source !== undefined) out.push(`/${d.id}/${n.id}`);
  RegionMaps.set(spec, (maps = out));
  return maps;
}

/**
 * `data` (over an assembly, not reaching `system`) lifted onto `system` through a map read by region that joins the
 * data's assembly to the one `system` is a sequence of; undefined when none does.
 */
export function liftOf(model: StaticModel, data: LayerData, system: string): LayerData | undefined {
  const over = model.assemblies.get(data.system);
  if (!over) return undefined;
  for (const address of regionMaps(model.spec)) {
    const m = model.annotations.get(address);
    if (!m?.bySystem || m.shape !== 'pairs') continue;
    // the data on the map's `from`, the section on its `to` (old genes on the new assembly), or the other way
    const side = m.target === over.address ? 'b' : m.ySystem === over.address ? 'a' : undefined;
    const near = side && m.bySystem.get(system);
    const onOther = side && model.assemblies.get(side === 'b' ? m.ySystem! : m.target)?.systems.includes(system);
    if (!near || !onOther || m.target === m.ySystem) continue;
    return { ...data, system, lift: { map: address, side: side!, near } };
  }
  return undefined;
}
