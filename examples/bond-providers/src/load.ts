/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Task } from '@molstar/core/task';
import type { StateObjectSelector } from '@molstar/core/state';
import { CCDFormat } from '@molstar/model/formats/structure/mmcif';
import { ComponentBond } from '@molstar/model/formats/structure/property/bonds/chem_comp';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import type { Model } from '@molstar/model/model/structure/model';
import type { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { TrajectoryFromMmCif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { ParseCifForPerception } from './removals.js';

// Loading the two sides (the bond-order benchmark's plan §7): the file is downloaded once and its bytes given to both
// viewers. Left, the file as served with its built-in bonds; right, the same file through the removals, with the bonds
// of the provider under test.

/** A PDB entry (`models.rcsb.org`'s BinaryCIF) or a CCD component at its ideal or model coordinates. */
export interface Source {
  kind: 'pdb' | 'ccd';
  id: string;
  coords: CCDFormat.CoordinateType;
}

/** A bond provider and its params (missing ones get defaults). */
export interface Bonds {
  provider: string;
  params: Record<string, unknown>;
}

export interface Removals {
  removeHydrogens: boolean;
  removeCharges: boolean;
}

/** The right side's state cells: the parse (its removals updated in place) and the trajectory. */
export interface Tested {
  parse: StateObjectSelector<SO.Format.Cif>;
  trajectory: StateObjectSelector<SO.Molecule.Trajectory>;
  modelIndex: number;
}

/** The file's bytes, BinaryCIF for an entry, mmCIF text for a component. */
export async function download(source: Source): Promise<Uint8Array<ArrayBuffer> | string> {
  const url =
    source.kind === 'pdb'
      ? `https://models.rcsb.org/${source.id.toLowerCase()}.bcif`
      : `https://files.rcsb.org/ligands/download/${source.id}.cif`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${source.id}: ${response.status} ${response.statusText} (${url})`);
  return source.kind === 'pdb' ? new Uint8Array(await response.arrayBuffer()) : response.text();
}

/** Left: the file as served, its built-in bonds. */
export async function showServed(viewer: PrismViewer, source: Source, data: Uint8Array<ArrayBuffer> | string) {
  await viewer.api.session.clear();
  await viewer.run(`Load ${source.id}`, async (plugin) => {
    const raw = await plugin.builders.data.rawData({ data, label: source.id }, { state: { isGhost: true } });
    const trajectory = await plugin.builders.structure.parseTrajectory(raw, 'mmcif');
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
      model: { modelIndex: await modelIndexOf(trajectory, source) },
      representationPreset: representationOf(source),
    });
  });
}

/** Right: the file through the removals, with `bonds`. */
export async function showTested(
  viewer: PrismViewer,
  source: Source,
  data: Uint8Array<ArrayBuffer> | string,
  removals: Removals,
  bonds: Bonds,
): Promise<Tested> {
  await viewer.api.session.clear();
  let tested: Tested | undefined;
  await viewer.run(`Load ${source.id}`, async (plugin) => {
    const raw = await plugin.builders.data.rawData({ data, label: source.id }, { state: { isGhost: true } });
    const parse = await plugin.state.data
      .build()
      .to(raw)
      .apply(ParseCifForPerception, removals, { state: { isGhost: true } })
      .commit({ revertOnError: true });
    const trajectory = await plugin.state.data
      .build()
      .to(parse)
      .apply(TrajectoryFromMmCif)
      .commit({ revertOnError: true });
    const modelIndex = await modelIndexOf(trajectory, source);
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
      model: { modelIndex },
      structure: { name: 'model', params: { bonds: { name: bonds.provider, params: bonds.params } } },
      representationPreset: representationOf(source),
    });
    tested = { parse, trajectory, modelIndex };
  });
  return tested!;
}

/** Applies other removals to the right side: parsed again, what depends on it updated (the camera kept). */
export async function updateRemovals(viewer: PrismViewer, tested: Tested, removals: Removals) {
  await viewer.run('Removals', async (plugin) => {
    await plugin.state.data.build().to(tested.parse).update(removals).commit({ revertOnError: true });
  });
}

/**
 * Whether the right side's model states its bond orders, so that `auto` keeps them (plan §7): a component bond table
 * (`chem_comp_bond`) or an index-pair one (`perceive.ts`'s test).
 */
export async function statesOrders(tested: Tested) {
  const t = tested.trajectory.cell?.obj?.data;
  if (!t) return false;
  const model: Model = await Task.resolveInContext(t.getFrameAtIndex(tested.modelIndex));
  return !!(ComponentBond.Provider.get(model) || IndexPairBonds.Provider.get(model));
}

function representationOf(source: Source) {
  return source.kind === 'pdb' ? 'polymer-and-ligand' : 'atomic-detail';
}

/**
 * An entry's first model; a CCD file's models are its ideal coordinates then its model ones, each when present
 * (`trajectoryFromCCD`).
 */
async function modelIndexOf(trajectory: StateObjectSelector<SO.Molecule.Trajectory>, source: Source) {
  const t = trajectory.cell?.obj?.data;
  if (!t || source.kind === 'pdb') return 0;
  for (let i = 0; i < t.frameCount; i++) {
    const model: Model = await Task.resolveInContext(t.getFrameAtIndex(i));
    if (CCDFormat.CoordinateType.get(model) === source.coords) return i;
  }
  throw new Error(`${source.id} has no ${source.coords} coordinates.`);
}
