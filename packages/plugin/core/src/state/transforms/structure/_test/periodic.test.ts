/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Model, type Structure } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import {
  CoordinatesWithBoxes,
  ModelFromTrajectory,
  StructureFromTrajectory,
  TrajectoryFromModelAndCoordinates,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { StateSelection } from '@molstar/core/state';
import { getFileNameInfo } from '@molstar/core/util/file-info';
import { PolymerCenter, type PeriodicProps } from '@molstar/plugin/state/helpers/periodic';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { createPlaybackPlugin } from '../../../animation/_test/harness.js';

// The fixture: data/examples/trajectory/periodic/split-peptide.gro, a peptide split across a dodecahedral box in its
// first frame (between residues 2 and 3, a false polymer gap where it was read).
const gro = fs.readFileSync(
  path.resolve(import.meta.dirname, '../../../../../../../../data/examples/trajectory/periodic/split-peptide.gro'),
  'utf8',
);

async function load(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: gro, label: 'split-peptide.gro' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'gro');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { model, structure };
}

const setModel = (plugin: PluginContext, ref: string, params: { modelIndex: number } & PeriodicProps) =>
  plugin
    .build()
    .to(ref)
    .update(ModelFromTrajectory, () => params)
    .commit();

/** The CA–CA distance of residues 2 and 3, across the face in the first frame as read. */
function caca(model: Model) {
  const { x, y, z } = model.atomicConformation;
  // five atoms a residue: N CA C O CB
  return Math.hypot(x[6] - x[11], y[6] - y[11], z[6] - z[11]);
}

const gaps = (s: Structure) => s.units.reduce((n, u) => n + u.gapElements.length, 0);

describe('ModelFromTrajectory in a periodic box', () => {
  it('makes the molecules whole by default, and shows the file as read when asked', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model, structure } = await load(plugin);
    const cell = () => plugin.state.data.cells.get(model.ref)!.obj!.data as Model;
    const units = () => plugin.state.data.cells.get(structure.ref)!.obj!.data as Structure;

    expect(caca(cell())).toBeLessThan(4);
    expect(gaps(units())).toBe(0);

    await setModel(plugin, model.ref, { modelIndex: 0, periodic: 'off' });
    expect(caca(cell())).toBeGreaterThan(10);
    // the structure was remapped onto the frame as read: its false gap is back
    expect(gaps(units())).toBeGreaterThan(0);

    await setModel(plugin, model.ref, { modelIndex: 0, periodic: 'whole' });
    expect(caca(cell())).toBeLessThan(4);
    expect(gaps(units())).toBe(0);
  });

  it('keeps every frame whole, the same model for the same frame', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await load(plugin);
    const first = plugin.state.data.cells.get(model.ref)!.obj!.data as Model;
    for (let i = 1; i < 8; i++) {
      // as the animation updates it: the model index only
      await setModel(plugin, model.ref, { modelIndex: i });
      expect(caca(plugin.state.data.cells.get(model.ref)!.obj!.data as Model)).toBeLessThan(4);
    }
    await setModel(plugin, model.ref, { modelIndex: 0 });
    expect(plugin.state.data.cells.get(model.ref)!.obj!.data).toBe(first);
  });

  it('moves the molecules into the box, around the polymers when asked', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await load(plugin);
    const cell = () => plugin.state.data.cells.get(model.ref)!.obj!.data as Model;
    /** The peptide's mean (its first 20 atoms) in fractional coordinates of the box. */
    const peptide = (m: Model) => {
      const { x, y, z } = m.atomicConformation;
      const c = Vec3();
      for (let i = 0; i < 20; i++) Vec3.add(c, c, Vec3.create(x[i], y[i], z[i]));
      return Vec3.transformMat3(c, Vec3.scale(c, c, 1 / 20), Model.PeriodicBox.get(m)!.inverse);
    };

    for (let i = 0; i < 8; i++) {
      await setModel(plugin, model.ref, { modelIndex: i, periodic: 'box' });
      expect(caca(cell())).toBeLessThan(4);
      for (const v of peptide(cell())) expect(v >= 0 && v < 1).toBe(true);

      await setModel(plugin, model.ref, {
        modelIndex: i,
        periodic: 'box',
        periodicCenter: { name: 'selection', params: { script: PolymerCenter } },
      });
      expect(caca(cell())).toBeLessThan(4);
      for (const v of peptide(cell())) expect(v).toBeCloseTo(0.5, 4);
    }
  });

  it('keeps the molecules from jumping as the frames play', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await load(plugin);
    const center = () => {
      const { x, y, z } = (plugin.state.data.cells.get(model.ref)!.obj!.data as Model).atomicConformation;
      const c = Vec3();
      for (let i = 0; i < 20; i++) Vec3.add(c, c, Vec3.create(x[i], y[i], z[i]));
      return Vec3.scale(c, c, 1 / 20);
    };
    let previous: Vec3 | undefined;
    for (const i of [0, 1, 2, 3, 4, 5, 6, 7, 6, 5]) {
      await setModel(plugin, model.ref, { modelIndex: i, periodic: 'nojump' });
      expect(caca(plugin.state.data.cells.get(model.ref)!.obj!.data as Model)).toBeLessThan(4);
      const c = center();
      if (previous) expect(Vec3.distance(c, previous)).toBeLessThan(8);
      previous = c;
    }
  });

  it('builds an ensemble of all frames from whole frames, or from the file as read', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const data = await plugin.builders.data.rawData({ data: gro, label: 'split-peptide.gro' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'gro');
    /** The longest CA–CA distance of residues 2 and 3 over the ensemble's frames. */
    const longest = async (periodic: 'whole' | 'off') => {
      const cell = await plugin.build().to(trajectory).apply(StructureFromTrajectory, { periodic }).commit();
      const s = cell.obj!.data as Structure;
      let max = 0;
      for (const u of s.units) {
        if (u.elements.length < 20) continue;
        const a = Vec3(),
          b = Vec3();
        u.conformation.position(u.elements[6], a);
        u.conformation.position(u.elements[11], b);
        max = Math.max(max, Vec3.distance(a, b));
      }
      await plugin.build().delete(cell.ref).commit();
      return max;
    };
    expect(await longest('whole')).toBeLessThan(4);
    expect(await longest('off')).toBeGreaterThan(10);
  });

  // Unit.remapModel places a frame's copies from its cell when its symmetry differs: the images follow an NPT box
  it("places the periodic images from each frame's box (NPT: the box changes)", async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    // two atoms, a 20 Å box, then 22 Å
    const frame = (t: number, l: number) =>
      `ITEM: TIMESTEP\n${t}\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n0 ${l}\n0 ${l}\n0 ${l}\nITEM: ATOMS id mol type x y z\n1 1 1 5 5 5\n2 1 1 6 5 5\n`;
    const data = await plugin.builders.data.rawData({ data: frame(0, 20) + frame(1, 22), label: 'npt.lammpstrj' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'lammps_traj_data');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model, {
      name: 'symmetry',
      params: { ijkMin: Vec3.create(-1, -1, -1), ijkMax: Vec3.create(1, 1, 1) },
    });
    /** The largest x translation of the copies: the box's a. */
    const a = () =>
      Math.max(
        ...(plugin.state.data.cells.get(structure.ref)!.obj!.data as Structure).units.map(
          (u) => u.conformation.operator.matrix[12],
        ),
      );
    expect((plugin.state.data.cells.get(structure.ref)!.obj!.data as Structure).units.length).toBe(27);
    expect(a()).toBeCloseTo(20, 3);
    await setModel(plugin, model.ref, { modelIndex: 1 });
    expect(a()).toBeCloseTo(22, 3);
  });

  it('gives a dropped NAMD box file to the last trajectory from coordinates, replacing boxes it had', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const example = (name: string) =>
      fs.readFileSync(path.resolve(import.meta.dirname, `../../../../../../../../data/examples/trajectory/${name}`));
    const pdb = await plugin.builders.data.rawData({
      data: example('protein.pdb').toString('utf8'),
      label: 'protein.pdb',
    });
    const model = await plugin.builders.structure.createModel(
      await plugin.builders.structure.parseTrajectory(pdb, 'pdb'),
    );
    const dcd = await plugin.builders.data.rawData({
      data: new Uint8Array(example('protein.dcd')),
      label: 'protein.dcd',
    });
    const coordinates = await plugin.dataFormats.get('dcd')!.parse(plugin, dcd);
    const trajectory = await plugin
      .build()
      .toRoot()
      .apply(TrajectoryFromModelAndCoordinates, { modelRef: model.ref, coordinatesRef: coordinates.ref })
      .commit();
    const frame = await plugin.build().to(trajectory).apply(ModelFromTrajectory, { modelIndex: 0 }).commit();

    const drop = async (text: string, name: string) => {
      const data = await plugin.builders.data.rawData({ data: text, label: name });
      const provider = plugin.dataFormats.auto(getFileNameInfo(name), data.cell!.obj!)!;
      expect(provider.name).toBe('xsc');
      await provider.visuals!(plugin, await provider.parse(plugin, data));
    };
    const origin = () =>
      Array.from(Model.PeriodicBox.get(plugin.state.data.cells.get(frame.ref)!.obj!.data as Model)!.origin);

    await drop('0 40 0 0 0 40 0 0 0 40 5 -3 2\n', 'a.xsc');
    expect(origin()).toEqual([5 - 20, -3 - 20, 2 - 20]);
    await drop('0 50 0 0 0 50 0 0 0 50 0 0 0\n', 'b.xsc');
    expect(origin()).toEqual([-25, -25, -25]);
    // one boxes step, the second file's
    expect(plugin.state.data.select(StateSelection.Generators.ofTransformer(CoordinatesWithBoxes))).toHaveLength(1);
  });
});
