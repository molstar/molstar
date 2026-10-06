/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import { setFSModule } from '@molstar/core/util/data-source';
import { PluginContext } from '@molstar/plugin/context';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { PluginSpec, type PluginRegistryEntry } from '@molstar/plugin/spec';
import { Ccp4 } from '@molstar/plugin/state/formats/volume/ccp4';
import { Dx } from '@molstar/plugin/state/formats/volume/dx';
import { Dcd } from '@molstar/plugin/state/formats/coordinates/dcd';
import { LammpsTrajectory } from '@molstar/plugin/state/formats/coordinates/lammps';
import { Nctraj } from '@molstar/plugin/state/formats/coordinates/nctraj';
import { Trr } from '@molstar/plugin/state/formats/coordinates/trr';
import { Xtc } from '@molstar/plugin/state/formats/coordinates/xtc';
import { Obj } from '@molstar/plugin/state/formats/shape/obj';
import { Ply } from '@molstar/plugin/state/formats/shape/ply';
import { Vtp } from '@molstar/plugin/state/formats/shape/vtp';
import { Prmtop } from '@molstar/plugin/state/formats/topology/prmtop';
import { Psf } from '@molstar/plugin/state/formats/topology/psf';
import { Top } from '@molstar/plugin/state/formats/topology/top';
import { Gro } from '@molstar/plugin/state/formats/trajectory/gro';
import { Mmcif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import { Mol } from '@molstar/plugin/state/formats/trajectory/mol';
import { Mol2 } from '@molstar/plugin/state/formats/trajectory/mol2';
import { Pdb, Pdbqt } from '@molstar/plugin/state/formats/trajectory/pdb';
import { Sdf } from '@molstar/plugin/state/formats/trajectory/sdf';
import { Xyz } from '@molstar/plugin/state/formats/trajectory/xyz';
import { BuiltInMarkdownExtension } from '@molstar/plugin/state/markdown/catalog';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { loadMVS } from '@molstar/mvs/load';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';

setFSModule(fs);

/** The format entries for every `parse`, `coordinates`, `topology`, and `trajectory` format the MVS loader supports. */
const MVSFormats: PluginRegistryEntry[] = [
  Mmcif, // mmcif, bcif, cif
  Pdb,
  Pdbqt,
  Gro,
  Xyz,
  Mol,
  Sdf,
  Mol2,
  Xtc,
  Trr,
  Dcd,
  Nctraj,
  LammpsTrajectory, // lammpstrj
  Psf,
  Prmtop,
  Top,
  Ccp4, // map
  Dx,
  Vtp,
  Ply,
  Obj,
];

/** The markdown commands of snapshot descriptions: MVS documents can use any of the built-in extensions. */
const MVSMarkdownExtensions: PluginRegistryEntry = { markdownExtensions: BuiltInMarkdownExtension };

/**
 * `examples/mvs/1cbs.mvsj` with its structure read from the local copy in `examples`. The label nodes are removed:
 * labels render text through the canvas API, which Node does not have.
 */
function createLocalDocument() {
  const mvsj = JSON.parse(fs.readFileSync('examples/mvs/1cbs.mvsj', 'utf8'));
  const visit = (node: any) => {
    if (node.kind === 'download') node.params.url = `file://${process.cwd()}/examples/1cbs_full.bcif`;
    node.children = node.children?.filter((c: any) => c.kind !== 'label');
    node.children?.forEach(visit);
  };
  visit(mvsj.root);
  return MVSData.fromMVSJ(JSON.stringify(mvsj));
}

async function load(registry: PluginRegistryEntry[]) {
  const plugin = new PluginContext({ registry, behaviors: [PluginSpec.Behavior(MolViewSpec)] });
  await plugin.init();
  const warn = jest.spyOn(console, 'warn').mockImplementation(() => {});
  const error = jest.spyOn(console, 'error').mockImplementation(() => {});
  let consoleWarnings: unknown[][] = [];
  let consoleErrors: unknown[][] = [];
  try {
    // `keepCamera`: the document has a camera node and the plugin has no canvas.
    await loadMVS(plugin, createLocalDocument(), { sanityChecks: true, keepCamera: true });
    consoleWarnings = [...warn.mock.calls];
    consoleErrors = [...error.mock.calls];
  } finally {
    warn.mockRestore();
    error.mockRestore();
  }
  const cells = [...plugin.state.data.cells.values()];
  return {
    plugin,
    consoleWarnings,
    consoleErrors,
    problems: plugin.log.entries.toArray().filter((e) => e.type === 'warning' || e.type === 'error'),
    cellErrors: cells.filter((c) => c.status === 'error').map((c) => c.errorText),
    representations: cells
      .filter((c) => c.transform.transformer.id === 'ms-plugin.structure-representation-3d')
      .map((c) => ({
        status: c.status,
        type: (c.transform.params as any).type.name,
        color: (c.transform.params as any).colorTheme.name,
        size: (c.transform.params as any).sizeTheme.name,
      })),
  };
}

describe('MVS-only plugin', () => {
  const expectedRepresentations = [
    { status: 'ok', type: 'cartoon', color: 'mvs-multilayer', size: 'uniform' },
    { status: 'ok', type: 'ball-and-stick', color: 'uniform', size: 'uniform' },
  ];

  it('loads an MVSJ document with the MVS registry, the MVS formats, and the markdown extensions', async () => {
    const result = await load([MVSRuntimeRegistry, ...MVSFormats, MVSMarkdownExtensions]);
    expect(result.consoleWarnings).toEqual([]);
    expect(result.consoleErrors).toEqual([]);
    expect(result.problems).toEqual([]);
    expect(result.cellErrors).toEqual([]);
    // The representations are the requested ones: an unregistered name would be replaced by the registry default.
    expect(result.representations).toEqual(expectedRepresentations);
  });

  it('does not need format entries or markdown extensions to load: the loader names transformers directly', async () => {
    const result = await load([MVSRuntimeRegistry]);
    expect(result.consoleWarnings).toEqual([]);
    expect(result.consoleErrors).toEqual([]);
    expect(result.problems).toEqual([]);
    expect(result.cellErrors).toEqual([]);
    expect(result.representations).toEqual(expectedRepresentations);
    // Except the formats and the drag-and-drop handler of MVS itself, which its behavior registers.
    expect(result.plugin.dataFormats.list.map((f) => f.name)).toEqual(['MVSJ', 'MVSX']);
  });

  it('without MVSRuntimeRegistry the requested representations are replaced by the registry default', async () => {
    const result = await load([]);
    expect(result.representations.map((r) => r.type)).not.toContain('cartoon');
    expect(result.representations.map((r) => r.type)).not.toContain('ball-and-stick');
  });
});
