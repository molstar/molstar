import assert from 'node:assert/strict';
import { Color } from '@molstar/core/util/color';
import { Task } from '@molstar/core/task';
import { CIF } from '@molstar/io/reader/cif';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/spec';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { unzipSync } from 'fflate';
import fs from 'node:fs';

assert.equal(typeof Color(0xff0000), 'number');
const taskResult = await Task.create('smoke task', async () => 42).run();
assert.equal(taskResult, 42);

const cifText = `data_smoke\n_entry.id smoke\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA ALA A 1 0 0 0 1 10 1\n#`;
const parsed = await CIF.parse(cifText).run();
assert.equal(parsed.isError, false, parsed.message);
assert.equal(parsed.result.blocks[0].header, 'smoke');

const fixturePdb = fs.readFileSync('tiny.pdb', 'utf8');
const plugin = new PluginContext(DefaultPluginSpec());
try {
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: fixturePdb, label: 'local smoke PDB' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  assert(structure.cell.obj.data.elementCount > 0, 'PluginContext parsed no atoms from the local PDB');
} finally {
  plugin.dispose();
}

const mvsData = MVSData.fromMVSJ(fs.readFileSync('tiny.mvsj', 'utf8'));
assert.equal(MVSData.validationIssues(mvsData), undefined, 'MVS builder rejected the local smoke document');
const mvsx = await MVSData.toMVSX(mvsData, { assets: { '/fixtures/tiny.pdb': fixturePdb } });
const archive = unzipSync(mvsx);
assert(archive['index.mvsj'], 'MVSX archive has no index.mvsj');
const archivedMvs = MVSData.fromMVSJ(new TextDecoder().decode(archive['index.mvsj']));
assert.equal(MVSData.validationIssues(archivedMvs), undefined, 'MVSX round trip produced invalid MVS data');
assert.match(
  archivedMvs.root.children[0].params.url,
  /^\.\/assets\//,
  'MVSX did not rewrite its local structure asset',
);
assert(
  Object.keys(archive).some((name) => name.startsWith('./assets/')),
  'MVSX archive has no local structure asset',
);
console.log('Node ESM plugin and MVS consumers passed');
