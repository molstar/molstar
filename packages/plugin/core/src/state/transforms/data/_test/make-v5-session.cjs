/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Writes `v5-session.molx`, the fixture of `v5-session.test.ts`: a session saved by Mol* 5.x with its downloads and an
// opened file as assets (a few KB). Run by hand with a built 5.x checkout:
//   node packages/plugin/core/src/state/transforms/data/_test/make-v5-session.cjs <path to molstar 5.x>
// The inputs are generated here: a 2-atom mmCIF (as text, and gzipped as binary through `DeflateData`, the EMDB
// loader's way), a 1-atom PDB, and the mmCIF opened as a file.

const fs = require('node:fs');
const path = require('node:path');
const zlib = require('node:zlib');

const v5 = path.resolve(process.argv[2] ?? '../molstar-v5');
const lib = (p) => require(path.join(v5, 'lib/commonjs', p));
const { PluginContext } = lib('mol-plugin/context');
const { DefaultPluginSpec } = lib('mol-plugin/spec');
const { Download, DeflateData, ReadFile } = lib('mol-plugin-state/transforms/data');
const { Asset } = lib('mol-util/assets');
const { File_ } = lib('mol-util/nodejs-shims');
const version = require(path.join(v5, 'package.json')).version;

const Cif = `data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ALA A 1 1.5 2.5 3.5 A 1 1
ATOM 2 N N ALA A 1 0.5 0.5 0.5 A 1 1
`;
const Pdb = 'ATOM      1  CA  ALA A   1       1.500   2.500   3.500  1.00  0.00           C\nEND\n';

async function main() {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  // the downloads' bytes held by the asset manager, as after a download: nothing reaches the network
  const held = (url, data) => {
    const asset = Asset.Url(url);
    plugin.managers.asset.set(asset, new File_([data], 'raw-data'));
    return asset;
  };
  const cifUrl = held('https://example.org/a.cif', Cif);
  const pdbUrl = held('https://example.org/b.pdb', Pdb);
  const gzUrl = held('https://example.org/c.cif.gz', zlib.gzipSync(Cif));
  const b = plugin.builders;
  const text = await b.data.download({ url: cifUrl, isBinary: false });
  await b.structure.parseTrajectory(text, 'mmcif');
  const pdb = await b.data.download({ url: pdbUrl, isBinary: false });
  await b.structure.parseTrajectory(pdb, 'pdb');
  const gz = await plugin
    .build()
    .toRoot()
    .apply(Download, { url: gzUrl, isBinary: true })
    // 5.x read binary as BinaryCIF: a gzipped text CIF is inflated to a string
    .apply(DeflateData, { isString: true })
    .commit();
  await b.structure.parseTrajectory(gz, 'mmcif');
  // 5.x's own `File` (its asset manager checks `instanceof` it)
  const file = Asset.File(new File_([Cif], 'd.cif'));
  const read = await plugin.build().toRoot().apply(ReadFile, { file, isBinary: false }).commit();
  await b.structure.parseTrajectory(read, 'mmcif');

  const params = { camera: false, canvas3d: false, canvas3dContext: false, animation: false, startAnimation: false };
  const session = await plugin.managers.snapshot.serialize({ type: 'zip', params });
  const out = path.join(__dirname, 'v5-session.molx');
  fs.writeFileSync(out, Buffer.from(await session.arrayBuffer()));
  console.log(`Mol* ${version}: wrote ${out} (${fs.statSync(out).size} bytes)`);
  plugin.dispose();
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
