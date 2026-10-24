/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { StateTransformer } from '@molstar/core/state';
import { RuntimeContext } from '@molstar/core/task';
import { Asset } from '@molstar/core/util/assets';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultFormats } from '@molstar/plugin/default-registry';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParseCif } from '../../../formats/cif.js';
import { DownloadFile, OpenFiles } from '../../../actions/file.js';
import { gzipSync } from 'node:zlib';
import { ParseJson } from '../json.js';
import { DeflateData, Download, ReadFile, ReadResponse } from '../fetch.js';
import { guessCifVariant } from '../../../formats/provider.js';
import { looksLikeMmcifParticles } from '@molstar/model/formats/particles/mmcif-sniff';
import { getFileNameInfo } from '@molstar/core/util/file-info';

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

async function createPlugin() {
  const plugin = new PluginContext({ behaviors: [], registry: [DefaultFormats] });
  await plugin.init();
  return plugin;
}

const cifFile = (name = 'x.cif') => Asset.File(new File([Cif], name));

describe('Data.Response (.v6/data-response)', () => {
  it('ReadFile with asResponse holds the bytes, shared with the asset manager', async () => {
    const plugin = await createPlugin();
    const file = cifFile();
    const { data } = await plugin.builders.data.readFile({ file, asResponse: true });
    const obj = data.cell!.obj!;
    expect(obj).toBeInstanceOf(SO.Data.Response);
    expect(obj.label).toBe('x.cif');
    const held = await plugin.managers.asset.resolve(file, 'binary', { shared: true }).run();
    expect(held.data).toBe(obj.data.bytes);
    held.dispose();
    expect(`${await obj.data.text(RuntimeContext.Synchronous)}`).toBe(Cif);
    // reads repeat, in any form
    expect(`${await obj.data.read(RuntimeContext.Synchronous, 'string')}`).toBe(Cif);

    await plugin.build().delete(data).commit();
    expect(plugin.managers.asset.has(file)).toBe(false);
    plugin.dispose();
  });

  it('ParseCif reads a response as it reads the string', async () => {
    const plugin = await createPlugin();
    const { data } = await plugin.builders.data.readFile({ file: cifFile(), asResponse: true });
    const fromResponse = await plugin.build().to(data).apply(ParseCif).commit();
    const raw = await plugin.builders.data.rawData({ data: Cif });
    const fromString = await plugin.build().to(raw).apply(ParseCif).commit();
    const a = fromResponse.cell!.obj!.data.blocks[0].categories.atom_site;
    const b = fromString.cell!.obj!.data.blocks[0].categories.atom_site;
    expect(a.rowCount).toBe(2);
    for (const name of b.fieldNames)
      expect(a.getField(name)!.toStringArray()).toEqual(b.getField(name)!.toStringArray());
    plugin.dispose();
  });

  it('ParseJson reads a response', async () => {
    const plugin = await createPlugin();
    const file = Asset.File(new File(['{"a":[1,2]}'], 'x.json'));
    const { data } = await plugin.builders.data.readFile({ file, asResponse: true });
    const json = await plugin.build().to(data).apply(ParseJson).commit();
    expect(json.cell!.obj!.data).toEqual({ a: [1, 2] });
    plugin.dispose();
  });

  it('ReadResponse gives String Data, or Binary Data as a copy', async () => {
    const plugin = await createPlugin();
    const { data } = await plugin.builders.data.readFile({ file: cifFile(), asResponse: true });
    const text = await plugin.build().to(data).apply(ReadResponse).commit();
    expect(text.cell!.obj).toBeInstanceOf(SO.Data.String);
    expect(`${text.cell!.obj!.data}`).toBe(Cif);
    const bytes = await plugin.build().to(data).apply(ReadResponse, { isBinary: true }).commit();
    expect(bytes.cell!.obj).toBeInstanceOf(SO.Data.Binary);
    expect(bytes.cell!.obj!.data).not.toBe(data.cell!.obj!.data.bytes);
    expect(Array.from(bytes.cell!.obj!.data as Uint8Array)).toEqual(Array.from(data.cell!.obj!.data.bytes));
    plugin.dispose();
  });

  it('without asResponse nothing changes', async () => {
    const plugin = await createPlugin();
    const { data } = await plugin.builders.data.readFile({ file: cifFile() });
    expect(data.cell!.obj).toBeInstanceOf(SO.Data.String);
    plugin.dispose();
  });

  it('Download and ReadFile recreate when asResponse or isBinary change', () => {
    for (const t of [Download, ReadFile]) {
      const update = (oldParams: any, newParams: any) =>
        t.definition.update!(
          { oldParams, newParams, a: undefined as any, b: undefined as any, cache: {} } as any,
          undefined as any,
        );
      expect(update({ asResponse: false }, { asResponse: true })).toBe(StateTransformer.UpdateResult.Recreate);
      expect(update({ isBinary: false }, { isBinary: true })).toBe(StateTransformer.UpdateResult.Recreate);
    }
  });

  it('OpenFiles reads CIF as a response and the other formats as before', async () => {
    const plugin = await createPlugin();
    expect(plugin.dataFormats.responseExtensions.has('cif')).toBe(true);
    expect(plugin.dataFormats.responseExtensions.has('bcif')).toBe(true);
    expect(plugin.dataFormats.responseExtensions.has('pdb')).toBe(false);

    const pdb = 'ATOM      1  CA  ALA A   1       1.500   2.500   3.500  1.00  0.00           C\nEND\n';
    const files = [cifFile('a.cif'), Asset.File(new File([pdb], 'b.pdb'))];
    await plugin.runTask(
      plugin.state.data.applyAction(OpenFiles, { files, format: { name: 'auto', params: {} }, visuals: false }),
    );
    const types = new Map<string, string>();
    plugin.state.data.cells.forEach((cell) => {
      const t = cell.transform.transformer;
      if (t === ReadFile) types.set(cell.obj!.label, cell.obj!.type.name);
    });
    expect(types.get('a.cif')).toBe('Response Data');
    expect(types.get('b.pdb')).toBe('String Data');
    let trajectories = 0;
    plugin.state.data.cells.forEach((cell) => {
      if (SO.Molecule.Trajectory.is(cell.obj)) trajectories++;
    });
    expect(trajectories).toBe(2);
    plugin.dispose();
  });

  it('the CIF sniffers tell text from its bytes as from the string', () => {
    const info = getFileNameInfo('x.cif');
    const texts = [
      Cif,
      'data_SERVER\n#\n_density_server_result.x 1\n',
      'data_SERVER\n#\ndata_SEGMENTATION_DATA\n',
      'data_x\n_volume_data_3d_info.name a\n_segmentation_data_3d.v 1\n',
      'data_x\n_refln.pdbx_FWT 1\n',
      'data_x\n_atom_site_fract_x 0.5\n',
      'data_cellPACK\n_atom_site.id 1\n',
      'data_x\n_pdbx_struct_assembly_gen.PDB_model_num 1\n',
    ];
    const enc = new TextEncoder();
    for (const t of texts) {
      expect(guessCifVariant(info, enc.encode(t))).toBe(guessCifVariant(info, t));
      expect(looksLikeMmcifParticles(info, enc.encode(t))).toBe(looksLikeMmcifParticles(info, t));
    }
  });
});

describe('decompress (util/_spec/assets.md §5)', () => {
  const gz = (text: string) => new Uint8Array(gzipSync(new TextEncoder().encode(text)));

  it('ReadFile inflates gzip by its bytes unless decompress is none', async () => {
    const plugin = await createPlugin();
    const file = () => Asset.File(new File([gz(Cif)], 'x.cif'));
    const text = (await plugin.builders.data.readFile({ file: file() })).data;
    expect(`${text.cell!.obj!.data}`).toBe(Cif);
    const response = (await plugin.builders.data.readFile({ file: file(), asResponse: true })).data;
    expect(`${await response.cell!.obj!.data.text(RuntimeContext.Synchronous)}`).toBe(Cif);
    const raw = (await plugin.builders.data.readFile({ file: file(), isBinary: true, decompress: 'none' })).data;
    expect((raw.cell!.obj!.data as Uint8Array)[0]).toBe(0x1f);
    plugin.dispose();
  });

  it('a Download of a held gzipped URL inflates it, and DeflateData after it passes the bytes through', async () => {
    const plugin = await createPlugin();
    const url = 'https://x.org/emd_1.map.gz';
    plugin.managers.asset.set(Asset.Url(url), new File([gz('map bytes')], 'emd_1.map.gz'), { isStatic: true });
    const data = await plugin.build().toRoot().apply(Download, { url, isBinary: true }).apply(DeflateData).commit();
    expect(new TextDecoder().decode(data.cell!.obj!.data as Uint8Array)).toBe('map bytes');
    plugin.dispose();
  });

  it('DownloadFile with the gzip format reads the inflated file by its inner name', async () => {
    const plugin = await createPlugin();
    const url = 'https://x.org/a.cif.gz';
    // a held asset: the download reads it, nothing reaches the network
    plugin.managers.asset.set(Asset.Url(url), new File([gz(Cif)], 'raw-data'), { isStatic: true });
    await plugin.runTask(
      plugin.state.data.applyAction(DownloadFile, { url, format: 'gzip', isBinary: true, visuals: false }),
    );
    let trajectories = 0;
    plugin.state.data.cells.forEach((cell) => {
      if (SO.Molecule.Trajectory.is(cell.obj)) trajectories++;
    });
    expect(trajectories).toBe(1);
    plugin.dispose();
  });

  it('Download and ReadFile recreate when decompress changes', () => {
    for (const t of [Download, ReadFile]) {
      const update = (oldParams: any, newParams: any) =>
        t.definition.update!(
          { oldParams, newParams, a: undefined as any, b: undefined as any, cache: {} } as any,
          undefined as any,
        );
      expect(update({}, { decompress: 'none' })).toBe(StateTransformer.UpdateResult.Recreate);
      expect(update({}, { decompress: 'auto' })).not.toBe(StateTransformer.UpdateResult.Recreate);
    }
  });
});
