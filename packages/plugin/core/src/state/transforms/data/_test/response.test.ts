/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
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
import { DeflateData, Download, DownloadBlob, ReadFile, ReadResponse } from '../fetch.js';
import { ParseBlob } from '../../../formats/cif.js';
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

  it('reads as a response by default; asResponse: false gives String Data as before', async () => {
    const plugin = await createPlugin();
    const byDefault = (await plugin.builders.data.readFile({ file: cifFile() })).data;
    expect(byDefault.cell!.obj).toBeInstanceOf(SO.Data.Response);
    const { data } = await plugin.builders.data.readFile({ file: cifFile(), asResponse: false });
    expect(data.cell!.obj).toBeInstanceOf(SO.Data.String);
    plugin.dispose();
  });

  it('every parser reads String, Binary and Response Data alike', async () => {
    const plugin = await createPlugin();
    const pdb = 'ATOM      1  CA  ALA A   1       1.500   2.500   3.500  1.00  0.00           C\nEND\n';
    const file = () => Asset.File(new File([pdb], 'x.pdb'));
    const counts: number[] = [];
    for (const asResponse of [true, false]) {
      for (const isBinary of [false, true]) {
        const { data } = await plugin.builders.data.readFile({ file: file(), asResponse, isBinary });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        counts.push(trajectory.cell!.obj!.data.representative.atomicHierarchy.atoms._rowCount);
      }
    }
    expect(counts).toEqual([1, 1, 1, 1]);
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
      expect(update({}, { asResponse: false })).toBe(StateTransformer.UpdateResult.Recreate);
      expect(update({}, { asResponse: true })).not.toBe(StateTransformer.UpdateResult.Recreate);
      expect(update({ isBinary: false }, { isBinary: true })).toBe(StateTransformer.UpdateResult.Recreate);
    }
  });

  it('OpenFiles reads every format as a response', async () => {
    const plugin = await createPlugin();

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
    expect(types.get('b.pdb')).toBe('Response Data');
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
      const bytes = enc.encode(t);
      expect(guessCifVariant(info, bytes)).toBe(guessCifVariant(info, t));
      // asked again for the same bytes (as by every CIF provider in turn): the cached answer
      expect(guessCifVariant(info, bytes)).toBe(guessCifVariant(info, t));
      expect(looksLikeMmcifParticles(info, enc.encode(t))).toBe(looksLikeMmcifParticles(info, t));
    }
  });
});

describe('decompress (util/_spec/assets.md §5)', () => {
  const gz = (text: string) => new Uint8Array(gzipSync(new TextEncoder().encode(text)));

  it('ReadFile inflates gzip by its bytes unless decompress is none', async () => {
    const plugin = await createPlugin();
    const file = () => Asset.File(new File([gz(Cif)], 'x.cif'));
    const text = (await plugin.builders.data.readFile({ file: file(), asResponse: false })).data;
    expect(`${text.cell!.obj!.data}`).toBe(Cif);
    const response = (await plugin.builders.data.readFile({ file: file(), asResponse: true })).data;
    expect(`${await response.cell!.obj!.data.text(RuntimeContext.Synchronous)}`).toBe(Cif);
    const raw = (
      await plugin.builders.data.readFile({ file: file(), isBinary: true, asResponse: false, decompress: 'none' })
    ).data;
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

  it('ReadFile recreates when its file changes', () => {
    const update = (oldParams: any, newParams: any) =>
      ReadFile.definition.update!(
        { oldParams, newParams, a: undefined as any, b: undefined as any, cache: {} } as any,
        undefined as any,
      );
    const a = Asset.File(new File(['a'], 'a.cif'));
    expect(update({ file: a }, { file: Asset.File(new File(['b'], 'b.cif')) })).toBe(
      StateTransformer.UpdateResult.Recreate,
    );
    expect(update({ file: a }, { file: a })).toBe(StateTransformer.UpdateResult.Unchanged);
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

describe('Data.ResponseBlob', () => {
  const gz = (text: string) => new Uint8Array(gzipSync(new TextEncoder().encode(text)));

  it('DownloadBlob with asResponse holds each entry as bytes (inflated), and ParseBlob parses them', async () => {
    const plugin = await createPlugin();
    const urls = ['https://x.org/a.cif', 'https://x.org/b.cif.gz'];
    // held assets: the downloads read them, nothing reaches the network
    plugin.managers.asset.set(Asset.Url(urls[0]), new File([Cif], 'raw-data'), { isStatic: true });
    plugin.managers.asset.set(Asset.Url(urls[1]), new File([gz(Cif)], 'raw-data'), { isStatic: true });
    const sources = urls.map((url, i) => ({ id: `${i}`, url }));
    const blob = await plugin.builders.data.downloadBlob({ sources });
    const obj = blob.cell!.obj!;
    expect(obj).toBeInstanceOf(SO.Data.ResponseBlob);
    expect(obj.data.map((e) => e.id)).toEqual(['0', '1']);
    for (const e of obj.data) expect(`${await e.data.text(RuntimeContext.Synchronous)}`).toBe(Cif);

    const formats = sources.map((s) => ({ id: s.id, format: 'cif' as const }));
    const parsed = await plugin.build().to(blob).apply(ParseBlob, { formats }).commit();
    expect(parsed.cell!.obj!.data.map((e) => e.data.blocks[0].categories.atom_site.rowCount)).toEqual([2, 2]);
    plugin.dispose();
  });

  it('DownloadBlob with asResponse: false gives Data.Blob as before, gzip inflated', async () => {
    const plugin = await createPlugin();
    const url = 'https://x.org/c.cif.gz';
    plugin.managers.asset.set(Asset.Url(url), new File([gz(Cif)], 'raw-data'), { isStatic: true });
    const blob = await plugin
      .build()
      .toRoot()
      .apply(DownloadBlob, { sources: [{ id: '0', url }], asResponse: false })
      .commit();
    const e = blob.cell!.obj!.data[0];
    expect(blob.cell!.obj).toBeInstanceOf(SO.Data.Blob);
    expect(e.kind).toBe('string');
    expect(`${e.data}`).toBe(Cif);
    plugin.dispose();
  });
});

describe('sessions (.molx) through the asset manager', () => {
  const gz = (text: string) => new Uint8Array(gzipSync(new TextEncoder().encode(text)));

  /** The trajectories of a plugin: their source cell's type and atom count */
  function trajectories(plugin: PluginContext) {
    const out: { data: string; atoms: number }[] = [];
    plugin.state.data.cells.forEach((cell) => {
      if (!SO.Molecule.Trajectory.is(cell.obj)) return;
      let ref = cell.transform.parent;
      while (ref && !(plugin.state.data.cells.get(ref)?.obj instanceof SO.Data.Response)) {
        ref = plugin.state.data.cells.get(ref)!.transform.parent;
      }
      out.push({
        data: plugin.state.data.cells.get(ref)?.obj?.type.name ?? '',
        atoms: cell.obj.data.representative.atomicHierarchy.atoms._rowCount,
      });
    });
    return out;
  }

  async function roundTrip(load: (plugin: PluginContext) => Promise<void>) {
    const plugin = await createPlugin();
    await load(plugin);
    const before = trajectories(plugin);
    const session = await plugin.managers.snapshot.serialize({ type: 'zip', params: { camera: false } });
    plugin.dispose();

    const restored = await createPlugin();
    await restored.managers.snapshot.open(new File([session], 'session.molx'));
    const after = trajectories(restored);
    restored.dispose();
    return { before, after };
  }

  it('an opened file, a download and a gzipped download come back from a session, without the network', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch');
    const { before, after } = await roundTrip(async (plugin) => {
      const file = Asset.File(new File([Cif], 'a.cif'));
      const read = (await plugin.builders.data.readFile({ file })).data;
      await plugin.builders.structure.parseTrajectory(read, 'mmcif');
      for (const [url, bytes] of [
        ['https://x.org/b.cif', new TextEncoder().encode(Cif)],
        ['https://x.org/c.cif.gz', gz(Cif)],
      ] as const) {
        // held as a download would be: the session stores what the manager holds
        const asset = Asset.Url(url);
        plugin.managers.asset.set(asset, new File([bytes], 'raw-data'));
        const data = await plugin.builders.data.download({ url: asset });
        await plugin.builders.structure.parseTrajectory(data, 'mmcif');
      }
    });
    expect(before).toEqual([
      { data: 'Response Data', atoms: 2 },
      { data: 'Response Data', atoms: 2 },
      { data: 'Response Data', atoms: 2 },
    ]);
    expect(after).toEqual(before);
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it('a session whose Download predates asResponse loads it as a response', async () => {
    const { after } = await roundTrip(async (plugin) => {
      const asset = Asset.Url('https://x.org/old.cif');
      plugin.managers.asset.set(asset, new File([Cif], 'raw-data'));
      // the params a session saved before this change: no asResponse, no decompress
      const data = await plugin.build().toRoot().apply(Download, { url: asset, isBinary: false }).commit();
      expect(data.cell!.obj).toBeInstanceOf(SO.Data.Response);
      await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    });
    expect(after).toEqual([{ data: 'Response Data', atoms: 2 }]);
  });

  it('a real download (held by the URL cache, gzip without Content-Encoding) comes back from a session', async () => {
    const fetchSpy = vi.spyOn(globalThis, 'fetch').mockImplementation(async () => new Response(gz(Cif)));
    const plugin = await createPlugin();
    const data = await plugin.builders.data.download({ url: 'https://x.org/d.cif.gz' });
    await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    expect(fetchSpy).toHaveBeenCalledTimes(1);
    const before = trajectories(plugin);
    const session = await plugin.managers.snapshot.serialize({ type: 'zip', params: { camera: false } });
    plugin.dispose();

    fetchSpy.mockClear();
    const restored = await createPlugin();
    await restored.managers.snapshot.open(new File([session], 'session.molx'));
    expect(trajectories(restored)).toEqual(before);
    expect(before).toEqual([{ data: 'Response Data', atoms: 2 }]);
    expect(fetchSpy).not.toHaveBeenCalled();
    restored.dispose();
    fetchSpy.mockRestore();
  });
});

describe('one form at a time (ResponseBody)', () => {
  it('turns to text once, lets go of its bytes, and keeps its header and sniffs', () => {
    let calls = 0;
    const body = new SO.Data.ResponseBody(new TextEncoder().encode(Cif) as Uint8Array<ArrayBuffer>, {
      onText: () => calls++,
    });
    expect(body.form).toBe('bytes');
    expect(body.sniffText('_atom_site.Cartn_x')).toBe(true);
    expect(body.form).toBe('bytes');
    const text = body.textNow();
    expect(`${text}`).toBe(Cif);
    expect(body.form).toBe('text');
    expect(body.heldBytes).toBeUndefined();
    expect(body.textNow()).toBe(text);
    expect(calls).toBe(1);
    expect(body.sniffText('_atom_site.Cartn_x')).toBe(true);
    expect(body.sniffBytes('data_test')).toBe(true);
    expect(body.form).toBe('text');
    // bytes again: encoded from the text, and held instead of it
    expect(new TextDecoder().decode(body.bytes)).toBe(Cif);
    expect(body.form).toBe('bytes');
  });

  it('a text format opened as a file holds one copy: the text, the file keeps the rest; CIF stays on bytes', async () => {
    const plugin = await createPlugin();
    const pdb = 'ATOM      1  CA  ALA A   1       1.500   2.500   3.500  1.00  0.00           C\nEND\n';
    const pdbFile = Asset.File(new File([pdb], 'x.pdb'));
    const cifFile = Asset.File(new File([Cif], 'y.cif'));
    await plugin.runTask(
      plugin.state.data.applyAction(OpenFiles, {
        files: [pdbFile, cifFile],
        format: { name: 'auto', params: {} },
        visuals: false,
      }),
    );
    const bodies = new Map<string, SO.Data.ResponseBody>();
    plugin.state.data.cells.forEach((cell) => {
      if (cell.obj instanceof SO.Data.Response) bodies.set(cell.obj.label, cell.obj.data);
    });
    expect(bodies.get('x.pdb')!.form).toBe('text');
    expect(bodies.get('y.cif')!.form).toBe('bytes');
    // the asset manager let go of the PDB's bytes: held as text by its body, as a file on disk by the manager
    expect((plugin.managers.asset.get(pdbFile) as any).bytes).toBeUndefined();
    expect((plugin.managers.asset.get(cifFile) as any).bytes).toBe(bodies.get('y.cif')!.heldBytes);
    plugin.dispose();
  });

  it('a session keeps a download whose body turned to text', async () => {
    const fetchSpy = vi
      .spyOn(globalThis, 'fetch')
      .mockImplementation(
        async () =>
          new Response('ATOM      1  CA  ALA A   1       1.500   2.500   3.500  1.00  0.00           C\nEND\n'),
      );
    const plugin = await createPlugin();
    const data = await plugin.builders.data.download({ url: 'https://x.org/e.pdb' });
    await plugin.builders.structure.parseTrajectory(data, 'pdb');
    expect(data.cell!.obj!.data.form).toBe('text');
    const session = await plugin.managers.snapshot.serialize({ type: 'zip', params: { camera: false } });
    plugin.dispose();
    fetchSpy.mockClear();
    const restored = await createPlugin();
    await restored.managers.snapshot.open(new File([session], 'session.molx'));
    let atoms = 0;
    restored.state.data.cells.forEach((cell) => {
      if (SO.Molecule.Trajectory.is(cell.obj)) atoms += cell.obj.data.representative.atomicHierarchy.atoms._rowCount;
    });
    expect(atoms).toBe(1);
    expect(fetchSpy).not.toHaveBeenCalled();
    restored.dispose();
    fetchSpy.mockRestore();
  });
});
