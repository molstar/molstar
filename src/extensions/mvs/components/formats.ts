/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { hashFnv32a } from '../../../mol-data/util';
import { DataFormatProvider } from '../../../mol-plugin-state/formats/provider';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { Download } from '../../../mol-plugin-state/transforms/data';
import { PluginContext } from '../../../mol-plugin/context';
import { StateAction, StateObjectRef } from '../../../mol-state';
import { RuntimeContext, Task } from '../../../mol-task';
import { Asset, AssetManager } from '../../../mol-util/assets';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { unzip } from '../../../mol-util/zip/zip';
import { loadMVS } from '../load';
import { MVSData } from '../mvs-data';
import { MVSTransform } from './annotation-structure-component';


/** Plugin state object storing `MVSData` */
export class Mvs extends SO.Create<{ mvsData: MVSData, sourceUrl?: string }>({ name: 'MVS Data', typeClass: 'Data' }) { }

/** Transformer for parsing data in MVSJ format */
export const ParseMVSJ = MVSTransform({
    name: 'mvs-parse-mvsj',
    display: { name: 'MolViewSpec from MVSJ', description: 'Create MolViewSpec view from MVSJ data' },
    from: SO.Data.String,
    to: Mvs,
})({
    apply({ a }, plugin: PluginContext) {
        const mvsData = MVSData.fromMVSJ(a.data);
        const sourceUrl = tryGetDownloadUrl(a, plugin);
        return new Mvs({ mvsData, sourceUrl });
    },
});

/** Transformer for parsing data in MVSX format (= zipped MVSJ + referenced files like structures and annotations) */
export const ParseMVSX = MVSTransform({
    name: 'mvs-parse-mvsx',
    display: { name: 'MolViewSpec from MVSX', description: 'Create MolViewSpec view from MVSX data' },
    from: SO.Data.Binary,
    to: Mvs,
    params: {
        mainFilePath: PD.Text('index.mvsj'),
    },
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Parse MVSX file', async ctx => {
            const data = await loadMVSX(plugin, ctx, a.data, params.mainFilePath);
            return new Mvs(data);
        });
    },
});


/** Params for the `LoadMvsData` action */
export const LoadMvsDataParams = {
    replaceExisting: PD.Boolean(false, { description: 'If true, the loaded MVS view will replace the current state; if false, the MVS view will be added to the current state.' }),
    keepCamera: PD.Boolean(false, { description: 'If true, any camera positioning from the MVS state will be ignored and the current camera position will be kept.' }),
    applyExtensions: PD.Boolean(true, { description: 'If true, apply builtin MVS-loading extensions (not a part of standard MVS specification).' }),
};

/** State action which loads a MVS view into Mol* */
export const LoadMvsData = StateAction.build({
    display: { name: 'Load MVS Data' },
    from: Mvs,
    params: LoadMvsDataParams,
})(({ a, params }, plugin: PluginContext) => Task.create('Load MVS Data', async () => {
    const { mvsData, sourceUrl } = a.data;
    await loadMVS(plugin, mvsData, { replaceExisting: params.replaceExisting, keepCamera: params.keepCamera, sourceUrl: sourceUrl, extensions: params.applyExtensions ? undefined : [] });
}));


/** Data format provider for MVSJ format.
 * If Visuals:On, it will load the parsed MVS view;
 * otherwise it will just create a plugin state object with parsed data. */
export const MVSJFormatProvider: DataFormatProvider<{}, StateObjectRef<Mvs>, any> = DataFormatProvider({
    label: 'MVSJ',
    description: 'MVSJ',
    category: 'Miscellaneous',
    stringExtensions: ['mvsj'],
    parse: async (plugin, data) => {
        return plugin.state.data.build().to(data).apply(ParseMVSJ).commit();
    },
    visuals: async (plugin, data) => {
        const ref = StateObjectRef.resolveRef(data);
        const params = PD.getDefaultValues(LoadMvsDataParams);
        return await plugin.state.data.applyAction(LoadMvsData, params, ref).run();
    },
});

/** Data format provider for MVSX format.
 * If Visuals:On, it will load the parsed MVS view;
 * otherwise it will just create a plugin state object with parsed data. */
export const MVSXFormatProvider: DataFormatProvider<{}, StateObjectRef<Mvs>, any> = DataFormatProvider({
    label: 'MVSX',
    description: 'MVSX',
    category: 'Miscellaneous',
    binaryExtensions: ['mvsx'],
    parse: async (plugin, data) => {
        return plugin.state.data.build().to(data).apply(ParseMVSX).commit();
    },
    visuals: MVSJFormatProvider.visuals,
});


/** Parse binary data `data` as MVSX archive,
 * add all contained files to `plugin`'s asset manager,
 * and parse the main file in the archive as MVSJ.
 * Return parsed MVS data and `sourceUrl` for resolution of relative URIs.  */
export async function loadMVSX(plugin: PluginContext, runtimeCtx: RuntimeContext, data: Uint8Array, mainFilePath: string = 'index.mvsj'): Promise<{ mvsData: MVSData, sourceUrl: string }> {
    const archiveId = `ni,fnv1a;${hashFnv32a(data)}`;
    let files: { [path: string]: Uint8Array };
    try {
        files = await unzip(runtimeCtx, data) as typeof files;
    } catch (err) {
        plugin.log.error('Invalid MVSX file');
        throw err;
    }
    for (const path in files) {
        const url = arcpUri(archiveId, path);
        ensureUrlAsset(plugin.managers.asset, url, files[path]);
    }
    const mainFile = files[mainFilePath];
    if (!mainFile) throw new Error(`File ${mainFilePath} not found in the MVSX archive`);
    const mvsData = MVSData.fromMVSJ(decodeUtf8(mainFile));
    const sourceUrl = arcpUri(archiveId, mainFilePath);
    return { mvsData, sourceUrl };
}

/** If the PluginStateObject `pso` comes from a Download transform, try to get its `url` parameter. */
function tryGetDownloadUrl(pso: SO.Data.String, plugin: PluginContext): string | undefined {
    const theCell = plugin.state.data.selectQ(q => q.ofTransformer(Download)).find(cell => cell.obj === pso);
    const urlParam = theCell?.transform.params?.url;
    return urlParam ? Asset.getUrl(urlParam) : undefined;
}

/** Return a URI referencing a file within an archive, using ARCP scheme (https://arxiv.org/pdf/1809.06935.pdf).
 * `archiveId` corresponds to the `authority` part of URI (e.g. 'uuid,EYVwjDiZhM20PWbF1OWWvQ' or 'ni,fnv1a;938511930')
 * `path` corresponds to the path to a file within the archive */
function arcpUri(archiveId: string, path: string): string {
    return new URL(path, `arcp://${archiveId}/`).href;
}

/** Add a URL asset to asset manager.
 * Skip if an asset with the same URL already exists. */
function ensureUrlAsset(manager: AssetManager, url: string, data: Uint8Array) {
    const asset = Asset.getUrlAsset(manager, url);
    if (!manager.has(asset)) {
        const filename = url.split('/').pop() ?? 'file';
        manager.set(asset, new File([data], filename));
    }
}

/** Decode bytes to text using UTF-8 encoding */
function decodeUtf8(bytes: Uint8Array): string {
    _decoder ??= new TextDecoder();
    return _decoder.decode(bytes);
}
let _decoder: TextDecoder | undefined;
