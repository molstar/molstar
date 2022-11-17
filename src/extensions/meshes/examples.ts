/** Testing examples for using mesh-extension.ts. */

import * as MS from './molstar-lib-imports';

import { ParseMeshlistTransformer, MeshShapeTransformer, MeshlistData } from './mesh-extension';
import * as MeshUtils from './mesh-utils';
import { BACKGROUND_OPACITY, FOREROUND_OPACITY, InitMeshStreaming } from './mesh-streaming/transformers';
import { MeshServerInfo } from './mesh-streaming/server-info';


export const DB_URL = '/db'; // local
// DB_URL = 'http://sestra.ncbr.muni.cz/data/cellstar-sample-data/db'; // download


export async function runMeshExtensionExamples(plugin: MS.PluginUIContext, db_url: string = DB_URL) {
    console.time('TIME MESH EXAMPLES');
    // await runIsosurfaceExample(plugin, db_url);
    // await runMolsurfaceExample(plugin);

    // Focused Ion Beam-Scanning Electron Microscopy of mitochondrial reticulum in murine skeletal muscle: https://www.ebi.ac.uk/empiar/EMPIAR-10070/
    // await runMeshExample(plugin, 'all', db_url);
    // await runMeshExample(plugin, 'fg', db_url);
    // await runMultimeshExample(plugin, 'fg', 'worst', db_url);
    // await runCifMeshExample(plugin);
    // await runMeshExample2(plugin, 'fg');
    await runMeshStreamingExample(plugin);

    console.timeEnd('TIME MESH EXAMPLES');
}

/** Example for downloading multiple separate segments, each containing 1 mesh. */
export async function runMeshExample(plugin: MS.PluginUIContext, segments: 'fg' | 'all', db_url: string = DB_URL) {
    const detail = 2;
    const segmentIds = (segments === 'all') ?
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17] // segment-16 has no detail-2
        : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 17]; // segment-13 and segment-15 are quasi background

    for (const segmentId of segmentIds) {
        await createMeshFromUrl(plugin, `${db_url}/empiar-10070-mesh-rounded/segment-${segmentId}/detail-${detail}`, segmentId, detail, true, true, undefined);
    }
}

/** Example for downloading multiple separate segments, each containing 1 mesh. */
export async function runMeshExample2(plugin: MS.PluginUIContext, segments: 'one' | 'few' | 'fg' | 'all') {
    const detail = 1;
    const segmentIds = (segments === 'one') ? [15]
        : (segments === 'few') ? [1, 4, 7, 10, 16]
            : (segments === 'all') ? [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17] // segment-16 has no detail-2
                : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 17]; // segment-13 and segment-15 are quasi background

    for (const segmentId of segmentIds) {
        await createMeshFromUrl(plugin, `http://localhost:9000/v2/empiar/empiar-10070/mesh_bcif/${segmentId}/${detail}`, segmentId, detail, false, true, undefined);
    }
}

/** Example for downloading a single segment containing multiple meshes. */
export async function runMultimeshExample(plugin: MS.PluginUIContext, segments: 'fg' | 'all', detailChoice: 'best' | 'worst', db_url: string = DB_URL) {
    const urlDetail = (detailChoice === 'best') ? '2' : 'worst';
    const numDetail = (detailChoice === 'best') ? 2 : 1000;
    await createMeshFromUrl(plugin, `${db_url}/empiar-10070-multimesh-rounded/segments-${segments}/detail-${urlDetail}`, 0, numDetail, false, true, undefined);
}

/** Download data and create state tree hierarchy down to visual representation. */
export async function createMeshFromUrl(plugin: MS.PluginContext, meshDataUrl: string, segmentId: number, detail: number, 
        collapseTree: boolean, log: boolean, color?: MS.Color, parent?: MS.StateObjectSelector, transparentIfBboxAbove?: number) {

    // PARAMS - Depend on the type of transformer T -> Params<T>
    // 1st argument to plugin.builders.data.rawData, 2nd argument to .apply

    // OPTIONS - Same for each type of transformer
    // Last argument to plugin.builders.data.rawData, plugin.builders.data.download, .apply
    // interface StateTransform.Options {
    //     ref?: string,  // State tree node ID (default: auto-generated ID)
    //     tags?: string | string[],  // I don't know what this is for
    //     state?: {
    //         isGhost?: boolean,  // is the cell shown in the UI
    //         isLocked?: boolean,  // can the corresponding be deleted by the user.
    //         isHidden?: boolean,  // is the representation associated with the cell hidden
    //         isCollapsed?: boolean,  // is the tree node collapsed?
    //     },
    //     dependsOn?: string[]  // references to other nodes, I think
    // }

    const update = parent ? plugin.build().to(parent) : plugin.build().toRoot();

    // RAW DATA NODE
    const rawDataNode = await update.apply(MS.Download,
        { url: meshDataUrl, isBinary: true, label: `Downloaded Data ${segmentId}` }, // params
        { tags: ['What', 'are', 'tags', 'good', 'for?'], state: { isCollapsed: collapseTree } } // options
    ).commit();
    if (log) console.log('rawDataNode:', rawDataNode);

    const cifNode = await plugin.build().to(rawDataNode).apply(MS.StateTransforms.Data.ParseCif).commit();
    if (log) console.log('cifNode:', rawDataNode);

    // PARSED DATA NODE
    const parsedDataNode = await plugin.build().to(cifNode).apply(
        ParseMeshlistTransformer,
        { label: undefined, segmentId: segmentId, segmentName: `Segment ${segmentId}`, detail: detail }, // params
        { ref: `ref-parsed-data-node-${segmentId}` } // options
    ).commit();
    if (log) console.log('parsedDataNode:', parsedDataNode);
    if (log) console.log('parsedDataNode.data:', parsedDataNode.data);
    if (log) console.log('parsedDataNode mesh list stats:\n', MeshlistData.stats(parsedDataNode.data!));

    // MESH SHAPE NODE
    const shapeNode = await plugin.build().to(parsedDataNode).apply(MeshShapeTransformer,
        { color: color }, // options
        { ref: `ref-shape-node-${segmentId}` }
    ).commit();
    if (log) console.log('shapeNode:', shapeNode);
    if (log) console.log('shapeNode.data:', shapeNode.data);

    
    let transparent = false;
    if (transparentIfBboxAbove !== undefined && parsedDataNode.data) {
        const bbox = MeshlistData.bbox(parsedDataNode.data) || MS.Box3D.zero();
        transparent = MS.Box3D.volume(bbox) > transparentIfBboxAbove;
    }

    // MESH REPR NODE
    const reprNode = await plugin.build().to(shapeNode).apply(MS.ShapeRepresentation3D,
        { alpha: transparent ? BACKGROUND_OPACITY : FOREROUND_OPACITY },
        { ref: `ref-repr-nod+e${segmentId}` }
    ).commit();
    if (log) console.log('reprNode:', reprNode);
    if (log) console.log('reprNode.data:', reprNode.data);

    return rawDataNode;

}

export async function runMeshStreamingExample(plugin: MS.PluginUIContext, source: MeshServerInfo.MeshSource = 'empiar', entryId: string = 'empiar-10070', serverUrl?: string, parent?: MS.StateObjectSelector) {
    const params = MS.ParamDefinition.getDefaultValues(MeshServerInfo.Params);
    if (serverUrl) params.serverUrl = serverUrl;
    params.source = source;
    params.entryId = entryId;
    await plugin.runTask(plugin.state.data.applyAction(InitMeshStreaming, params, parent?.ref), { useOverlay: false });
}

/** Example for downloading a protein structure and visualizing molecular surface. */
export async function runMolsurfaceExample(plugin: MS.PluginUIContext) {
    const entryId = 'pdb-7etq';

    // Node "https://www.ebi.ac.uk/pdbe/entry-files/download/7etq.bcif" ("transformer": "ms-plugin.download") -> var data
    const data = await plugin.builders.data.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/7etq.bcif', isBinary: true }, { state: { isGhost: false } });
    console.log('formats:', plugin.dataFormats.list);

    // Node "CIF File" ("transformer": "ms-plugin.parse-cif")
    // Node "7ETQ 1 model" ("transformer": "ms-plugin.trajectory-from-mmcif") -> var trajectory
    const parsed = await plugin.dataFormats.get('mmcif')!.parse(plugin, data, { entryId });
    const trajectory: MS.StateObjectSelector<MS.PluginStateObject.Molecule.Trajectory> = parsed.trajectory;
    console.log('parsed', parsed);
    console.log('trajectory', trajectory);

    // Node "Model 1" ("transformer": "ms-plugin.model-from-trajectory") -> var model
    const model = await plugin.build().to(trajectory).apply(MS.StateTransforms.Model.ModelFromTrajectory).commit();
    console.log('model:', model);

    // Node "Model 91 elements" ("transformer": "ms-plugin.structure-from-model") -> var structure
    const structure = await plugin.build().to(model).apply(MS.StateTransforms.Model.StructureFromModel,).commit();
    console.log('structure:', structure);

    // Node "Molecular Surface" ("transformer": "ms-plugin.structure-representation-3d") -> var repr
    const reprParams = MS.createStructureRepresentationParams(plugin, undefined, { type: 'molecular-surface' });
    const repr = await plugin.build().to(structure).apply(MS.StateTransforms.Representation.StructureRepresentation3D, reprParams).commit();
    console.log('repr:', repr);
}

/** Example for downloading an EMDB density data and visualizing isosurface. */
export async function runIsosurfaceExample(plugin: MS.PluginUIContext, db_url: string = DB_URL) {
    const entryId = 'emd-1832';
    const isoLevel = 2.73;

    let root = await plugin.build();
    const data = await plugin.builders.data.download({ url: `${db_url}/emd-1832-box`, isBinary: true }, { state: { isGhost: false } });
    const parsed = await plugin.dataFormats.get('dscif')!.parse(plugin, data, { entryId });

    const volume: MS.StateObjectSelector<MS.PluginStateObject.Volume.Data> = parsed.volumes?.[0] ?? parsed.volume;
    const volumeData = volume.cell!.obj!.data;
    console.log('data:', data);
    console.log('parsed:', parsed);
    console.log('volume:', volume);
    console.log('volumeData:', volumeData);

    root = await plugin.build();
    console.log('root:', root);
    console.log('to:', root.to(volume));
    console.log('toRoot:', root.toRoot());

    let volumeParams;
    volumeParams = MS.createVolumeRepresentationParams(plugin, volumeData, {
        type: 'isosurface',
        typeParams: {
            alpha: 0.5,
            isoValue: MS.Volume.adjustedIsoValue(volumeData, isoLevel, 'relative'),
            visuals: ['solid'],
            sizeFactor: 1,
        },
        color: 'uniform',
        colorParams: { value: MS.Color(0x00aaaa) },

    });
    root.to(volume).apply(MS.StateTransforms.Representation.VolumeRepresentation3D, volumeParams);

    volumeParams = MS.createVolumeRepresentationParams(plugin, volumeData, {
        type: 'isosurface',
        typeParams: {
            alpha: 1.0,
            isoValue: MS.Volume.adjustedIsoValue(volumeData, isoLevel, 'relative'),
            visuals: ['wireframe'],
            sizeFactor: 1,
        },
        color: 'uniform',
        colorParams: { value: MS.Color(0x8800aa) },

    });
    root.to(volume).apply(MS.StateTransforms.Representation.VolumeRepresentation3D, volumeParams);
    await root.commit();
}


export async function runCifMeshExample(plugin: MS.PluginUIContext, api: string = 'http://localhost:9000/v2',
    source: MeshServerInfo.MeshSource = 'empiar', entryId: string = 'empiar-10070',
    segmentId: number = 1, detail: number = 10,
) {
    const url = `${api}/${source}/${entryId}/mesh_bcif/${segmentId}/${detail}`;
    getMeshFromBcif(plugin, url);
}

async function getMeshFromBcif(plugin: MS.PluginUIContext, url: string) {
    const urlAsset = MS.Asset.getUrlAsset(plugin.managers.asset, url); // QUESTION how is urlAsset better than normal `fetch`
    const asset = await plugin.runTask(plugin.managers.asset.resolve(urlAsset, 'binary'));
    const parsed = await plugin.runTask(MS.CIF.parseBinary(asset.data));
    if (parsed.isError) {
        plugin.log.error('VolumeStreaming, parsing CIF: ' + parsed.toString());
        return;
    }
    console.log('blocks:', parsed.result.blocks);
    const mesh = await MeshUtils.meshFromCif(parsed.result);
    console.log(mesh);
}