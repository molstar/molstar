/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 */

import { Volume } from "../../mol-model/volume";
import { createVolumeRepresentationParams } from "../../mol-plugin-state/helpers/volume-representation-params";
import { PluginStateObject } from "../../mol-plugin-state/objects";
import { StateTransforms } from "../../mol-plugin-state/transforms";
import { VolumeRepresentation3D } from "../../mol-plugin-state/transforms/representation";
import { PluginContext } from "../../mol-plugin/context";
import { Asset } from "../../mol-util/assets";
import { Color } from "../../mol-util/color";
import { getFileNameInfo } from "../../mol-util/file-info";
import { Unzip } from "../../mol-util/zip/zip";
import { AnnotationMetadata } from "./volseg-api/data";


const DEFAULT_SEGMENT_COLOR = Color.fromNormalizedRgb(0.8, 0.8, 0.8);

// export function VolsegControls({ plugin, data, setData }: { plugin: PluginContext, data: VolsegUIData, setData: (d: VolsegUIData) => void }) {
//     const entryData = data.activeNode?.data;
//     if (!entryData) {
//         return <p>No data!</p>;
//     }
//     if (!data.globalState) {
//         return <p>No global state!</p>;
//     }

//     const params = {
//         /** Reference to the active VolsegEntry node */
//         entry: PD.Select(data.activeNode!.data.ref, data.availableNodes.map(entry => [entry.data.ref, entry.data.entryId]))
//     };
//     const values: PD.Values<typeof params> = {
//         entry: data.activeNode!.data.ref,
//     };

//     const globalState = useBehavior(data.globalState.currentState);

//     return <>
//         <ParameterControls params={params} values={values} onChangeValues={next => setData(VolsegUIData.changeActiveNode(data, next.entry))} />

//         <TimeFrameSlider entryData={entryData} />
//         <ExpandGroup header='Global options'>
//             <WaitingParameterControls params={VolsegGlobalStateParams} values={globalState} onChangeValues={async next => await data.globalState?.updateState(plugin, next)} />
//         </ExpandGroup>

//         <VolsegEntryControls entryData={entryData} key={entryData.ref} />
//     </>;
// }

export interface OutputByFormat {
    format: string,
    visuals: any[],
    parsed: any[]
}

function createPaletteCVSX(annotations: AnnotationMetadata) {
    const segmentIds = annotations.annotations.map(a => a.segment_id);
    const colorMap = new Map<number, Color>();
    // for (const segment of this.entryData.metadata.allSegments) {
    if (annotations.annotations) {
        for (const annotation of annotations.annotations) {
            if (annotation.color) {
                const color = Color.fromNormalizedArray(annotation.color, 0);
                colorMap.set(annotation.segment_id, color);
            }
        }
        if (colorMap.size === 0) return undefined;
        for (const segid of segmentIds) {
            colorMap.get(segid);
        }
        const colors = segmentIds.map(segid => colorMap.get(segid) ?? DEFAULT_SEGMENT_COLOR);
        return { name: 'colors' as const, params: { list: { kind: 'set' as const, colors: colors } } };
    }
}

export async function updateVisualsBasedOnAnnotations(annotations: AnnotationMetadata, plugin: PluginContext, outputByFormats: OutputByFormat[]) {
    for (const vs of outputByFormats) {
        const { format, visuals, parsed } = vs;
        debugger;
        // TODO: update visual
        if (format === 'dscif') {
            for (const visual of visuals) {
                const update = plugin.build().to(visual.cell.transform.parent);
                for (const visual of visuals) {
                    update.to(visual).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.alpha = 0.2; p.colorTheme.params.value = Color.fromHexStyle('#000000')});
                    // TODO: update volume channel annotations if any similar to labels
                    // TODO: update isolevel?
                    // const isoLevelPromise = ExternalAPIs.tryGetIsovalue(this.entryData.metadata.raw.annotation?.entry_id.source_db_id ?? this.entryData.entryId);
                    // const color = this.entryData.metadata.getVolumeChannelColor(channelId);
                    // const volumeData = volumeNode.cell!.obj!.data;
                    // let label = this.entryData.metadata.getVolumeChannelLabel(channelId);
                    // if (!label) label = channelId.toString();
                }
                await update.commit();
            }
            } else if (format === 'segcif') {
            for (let i = 0; i < visuals.length; i++) {
                const visual = visuals[i];
                const parsedOne: Volume = parsed.volumes[0].data;
                const update = plugin.build().to(visual.cell.transform.parent);
                const params = createVolumeRepresentationParams(plugin, parsedOne, {
                    type: 'segment',
                    typeParams: { tryUseGpu: false },
                    color: 'volume-segment',
                    colorParams: { palette: createPaletteCVSX(annotations) 
                });
                update.to(visual).update(StateTransforms.Representation.VolumeRepresentation3D, p => 
                    params
                    );
                    await update.commit();
            }

        }
    }
}

export async function processCvsxAnnotationsFile(file: Asset.File, plugin: PluginContext) {
    // Parse to interface
    // file.file
    console.log(file);
    const asset = plugin.managers.asset.resolve(file, 'string');
    const data = (await asset.run()).data;
    const parsedData: AnnotationMetadata = JSON.parse(data);
    console.log(parsedData);
    return parsedData;
}

export async function processCvsxFile(file: Asset.File, plugin: PluginContext, format: string, visuals: boolean) {
    // Need to select provider here
    const info = getFileNameInfo(file.file?.name ?? '');
    const isBinary = plugin.dataFormats.binaryExtensions.has(info.ext);
    const { data } = await plugin.builders.data.readFile({ file, isBinary });
    const provider = format === 'auto'
        ? plugin.dataFormats.auto(info, data.cell?.obj!)
        : plugin.dataFormats.get(format);

    if (!provider) {
        plugin.log.warn(`OpenFiles: could not find data provider for '${info.ext}'`);
        await plugin.state.data.build().delete(data).commit();
        return;
    }

    // need to await so that the enclosing Task finishes after the update is done.
    const parsed = await provider.parse(plugin, data);
    if (visuals) {
        const visuals = await provider.visuals?.(plugin, parsed);
        const visualsByFormats: OutputByFormat = {
            format: format,
            visuals: visuals,
            parsed: parsed
        }
        return visualsByFormats;
    }
};


