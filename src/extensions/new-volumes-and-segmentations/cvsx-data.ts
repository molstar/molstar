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
                    update.to(visual).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.alpha = 0.5; });
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

    // TODO: need to first parse annotations.json and pass this information
    // to processCvsxFile for parsing segmentation


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

// export interface CVSXData {
//     /** MolViewSpec tree */
//     // root: MVSTree,
//     /** Associated metadata */
//     // metadata: MVSMetadata,
//     volume: 
//     segmentation:
// }

// interface MVSMetadata {
//     /** Version of the spec used to write this tree */
//     version: string,
//     /** Name of this view */
//     title?: string,
//     /** Detailed description of this view */
//     description?: string,
//     /** Format of the description */
//     description_format?: 'markdown' | 'plaintext',
//     /** Timestamp when this view was exported */
//     timestamp: string,
// }

// export const CVSXData = {
//     /** Currently supported major version of MolViewSpec format (e.g. 1 for version '1.0.8') */
//     SupportedVersion: 1,

//     /** Parse MVSJ (MolViewSpec-JSON) format to `MVSData`. Does not include any validation. */
//     fromCVSX(data: {[k: string]: Uint8Array}): CVSXData {
//         // const result: MVSData = JSON.parse(mvsjString);
//         const result: CVSXData = Unzip(data)
//         const major = majorVersion(result?.metadata?.version);
//         if (major === undefined) {
//             console.error('Loaded MVS does not contain valid version info.');
//         } else if (major > (majorVersion(MVSData.SupportedVersion) ?? 0)) {
//             console.warn(`Loaded MVS is of higher version (${result.metadata.version}) than currently supported version (${MVSData.SupportedVersion}). Some features may not work as expected.`);
//         }
//         return result;
//     },

//     /** Encode `MVSData` to MVSJ (MolViewSpec-JSON) string. Use `space` parameter to control formatting (as with `JSON.stringify`). */
//     toMVSJ(mvsData: MVSData, space?: string | number): string {
//         return JSON.stringify(mvsData, undefined, space);
//     },

//     /** Validate `MVSData`. Return `true` if OK; `false` if not OK.
//      * If `options.noExtra` is true, presence of any extra node parameters is treated as an issue. */
//     isValid(mvsData: MVSData, options: { noExtra?: boolean } = {}): boolean {
//         return MVSData.validationIssues(mvsData, options) === undefined;
//     },

//     /** Validate `MVSData`. Return `undefined` if OK; list of issues if not OK.
//      * If `options.noExtra` is true, presence of any extra node parameters is treated as an issue. */
//     validationIssues(mvsData: MVSData, options: { noExtra?: boolean } = {}): string[] | undefined {
//         const version = mvsData?.metadata?.version;
//         if (typeof version !== 'string') return [`"version" in MVS must be a string, not ${typeof version}: ${version}`];
//         if (mvsData.root === undefined) return [`"root" missing in MVS`];
//         return treeValidationIssues(MVSTreeSchema, mvsData.root, options);
//     },

//     /** Return a human-friendly textual representation of `mvsData`. */
//     toPrettyString(mvsData: MVSData): string {
//         const title = mvsData.metadata.title !== undefined ? ` "${mvsData.metadata.title}"` : '';
//         return `MolViewSpec tree${title} (version ${mvsData.metadata.version}, created ${mvsData.metadata.timestamp}):\n${treeToString(mvsData.root)}`;
//     },

//     /** Create a new MolViewSpec builder containing only a root node. Example of MVS builder usage:
//      *
//      * ```
//      * const builder = MVSData.createBuilder();
//      * builder.canvas({ background_color: 'white' });
//      * const struct = builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1og2_updated.cif' }).parse({ format: 'mmcif' }).modelStructure();
//      * struct.component().representation().color({ color: '#3050F8' });
//      * console.log(MVSData.toPrettyString(builder.getState()));
//      * ```
//      */
//     createBuilder(): Root {
//         return createMVSBuilder();
//     },
// };


// /** Get the major version from a semantic version string, e.g. '1.0.8' -> 1 */
// function majorVersion(semanticVersion: string | number): number | undefined {
//     if (typeof semanticVersion === 'string') return parseInt(semanticVersion.split('.')[0]);
//     if (typeof semanticVersion === 'number') return Math.floor(semanticVersion);
//     console.error(`Version should be a string, not ${typeof semanticVersion}: ${semanticVersion}`);
//     return undefined;
// }
