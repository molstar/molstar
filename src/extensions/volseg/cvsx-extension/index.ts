/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 */

import { Volume } from "../../../mol-model/volume";
import { createVolumeRepresentationParams } from "../../../mol-plugin-state/helpers/volume-representation-params";
import { PluginStateObject } from "../../../mol-plugin-state/objects";
import { StateTransforms } from "../../../mol-plugin-state/transforms";
import { VolumeRepresentation3D } from "../../../mol-plugin-state/transforms/representation";
import { PluginBehavior } from "../../../mol-plugin/behavior/behavior";
import { PluginContext } from "../../../mol-plugin/context";
import { Asset } from "../../../mol-util/assets";
import { Color } from "../../../mol-util/color";
import { ColorNames } from "../../../mol-util/color/names";
import { getFileNameInfo } from "../../../mol-util/file-info";
import { Unzip } from "../../../mol-util/zip/zip";
import { objectToArray } from "../new-volumes-and-segmentations/helpers";
import { CreateShapePrimitiveProvider, CreateShapePrimitiveProviderCVSX, VolsegGeometricSegmentation, VolsegShapePrimitivesData } from "../new-volumes-and-segmentations/shape_primitives";
import { AnnotationMetadata, DescriptionData, ShapePrimitiveData } from "../new-volumes-and-segmentations/volseg-api/data";
import { CSVXUI, CVSX_ANNOTATIONS_FILE_TAG, CVSX_GEOMETRIC_SEGMENTATION_FILE, CVSX_LATTICE_SEGMENTATION_VISUAL_TAG, CVSX_VOLUME_VISUAL_TAG } from "./cvsx";
import { VisualizeStaticQueryZipUI } from "./ui";

// TODO: where VolsegEntryData is used
// it is used in index.ts of new-volumes-and-segmentations
// as entryNode.data
// how to instantiate it from app.ts?
// similar to MVS
// const mvsData = MVSData.fromMVSJ(data);

export const VisualizeStaticQueryZip = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'visualize-static-query-zip',
    category: 'misc',
    display: {
        name: 'Visualize Static Query Zip',
        description: 'Visualize Static Query Zip'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register() {
            this.ctx.customStructureControls.set('visualize-static-query-zip', CSVXUI as any);
        }
        unregister() {
            this.ctx.customStructureControls.delete('visualize-static-query-zip');
        }
    }
});


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

function findNodesByRef(plugin: PluginContext, ref: string) {
    return plugin.state.data.selectQ(q => q.byRef(ref).subtree())[0];
}

function findNodesByTags(plugin: PluginContext, ...tags: string[]) {
    return plugin.state.data.selectQ(q => {
        let builder = q.root.subtree();
        for (const tag of tags) builder = builder.withTag(tag);
        return builder;
    });
}

export async function processCvsxAnnotationsFile(file: Asset.File, plugin: PluginContext) {
    const info = getFileNameInfo(file.file?.name ?? '');
    const isBinary = plugin.dataFormats.binaryExtensions.has(info.ext);
    const { data } = await plugin.builders.data.readFile({ file, isBinary }, { tags: [CVSX_ANNOTATIONS_FILE_TAG] });
    const parsedData: AnnotationMetadata = JSON.parse(data.cell!.obj!.data as string);
    return parsedData;
}

export async function processCvsxGeometricSegmentationFile(file: Asset.File, plugin: PluginContext, annotations?: AnnotationMetadata) {
    const info = getFileNameInfo(file.file?.name ?? '');
    const isBinary = plugin.dataFormats.binaryExtensions.has(info.ext);
    const { data } = await plugin.builders.data.readFile({ file, isBinary }, { tags: [CVSX_GEOMETRIC_SEGMENTATION_FILE] });
    const parsedData: ShapePrimitiveData = JSON.parse(data.cell!.obj!.data as string);

    const segmentationId = annotations?.annotations[0].segmentation_id;
    const segmentAnnotations = annotations?.annotations;
    const descriptions = (objectToArray(annotations?.descriptions) as DescriptionData[]);
    const update = plugin.build().toRoot(); 
    for (const primitiveData of parsedData.shape_primitive_list) {
        update.to(data.ref)
            //     // TODO: can provide a single description and a single segment annotation
            .apply(CreateShapePrimitiveProviderCVSX, { segmentId: primitiveData.id, descriptions: descriptions, segmentAnnotations: segmentAnnotations, segmentationId: segmentationId })
            //         // TODO: shape representation 3d could have no alpha
            .apply(StateTransforms.Representation.ShapeRepresentation3D, { alpha: 0.5 }, { tags: ['geometric-segmentation-visual', '0', `segment-${primitiveData.id}`] })
            // .commit();
    }
    await update.commit();
    // debugger;
}

export async function processCvsxFile(file: Asset.File, plugin: PluginContext, format: string, annotations?: AnnotationMetadata) {
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
    if (format === 'dscif') {
        const parsedOne: Volume = parsed.volumes[0].data;
        const update = plugin.build().toRoot();

        const newParams = createVolumeRepresentationParams(plugin, parsedOne, {
            type: 'isosurface',
            typeParams: { isoValue: Volume.IsoValue.relative(1.5), alpha: 0.2 },
            color: 'uniform',
            colorParams: { value: ColorNames.black }
        })

        const volumeRepresentation3D = await update
            .to(parsed.volumes[0])
            .apply(StateTransforms.Representation.VolumeRepresentation3D, newParams, { tags: [CVSX_VOLUME_VISUAL_TAG] })
            .commit();
    } else if (format === 'segcif') {
        const parsedOne: Volume = parsed.volumes[0].data;
        const update = plugin.build().toRoot();
        if (annotations) {
            const params = createVolumeRepresentationParams(plugin, parsedOne, {
                type: 'segment',
                typeParams: { tryUseGpu: false },
                color: 'volume-segment',
                colorParams: {
                    palette: createPaletteCVSX(annotations)
                }
            })
            const volumeRepresentation3D = await update
                .to(parsed.volumes[0])
                .apply(StateTransforms.Representation.VolumeRepresentation3D, params, { tags: CVSX_LATTICE_SEGMENTATION_VISUAL_TAG })
                .commit();
        } else {
            const params = createVolumeRepresentationParams(plugin, parsedOne, {
                type: 'segment',
                typeParams: { tryUseGpu: false },
                color: 'volume-segment',
            })
            const volumeRepresentation3D = await update
                .to(parsed.volumes[0])
                .apply(StateTransforms.Representation.VolumeRepresentation3D, params, { tags: CVSX_LATTICE_SEGMENTATION_VISUAL_TAG })
                .commit();
        }
    }

    // if (visuals) {
    //     const visuals = await provider.visuals?.(plugin, parsed);
    //     const visualsByFormats: OutputByFormat = {
    //         format: format,
    //         visuals: visuals,
    //         parsed: parsed
    //     }
    //     return visualsByFormats;
    // }
};


