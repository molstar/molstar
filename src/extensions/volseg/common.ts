import { Volume } from '../../mol-model/volume';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { PluginContext } from '../../mol-plugin/context';
import { Asset } from '../../mol-util/assets';
import { Choice } from '../../mol-util/param-choice';
import { VolsegEntryData } from './new-volumes-and-segmentations/entry-root';
import { SEGMENT_VISUAL_TAG } from './new-volumes-and-segmentations/entry-segmentation';
import { DescriptionData, ParsedSegmentKey } from './new-volumes-and-segmentations/volseg-api/data';
import { createSegmentKey, parseSegmentKey } from './new-volumes-and-segmentations/volseg-api/utils';

export async function parseCVSXJSON(rawFile: [string, Uint8Array], plugin: PluginContext) {
    const [fn, filedata] = rawFile;
    const file = new File([filedata], fn);
    const f = Asset.File(file);
    const asset = plugin.managers.asset.resolve(f, 'string');
    const runnedAsset = (await asset.run()).data;
    const parsedData = JSON.parse(runnedAsset);
    return parsedData;
}

// export async function actionToggleAllSegments(model: VolsegEntryData, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive') {
//     const currentTimeframe = model.currentTimeframe.value;
//     const current = model.currentState.value.visibleSegments.map(seg => seg.segmentKey);
//     const currentForThisSegmentation = current.filter(c =>
//         parseSegmentKey(c).segmentationId === segmentationId &&
//         parseSegmentKey(c).kind === kind
//     );
//     const currentForOtherSegmentations = current.filter(item => currentForThisSegmentation.indexOf(item) < 0);
//     if (currentForThisSegmentation.length !== model.metadata.value!.getAllSegmentAnotationsForSegmentationAndTimeframe(segmentationId, kind, currentTimeframe).length) {
//         const allSegmentKeysFromSegmentation = model.metadata.value!.getAllSegmentAnotationsForSegmentationAndTimeframe(segmentationId, kind, currentTimeframe).map(a =>
//             createSegmentKey(a.segment_id, a.segmentation_id, a.segment_kind)
//         );
//         const allSegmentKeys = [...allSegmentKeysFromSegmentation, ...currentForOtherSegmentations];
//         await actionShowSegments(allSegmentKeys, model);
//     } else {
//         await actionShowSegments(currentForOtherSegmentations, model);
//     }
// }



export async function actionToggleAllFilteredSegments(model: VolsegEntryData, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive', filteredDescriptions: DescriptionData[]) {
    const currentTimeframe = model.currentTimeframe.value;
    // This is currently visible
    const current = model.currentState.value.visibleSegments.map(seg => seg.segmentKey);
    // this is currently visible for this segmentation
    const currentForThisSegmentation = current.filter(c =>
        parseSegmentKey(c).segmentationId === segmentationId &&
        parseSegmentKey(c).kind === kind
    );
    // const currentForThisSegmentation = filteredDescriptions;
    // this is currently visible for other segmentations
    const currentForOtherSegmentations = current.filter(item => currentForThisSegmentation.indexOf(item) < 0);
    // length of currently visible for this segmentation !== number of all segments for this segmentation

    // This checks if
    // if (currentForThisSegmentation.length !== model.metadata.value!.getAllSegmentAnotationsForSegmentationAndTimeframe(segmentationId, kind, currentTimeframe).length) {
    const allFilteredSegmentKeys = filteredDescriptions.map(d => {
        // there are three types that are acceptable
        // should add check that this is not entry
        if (d.target_kind !== 'entry') {
            return createSegmentKey(d.target_id!.segment_id, d.target_id!.segmentation_id, d.target_kind);
        }
    }).filter(Boolean);
    const allSegmentKeysFromSegmentation = model.metadata.value!.getAllDescriptionsForSegmentationAndTimeframe(segmentationId, kind, currentTimeframe).map(d => {
        if (d.target_kind !== 'entry') {
            return createSegmentKey(d.target_id!.segment_id, d.target_id!.segmentation_id, d.target_kind);
        }
    }).filter(Boolean);
    if (currentForThisSegmentation.length !== filteredDescriptions.length) {

        console.log(filteredDescriptions);
        const allSegmentKeys = [...allFilteredSegmentKeys, ...currentForOtherSegmentations];
        await actionShowSegments(allSegmentKeys as string[], model);
    } else {
        // TODO: add to currentForOtherSegmentations
        // segments from this segmentation, but not filtered
        const descriptionsFromThisSegmentationNotFiltered = allSegmentKeysFromSegmentation.filter(i => allFilteredSegmentKeys.indexOf(i) < 0);
        console.log(descriptionsFromThisSegmentationNotFiltered);
        debugger;
        await actionShowSegments(currentForOtherSegmentations, model);
    }
}

export async function actionShowSegments(segmentKeys: string[], model: VolsegEntryData) {
    const allExistingLatticeSegmentationIds = model.metadata.value!.raw.grid.segmentation_lattices!.segmentation_ids;
    const allExistingMeshSegmentationIds = model.metadata.value!.raw.grid.segmentation_meshes!.segmentation_ids;
    const allExistingGeometricSegmentationIds = model.metadata.value!.raw.grid.geometric_segmentation!.segmentation_ids;
    if (segmentKeys.length === 0) {
        for (const id of allExistingLatticeSegmentationIds) {
            // await model.latticeSegmentationData.showSegments([], id);
            await showSegments([], id, 'lattice', model);
        }
        for (const id of allExistingMeshSegmentationIds) {
            // await model.meshSegmentationData.showSegments([], id);
            await showSegments([], id, 'mesh', model);
        }
        for (const id of allExistingGeometricSegmentationIds) {
            // await model.geometricSegmentationData.showSegments([], id);
            await showSegments([], id, 'primitive', model);
        }
    }
    const parsedSegmentKeys = segmentKeys.map(
        k => parseSegmentKey(k)
    );
    // LATTICES PART
    _actionShowSegments(parsedSegmentKeys, allExistingLatticeSegmentationIds, 'lattice', model);
    // MESHES PART
    _actionShowSegments(parsedSegmentKeys, allExistingMeshSegmentationIds, 'mesh', model);
    // GEOMETRIC SEGMENTATION PAR
    _actionShowSegments(parsedSegmentKeys, allExistingGeometricSegmentationIds, 'primitive', model);

    await model.updateStateNode({ visibleSegments: segmentKeys.map(s => ({ segmentKey: s })) });
}

export async function _actionShowSegments(parsedSegmentKeys: ParsedSegmentKey[], existingSegmentationIds: string[], kind: 'mesh' | 'lattice' | 'primitive', model: VolsegEntryData) {
    const segmentKeys = parsedSegmentKeys.filter(k => k.kind === kind);
    // const segmentIds = segmentKeys.map(k => k.segmentId);
    const SegmentationIdsToSegmentIds = new Map<string, number[]>;
    for (const key of segmentKeys) {
        if (!SegmentationIdsToSegmentIds.has(key.segmentationId)) {
            SegmentationIdsToSegmentIds.set(key.segmentationId, [key.segmentId]);
        } else {
            // have this key, add segmentId
            const currentSegmentationIds = SegmentationIdsToSegmentIds.get(key.segmentationId);
            SegmentationIdsToSegmentIds.set(key.segmentationId, [...currentSegmentationIds!, key.segmentId]);
        }
    }
    for (const id of existingSegmentationIds) {
        if (!SegmentationIdsToSegmentIds.has(id)) {
            SegmentationIdsToSegmentIds.set(id, []);
        }
    }
    console.log(SegmentationIdsToSegmentIds);
    const promises: Promise<void>[] = [];
    SegmentationIdsToSegmentIds.forEach((value, key) => {
        // TODO: refactor these three showSegments - they seem to not use any submodel parameters
        // Take three
        promises.push(showSegments(value, key, kind, model));
        // if (kind === 'lattice') promises.push(showSegments(value, key, 'lattice'));
        // else if (kind === 'mesh') promises.push(showSegments(value, key));
        // else if (kind === 'primitive') promises.push(showSegments(value, key));
    });
    await Promise.all(promises);
}


export function findNodesByTags(plugin: PluginContext, ...tags: string[]) {
    return plugin.state.data.selectQ(q => {
        let builder = q.root.subtree();
        for (const tag of tags) builder = builder.withTag(tag);
        return builder;
    });
}

export function findNodesByRef(plugin: PluginContext, ref: string) {
    // return this.plugin.state.data.selectQ(q => q.byRef(ref).subtree())[0];
    return plugin.state.data.selectQ(q => q.byRef(ref).subtree())[0];
}

async function selectSegment(model: VolsegEntryData, segmentId?: number, segmentationId?: string) {
    if (segmentId === undefined || segmentId < 0 || !segmentationId) return;
    const segmentLoci = makeLoci([segmentId], segmentationId, model);
    if (!segmentLoci) return;
    model.plugin.managers.interactivity.lociSelects.select(segmentLoci, false);
}

function makeLoci(segments: number[], segmentationId: string, model: VolsegEntryData) {
    const vis = findNodesByTags(model.plugin, SEGMENT_VISUAL_TAG, segmentationId)[0];
    if (!vis) return undefined;
    const repr = vis.obj?.data.repr;
    const wholeLoci = repr.getAllLoci()[0];
    if (!wholeLoci || !Volume.Segment.isLoci(wholeLoci)) return undefined;
    return { loci: Volume.Segment.Loci(wholeLoci.volume, segments), repr: repr };
}

async function showSegments(segmentIds: number[], segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive', model: VolsegEntryData) {
    if (kind === 'lattice') {
        const repr = findNodesByTags(model.plugin, SEGMENT_VISUAL_TAG, segmentationId)[0];
        if (!repr) return;
        const selectedSegmentKey = model.currentState.value.selectedSegment;
        const parsedSelectedSegmentKey = parseSegmentKey(selectedSegmentKey);
        const selectedSegmentId = parsedSelectedSegmentKey.segmentId;
        const selectedSegmentSegmentationId = parsedSelectedSegmentKey.segmentationId;
        const mustReselect = segmentIds.includes(selectedSegmentId) && !repr.params?.values.type.params.segments.includes(selectedSegmentId);
        const update = model.plugin.build().toRoot();
        update.to(repr).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.segments = segmentIds; });
        await update.commit();
        if (mustReselect) {
            await selectSegment(model, selectedSegmentId, selectedSegmentSegmentationId);
        }
    } else if (kind === 'mesh') {
        const segmentsToShow = new Set(segmentIds);

        const visuals = findNodesByTags(model.plugin, 'mesh-segment-visual', segmentationId);
        // debugger;
        for (const visual of visuals) {
            const theTag = visual.obj?.tags?.find(tag => tag.startsWith('segment-'));
            if (!theTag) continue;
            const id = parseInt(theTag.split('-')[1]);
            const visibility = segmentsToShow.has(id);
            setSubtreeVisibility(model.plugin.state.data, visual.transform.ref, !visibility); // true means hide, ¯\_(ツ)_/¯
            segmentsToShow.delete(id);
        }
    } else if (kind === 'primitive') {
        const segmentsToShow = new Set(segmentIds);

        // This will select all segments of that segmentation
        const visuals = findNodesByTags(model.plugin, 'geometric-segmentation-visual', segmentationId);
        for (const visual of visuals) {
            const theTag = visual.obj?.tags?.find(tag => tag.startsWith('segment-'));
            if (!theTag) continue;
            const id = parseInt(theTag.split('-')[1]);
            const visibility = segmentsToShow.has(id);
            setSubtreeVisibility(model.plugin.state.data, visual.transform.ref, !visibility); // true means hide, ¯\_(ツ)_/¯
            segmentsToShow.delete(id);
        }
    } else {
        throw new Error('Not supported');
    }
}

export async function actionSelectSegment(model: VolsegEntryData, segmentKey?: string) {
    if (segmentKey !== undefined && model.currentState.value.visibleSegments.find(s => s.segmentKey === segmentKey) === undefined) {
        // first make the segment visible if it is not
        await actionToggleSegment(model, segmentKey);
    }
    await model.updateStateNode({ selectedSegment: segmentKey });
}

export async function actionToggleSegment(model: VolsegEntryData, segmentKey: string) {
    const current = model.currentState.value.visibleSegments.map(seg => seg.segmentKey);
    if (current.includes(segmentKey)) {
        await actionShowSegments(current.filter(s => s !== segmentKey), model);
    } else {
        await actionShowSegments([...current, segmentKey], model);
    }
}
export const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR', idr: 'IDR', pdbe: 'PDBe', custom: 'CUSTOM' }, 'emdb');
export type Source = Choice.Values<typeof SourceChoice>;

