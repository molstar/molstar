import { Button, IconButton } from "../../mol-plugin-ui/controls/common";
import { sleep } from "../../mol-util/sleep";
import { actionSelectSegment, actionToggleAllSegments, actionToggleSegment } from "./common";
import { WaitingButton } from "./new-volumes-and-segmentations/ui";
import { createSegmentKey, parseSegmentKey } from "./new-volumes-and-segmentations/volseg-api/utils";
import * as Icons from '../../mol-plugin-ui/controls/icons';
import { useBehavior } from "../../mol-plugin-ui/hooks/use-behavior";
import { VolsegEntryData } from "./new-volumes-and-segmentations/entry-root";
import { CVSXStateModel } from "./cvsx-extension/cvsx";

export function DescriptionsList({ model, targetSegmentationId, targetKind }: { model: VolsegEntryData | CVSXStateModel, targetSegmentationId: string, targetKind: 'lattice' | 'mesh' | 'primitive' }) {
    const state = useBehavior(model.currentState);
    const currentTimeframe = useBehavior(model.currentTimeframe);
    const metadata = useBehavior(model.metadata);
    // const allDescriptions = metadata!.allDescriptions;
    // is there a method in metadata that gets all descriptions for target segmentation id and kind?
    const allDescriptionsForSegmentationId = metadata!.getAllDescriptionsForSegmentationAndTimeframe(
        targetSegmentationId,
        targetKind,
        currentTimeframe
    );
    const parsedSelectedSegmentKey = parseSegmentKey(state.selectedSegment);
    const { segmentId, segmentationId, kind } = parsedSelectedSegmentKey;
    const selectedSegmentDescriptions = model.metadata.value!.getSegmentDescription(segmentId, segmentationId, kind);
    // NOTE: for now single description
    const selectedSegmentDescription = selectedSegmentDescriptions ? selectedSegmentDescriptions[0] : undefined;
    const visibleSegmentKeys = state.visibleSegments.map(seg => seg.segmentKey);
    console.log(visibleSegmentKeys);
    // const visibleModels = state.visibleModels.map(model => model.pdbId);
    // const allPdbs = model.pdbs;

    
    return <>{allDescriptionsForSegmentationId.length > 0 && <>
        <WaitingButton onClick={async () => { await sleep(20); await actionToggleAllSegments(model, targetSegmentationId, targetKind); }} style={{ marginTop: 1 }}>
            Toggle All segments
        </WaitingButton>
        <div style={{ maxHeight: 200, overflow: 'hidden', overflowY: 'auto', marginBlock: 1 }}>
            {allDescriptionsForSegmentationId.map(d => {
                if (d.target_kind === 'entry' || !d.target_id || d.is_hidden === true) return;
                // NOTE: if time is a single number
                if (d.time && Number.isFinite(d.time) && d.time !== currentTimeframe) return;
                // NOTE: if time is array
                if (d.time && Array.isArray(d.time) && d.time.every(i => Number.isFinite(i)) && !(d.time as number[]).includes(currentTimeframe)) return;
                const segmentKey = createSegmentKey(d.target_id.segment_id, d.target_id.segmentation_id, d.target_kind);
                // How it was before
                // return <div style={{ display: 'flex', marginBottom: 1 }} key={`${d.target_id?.segment_id}:${d.target_id?.segmentation_id}:${d.target_kind}`}
                return <div className='msp-flex-row' style={{ marginTop: '1px' }} key={`${d.target_id?.segment_id}:${d.target_id?.segmentation_id}:${d.target_kind}`}
                // onMouseEnter={() => model.actionHighlightSegment(segmentKey)}
                // onMouseLeave={() => model.actionHighlightSegment()}
                >

                    <Button noOverflow flex onClick={() => actionSelectSegment(model, d !== selectedSegmentDescription ? segmentKey : undefined)}
                        style={{
                            fontWeight: d.target_id.segment_id === selectedSegmentDescription?.target_id?.segment_id
                                && d.target_id.segmentation_id === selectedSegmentDescription?.target_id.segmentation_id
                                ? 'bold' : undefined, textAlign: 'left'
                        }}>
                        <div title={d.name ?? 'Unnamed segment'} style={{ maxWidth: 240, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
                            {d.name ?? 'Unnamed segment'} ({d.target_id?.segment_id}) ({d.target_id?.segmentation_id})
                        </div>
                    </Button>
                    {/* <IconButton svg={Icons.WarningSvg} title={'Remove description'}
                        onClick={() => {
                            model.removeDescription(d.id);
                            // NOTE: assumes single description per segment
                            // model.actionToggleSegment(segmentKey);
                        }} />
                    <IconButton svg={Icons.WarningSvg} title={'Remove segment annotation'}
                        onClick={() => {
                            // there is just one segment annotation for each segment
                            // pick it by id
                            model.removeSegmentAnnotation(d.target_id!.segment_id, d.target_id!.segmentation_id, d.target_kind!);
                            // NOTE: assumes single description per segment
                            // model.actionToggleSegment(segmentKey);
                        }} /> */}
                    <IconButton svg={visibleSegmentKeys.includes(segmentKey) ? Icons.VisibilityOutlinedSvg : Icons.VisibilityOffOutlinedSvg}
                        title={visibleSegmentKeys.includes(segmentKey) ? 'Hide segment' : 'Show segment'}
                        onClick={() => actionToggleSegment(model, segmentKey)} />
                </div>;
            }
            )}
        </div>
    </>}</>;
}