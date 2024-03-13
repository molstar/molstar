import { Button, ExpandGroup, IconButton, TextInput } from '../../mol-plugin-ui/controls/common';
import { sleep } from '../../mol-util/sleep';
import { actionSelectSegment, actionToggleAllFilteredSegments, actionToggleSegment } from './common';
import { WaitingButton } from './new-volumes-and-segmentations/ui';
import { createSegmentKey, parseSegmentKey } from './new-volumes-and-segmentations/volseg-api/utils';
import * as Icons from '../../mol-plugin-ui/controls/icons';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';
import { VolsegEntryData } from './new-volumes-and-segmentations/entry-root';
import { CVSXStateModel } from './cvsx-extension/cvsx';
import Markdown from 'react-markdown';
import { capitalize } from '../../mol-util/string';
import { useState } from 'react';
import { DescriptionData } from './new-volumes-and-segmentations/volseg-api/data';


export const MetadataTextFilter = ({ setFilteredDescriptions, descriptions, model }: { setFilteredDescriptions: any, descriptions: DescriptionData[], model: VolsegEntryData }) => {
    const [text, setText] = useState('');

    function filterDescriptions(keyword: string) {
        return model.metadata!.value!.filterDescriptionsBasedOnKeywordInMetadata(descriptions, keyword);
    }

    return (
        // <View style={{ padding: 10 }}>
        <TextInput
            style={{ order: 1, flex: '1 1 auto', minWidth: 0, marginBlock: 1 }} className='msp-form-control'
            value={text}
            placeholder="Type keyword to filter segments..."
            // this will just set text, need to filter metadata based on text
            onChange={newText => {
                setText(newText);
                const filteredDescriptions = filterDescriptions(newText);
                debugger;
                setFilteredDescriptions(filteredDescriptions);
            }}
        />
    );
};

export function DescriptionsList({ model, targetSegmentationId, targetKind }: { model: VolsegEntryData, targetSegmentationId: string, targetKind: 'lattice' | 'mesh' | 'primitive' }) {
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

    const [filteredDescriptions, setFilteredDescriptions] = useState(allDescriptionsForSegmentationId);

    const parsedSelectedSegmentKey = parseSegmentKey(state.selectedSegment);
    const { segmentId, segmentationId, kind } = parsedSelectedSegmentKey;
    const selectedSegmentDescriptions = model.metadata.value!.getSegmentDescription(segmentId, segmentationId, kind);
    // NOTE: for now single description
    const selectedSegmentDescription = selectedSegmentDescriptions ? selectedSegmentDescriptions[0] : undefined;
    const visibleSegmentKeys = state.visibleSegments.map(seg => seg.segmentKey);
    console.log(visibleSegmentKeys);
    // const visibleModels = state.visibleModels.map(model => model.pdbId);
    // const allPdbs = model.pdbs;

    return <>
        <MetadataTextFilter setFilteredDescriptions={setFilteredDescriptions} descriptions={allDescriptionsForSegmentationId} model={model}></MetadataTextFilter>
        {filteredDescriptions.length > 0 && <>
            <WaitingButton onClick={async () => { await sleep(20); await actionToggleAllFilteredSegments(model, targetSegmentationId, targetKind, filteredDescriptions); }} style={{ marginTop: 1 }}>
                Toggle All segments
            </WaitingButton>
            <div style={{ maxHeight: 200, overflow: 'hidden', overflowY: 'auto', marginBlock: 1 }}>
                {filteredDescriptions.map(d => {
                    if (d.target_kind === 'entry' || !d.target_id || d.is_hidden === true) return;
                    // NOTE: if time is a single number
                    if (d.time && Number.isFinite(d.time) && d.time !== currentTimeframe) return;
                    // NOTE: if time is array
                    if (d.time && Array.isArray(d.time) && d.time.every(i => Number.isFinite(i)) && !(d.time as number[]).includes(currentTimeframe)) return;
                    const segmentKey = createSegmentKey(d.target_id.segment_id, d.target_id.segmentation_id, d.target_kind);
                    // How it was before
                    // return <div style={{ display: 'flex', marginBottom: 1 }} key={`${d.target_id?.segment_id}:${d.target_id?.segmentation_id}:${d.target_kind}`}
                    return <div className='msp-flex-row' style={{ marginTop: '1px' }} key={`${d.target_id?.segment_id}:${d.target_id?.segmentation_id}:${d.target_kind}`}
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

export function SelectedSegmentDescription({ model, targetSegmentationId, targetKind }: { model: VolsegEntryData | CVSXStateModel, targetSegmentationId: string, targetKind: 'lattice' | 'mesh' | 'primitive' }) {
    const state = useBehavior(model.currentState);
    const currentTimeframe = useBehavior(model.currentTimeframe);
    const metadata = useBehavior(model.metadata);
    // const allDescriptions = metadata!.allDescriptions;
    // is there a method in metadata that gets all descriptions for target segmentation id and kind?
    const anyDescriptions = metadata!.allDescriptions.length > 0;
    const parsedSelectedSegmentKey = parseSegmentKey(state.selectedSegment);
    const { segmentId, segmentationId, kind } = parsedSelectedSegmentKey;
    const selectedSegmentDescriptions = model.metadata.value!.getSegmentDescription(segmentId, segmentationId, kind);
    // NOTE: for now single description
    const selectedSegmentDescription = selectedSegmentDescriptions ? selectedSegmentDescriptions[0] : undefined;
    const visibleSegmentKeys = state.visibleSegments.map(seg => seg.segmentKey);
    console.log(visibleSegmentKeys);
    return <>{
        anyDescriptions && <ExpandGroup header='Selected segment descriptions' initiallyExpanded>
            <div style={{ paddingTop: 4, paddingRight: 8, maxHeight: 300, overflow: 'hidden', overflowY: 'auto' }}>
                {!selectedSegmentDescription && 'No segment selected'}
                {selectedSegmentDescription &&
                    selectedSegmentDescription.target_kind !== 'entry' &&
                    selectedSegmentDescription.target_id &&
                    <b>Segment {selectedSegmentDescription.target_id.segment_id} from segmentation {selectedSegmentDescription.target_id.segmentation_id}:<br />{selectedSegmentDescription.name ?? 'Unnamed segment'}</b>}
                {selectedSegmentDescription && selectedSegmentDescription.description && selectedSegmentDescription.description.format === 'markdown' &&
                    <>
                        <br />
                        <br />
                        <b>Description: </b>
                        <Markdown skipHtml>{selectedSegmentDescription.description.text}</Markdown>
                    </>}
                {selectedSegmentDescription && selectedSegmentDescription.description && selectedSegmentDescription.description.format === 'text' &&
                    <>
                        <br />
                        <br />
                        <b>Description: </b>
                        <p>{selectedSegmentDescription.description.text}</p>
                    </>}
                {selectedSegmentDescription?.external_references?.map(ref => {
                    // if (description.target_kind === 'entry' || !description.target_id) return;
                    return <p key={ref.id} style={{ marginTop: 4 }}>
                        {/* <small>{ref.resource}:{ref.accession}</small><br /> */}
                        {ref.url ? <a href={ref.url}>{ref.resource}:{ref.accession}</a> :
                            <small>{ref.resource}:{ref.accession}</small>}
                        <br />
                        <b>{capitalize(ref.label ? ref.label : '')}</b><br />
                        {ref.description}
                    </p>;
                }
                )}
            </div>
        </ExpandGroup>
    }</>;

}