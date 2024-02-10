/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { useCallback, useEffect, useRef, useState } from 'react';

import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button, ControlRow, ExpandGroup, IconButton } from '../../mol-plugin-ui/controls/common';
import * as Icons from '../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import { Slider } from '../../mol-plugin-ui/controls/slider';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';
import { UpdateTransformControl } from '../../mol-plugin-ui/state/update-transform';
import { PluginContext } from '../../mol-plugin/context';
import { shallowEqualArrays } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { sleep } from '../../mol-util/sleep';

import { VolsegEntry, VolsegEntryData } from './entry-root';
import { SimpleVolumeParams, SimpleVolumeParamValues } from './entry-volume';
import { VolsegGlobalState, VolsegGlobalStateData, VolsegGlobalStateParams } from './global-state';
import { isDefined } from './helpers';
import { ProjectDataParamsValues } from './transformers';
import { StateObjectCell } from '../../mol-state';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { createSegmentKey, parseSegmentKey } from './volseg-api/utils';
import Markdown from 'react-markdown';


interface VolsegUIData {
    globalState?: VolsegGlobalStateData,
    availableNodes: VolsegEntry[],
    activeNode?: VolsegEntry,
}
namespace VolsegUIData {
    export function changeAvailableNodes(data: VolsegUIData, newNodes: VolsegEntry[]): VolsegUIData {
        const newActiveNode = newNodes.length > data.availableNodes.length ?
            newNodes[newNodes.length - 1]
            : newNodes.find(node => node.data.ref === data.activeNode?.data.ref) ?? newNodes[0];
        return { ...data, availableNodes: newNodes, activeNode: newActiveNode };
    }
    export function changeActiveNode(data: VolsegUIData, newActiveRef: string): VolsegUIData {
        const newActiveNode = data.availableNodes.find(node => node.data.ref === newActiveRef) ?? data.availableNodes[0];
        return { ...data, availableNodes: data.availableNodes, activeNode: newActiveNode };
    }
    export function equals(data1: VolsegUIData, data2: VolsegUIData) {
        return shallowEqualArrays(data1.availableNodes, data2.availableNodes) && data1.activeNode === data2.activeNode && data1.globalState === data2.globalState;
    }
}

export class VolsegUI extends CollapsableControls<{}, { data: VolsegUIData }> {
    protected defaultState(): CollapsableState & { data: VolsegUIData } {
        return {
            header: 'Volume & Segmentation',
            isCollapsed: true,
            brand: { accent: 'orange', svg: Icons.ExtensionSvg },
            data: {
                globalState: undefined,
                availableNodes: [],
                activeNode: undefined,
            }
        };
    }
    protected renderControls(): JSX.Element | null {
        return <VolsegControls plugin={this.plugin} data={this.state.data} setData={d => this.setState({ data: d })} />;
    }
    componentDidMount(): void {
        this.setState({ isHidden: true, isCollapsed: false });
        this.subscribe(this.plugin.state.data.events.changed, e => {
            const nodes = e.state.selectQ(q => q.ofType(VolsegEntry)).map(cell => cell?.obj).filter(isDefined);
            const isHidden = nodes.length === 0;
            const newData = VolsegUIData.changeAvailableNodes(this.state.data, nodes);
            if (!this.state.data.globalState?.isRegistered()) {
                const globalState = e.state.selectQ(q => q.ofType(VolsegGlobalState))[0]?.obj?.data;
                if (globalState) newData.globalState = globalState;
            }
            if (!VolsegUIData.equals(this.state.data, newData) || this.state.isHidden !== isHidden) {
                this.setState({ data: newData, isHidden: isHidden });
            }
        });
    }
}


function VolsegControls({ plugin, data, setData }: { plugin: PluginContext, data: VolsegUIData, setData: (d: VolsegUIData) => void }) {
    const entryData = data.activeNode?.data;
    if (!entryData) {
        return <p>No data!</p>;
    }
    if (!data.globalState) {
        return <p>No global state!</p>;
    }

    const params = {
        /** Reference to the active VolsegEntry node */
        entry: PD.Select(data.activeNode!.data.ref, data.availableNodes.map(entry => [entry.data.ref, entry.data.entryId]))
    };
    const values: PD.Values<typeof params> = {
        entry: data.activeNode!.data.ref,
    };

    const globalState = useBehavior(data.globalState.currentState);

    return <>
        <ParameterControls params={params} values={values} onChangeValues={next => setData(VolsegUIData.changeActiveNode(data, next.entry))} />

        <TimeFrameSlider entryData={entryData} />
        <ExpandGroup header='Global options'>
            <WaitingParameterControls params={VolsegGlobalStateParams} values={globalState} onChangeValues={async next => await data.globalState?.updateState(plugin, next)} />
        </ExpandGroup>

        <VolsegEntryControls entryData={entryData} key={entryData.ref} />
    </>;
}

function VolsegEntryControls({ entryData }: { entryData: VolsegEntryData }) {
    const state = useBehavior(entryData.currentState);

    const allDescriptions = entryData.metadata.allDescriptions;
    const parsedSelectedSegmentKey = parseSegmentKey(state.selectedSegment);
    const { segmentId, segmentationId, kind } = parsedSelectedSegmentKey;
    const selectedSegmentDescriptions = entryData.metadata.getSegment(segmentId, segmentationId, kind);
    // NOTE: for now single description
    const selectedSegmentDescription = selectedSegmentDescriptions ? selectedSegmentDescriptions[0] : undefined;
    const visibleSegmentKeys = state.visibleSegments.map(seg => seg.segmentKey);
    console.log(visibleSegmentKeys);
    const visibleModels = state.visibleModels.map(model => model.pdbId);
    const allPdbs = entryData.pdbs;

    const currentTimeframe = useBehavior(entryData.currentTimeframe);
    console.log('Current timframe is: ', currentTimeframe);
    return <>
        {/* Title */}
        <div style={{ fontWeight: 'bold', padding: 8, paddingTop: 6, paddingBottom: 4, overflow: 'hidden' }}>
            {entryData.metadata.raw.annotation?.name ?? 'Unnamed Annotation'}
        </div>

        {/* Fitted models */}
        {allPdbs.length > 0 && <ExpandGroup header='Fitted models in PDB' initiallyExpanded>
            {allPdbs.map(pdb =>
                <WaitingButton key={pdb} onClick={() => entryData.actionShowFittedModel(visibleModels.includes(pdb) ? [] : [pdb])}
                    style={{ fontWeight: visibleModels.includes(pdb) ? 'bold' : undefined, textAlign: 'left', marginTop: 1 }}>
                    {pdb}
                </WaitingButton>
            )}
        </ExpandGroup>}

        {/* Volume */}
        <VolumeControls entryData={entryData} />

        {allDescriptions.length > 0 && <ExpandGroup header='Segmentation data' initiallyExpanded>
            {/* Segment opacity slider */}
            <ControlRow label='Opacity' control={
                <WaitingSlider min={0} max={1} value={state.segmentOpacity} step={0.05} onChange={async v => await entryData.actionSetOpacity(v)} />
            } />

            {/* Segment toggles */}
            {allDescriptions.length > 0 && <>
                <WaitingButton onClick={async () => { await sleep(20); await entryData.actionToggleAllSegments(); }} style={{ marginTop: 1 }}>
                    Toggle All segments
                </WaitingButton>
                <div style={{ maxHeight: 200, overflow: 'hidden', overflowY: 'auto', marginBlock: 1 }}>
                    {allDescriptions.map(d => {
                        if (d.target_kind === 'entry' || !d.target_id) return;
                        // NOTE: if time is a single number
                        if (d.time && Number.isFinite(d.time) && d.time !== currentTimeframe) return;
                        // NOTE: if time is array
                        if (d.time && Array.isArray(d.time) && d.time.every(i => Number.isFinite(i)) && !(d.time as number[]).includes(currentTimeframe)) return;
                        const segmentKey = createSegmentKey(d.target_id.segment_id, d.target_id.segmentation_id, d.target_kind);
                        return <div style={{ display: 'flex', marginBottom: 1 }} key={`${d.target_id?.segment_id}:${d.target_id?.segmentation_id}:${d.target_kind}`}
                            onMouseEnter={() => entryData.actionHighlightSegment(segmentKey)}
                            onMouseLeave={() => entryData.actionHighlightSegment()}>

                            <Button onClick={() => entryData.actionSelectSegment(d !== selectedSegmentDescription ? segmentKey : undefined)}
                                style={{
                                    fontWeight: d.target_id.segment_id === selectedSegmentDescription?.target_id?.segment_id
                                    && d.target_id.segmentation_id === selectedSegmentDescription?.target_id.segmentation_id
                                        ? 'bold' : undefined, marginRight: 1, flexGrow: 1, textAlign: 'left'
                                }}>
                                <div title={d.name ?? 'Unnamed segment'} style={{ maxWidth: 240, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
                                    {d.name ?? 'Unnamed segment'} ({d.target_id?.segment_id}) ({d.target_id?.segmentation_id})
                                </div>
                            </Button>
                            <IconButton svg={visibleSegmentKeys.includes(segmentKey) ? Icons.VisibilityOutlinedSvg : Icons.VisibilityOffOutlinedSvg}
                                onClick={() => entryData.actionToggleSegment(segmentKey)} />
                        </div>;
                    }
                    )}
                </div>
            </>}
        </ExpandGroup>}

        {/* Segment annotations */}
        {allDescriptions.length > 0 && <ExpandGroup header='Selected segment annotation' initiallyExpanded>
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
                        <small>{ref.resource}:{ref.accession}</small><br />
                        <b>{capitalize(ref.label ? ref.label : '')}</b><br />
                        {ref.description}
                    </p>;
                }
                )}
            </div>
        </ExpandGroup>}
    </>;
}

function TimeFrameSlider({ entryData }: { entryData: VolsegEntryData }) {
    // gets time info from volume
    const timeInfo = entryData.metadata.raw.grid.volumes.time_info;
    const timeInfoStart = timeInfo.start;
    const timeInfoValue = useBehavior(entryData.currentTimeframe);
    const timeInfoEnd = timeInfo.end;
    if (timeInfoEnd === 0) return null;

    return <ControlRow label='Time Frame' control={
        <WaitingSlider min={timeInfoStart} max={timeInfoEnd} value={timeInfoValue} step={1}
            onChange={async v => {
                await entryData.updateProjectData(v);
            }}
        />}
    />;
}

function VolumeChannelControls({ entryData, volume }: { entryData: VolsegEntryData, volume: StateObjectCell<PluginStateObject.Volume.Data> }) {
    const projectDataTransform = volume.transform;

    if (!projectDataTransform) return null;
    const params: ProjectDataParamsValues = projectDataTransform.params;
    const channelId = params.channelId;
    const channelLabel = volume.obj!.label;
    const childRef = entryData.plugin.state.data.tree.children.get(projectDataTransform.ref).toArray()[0];
    const volumeRepresentation3DNode = entryData.findNodesByRef(childRef);
    const transform = volumeRepresentation3DNode.transform;
    if (!transform) return null;
    const volumeValues: SimpleVolumeParamValues = {
        volumeType: transform.state.isHidden ? 'off' : transform.params?.type.name as any,
        opacity: transform.params?.type.params.alpha,
    };

    return <ExpandGroup header={`${channelLabel}`}>
        <WaitingParameterControls params={SimpleVolumeParams} values={volumeValues} onChangeValues={async next => { await sleep(20); await entryData.actionUpdateVolumeVisual(next, channelId, transform); }} />
        <UpdateTransformControl state={entryData.plugin.state.data} transform={transform} customHeader='none' />
    </ExpandGroup>;
}

function VolumeControls({ entryData }: { entryData: VolsegEntryData }) {
    const h = useBehavior(entryData.state.hierarchy);
    if (!h) return null;
    return <>
        {/* <Button onClick={() => { console.log('volume cache, segmentation cache: ', entryData.cachedVolumeTimeframesData, entryData.cachedSegmentationTimeframesData); }}>Get volume and segmentation cache</Button> */}
        <ExpandGroup header='Volume data'>
            {h.volumes.map((v) => {
                const params: ProjectDataParamsValues = v.transform.params;
                return <VolumeChannelControls key={params.channelId} entryData={entryData} volume={v} />;
            })}
        </ExpandGroup>
    </>;
}


type ComponentParams<T extends React.Component<any, any, any> | ((props: any) => JSX.Element)> =
    T extends React.Component<infer P, any, any> ? P : T extends (props: infer P) => JSX.Element ? P : never;

function WaitingSlider({ value, onChange, ...etc }: { value: number, onChange: (value: number) => any } & ComponentParams<Slider>) {
    const [changing, sliderValue, execute] = useAsyncChange(value);

    return <Slider value={sliderValue} disabled={changing} onChange={newValue => execute(onChange, newValue)} {...etc} />;
}

function WaitingButton({ onClick, ...etc }: { onClick: () => any } & ComponentParams<typeof Button>) {
    const [changing, _, execute] = useAsyncChange(undefined);

    return <Button disabled={changing} onClick={() => execute(onClick, undefined)} {...etc}>
        {etc.children}
    </Button>;
}

function WaitingParameterControls<T extends PD.Params>({ values, onChangeValues, ...etc }: { values: PD.ValuesFor<T>, onChangeValues: (values: PD.ValuesFor<T>) => any } & ComponentParams<ParameterControls<T>>) {
    const [changing, currentValues, execute] = useAsyncChange(values);

    return <ParameterControls isDisabled={changing} values={currentValues} onChangeValues={newValue => execute(onChangeValues, newValue)} {...etc} />;
}

function capitalize(text: string) {
    const first = text.charAt(0);
    const rest = text.slice(1);
    return first.toUpperCase() + rest;
}

function useAsyncChange<T>(initialValue: T) {
    const [isExecuting, setIsExecuting] = useState(false);
    const [value, setValue] = useState(initialValue);
    const isMounted = useRef(false);

    useEffect(() => setValue(initialValue), [initialValue]);

    useEffect(() => {
        isMounted.current = true;
        return () => { isMounted.current = false; };
    }, []);

    const execute = useCallback(
        async (func: (val: T) => Promise<any>, val: T) => {
            setIsExecuting(true);
            setValue(val);
            try {
                await func(val);
            } catch (err) {
                if (isMounted.current) {
                    setValue(initialValue);
                }
                throw err;
            } finally {
                if (isMounted.current) {
                    setIsExecuting(false);
                }
            }
        },
        []
    );

    return [isExecuting, value, execute] as const;
}
