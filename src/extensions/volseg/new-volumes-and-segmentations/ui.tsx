/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */
import Popup from 'reactjs-popup';
import 'reactjs-popup/dist/index.css';
import { useCallback, useEffect, useRef, useState } from 'react';
// import { ParamDefinition as PD } from '../../mol-util/param-definition';

import { CollapsableControls, CollapsableState } from '../../../mol-plugin-ui/base';
import { Button, ControlRow, ExpandGroup, IconButton, TextInput } from '../../../mol-plugin-ui/controls/common';
import * as Icons from '../../../mol-plugin-ui/controls/icons';
import { ParameterControls, TextControl } from '../../../mol-plugin-ui/controls/parameters';
import { Slider } from '../../../mol-plugin-ui/controls/slider';
import { useBehavior } from '../../../mol-plugin-ui/hooks/use-behavior';
import { UpdateTransformControl } from '../../../mol-plugin-ui/state/update-transform';
import { PluginContext } from '../../../mol-plugin/context';
import { shallowEqualArrays } from '../../../mol-util';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { sleep } from '../../../mol-util/sleep';

import { VolsegEntry, VolsegEntryData } from './entry-root';
import { SimpleVolumeParams, SimpleVolumeParamValues } from './entry-volume';
import { VolsegGlobalState, VolsegGlobalStateData, VolsegGlobalStateParams } from './global-state';
import { isDefined } from './helpers';
import { ProjectDataParamsValues, ProjectGeometricSegmentationDataParamsValues, ProjectLatticeSegmentationDataParamsValues, ProjectMeshSegmentationDataParamsValues } from './transformers';
import { StateObjectCell } from '../../../mol-state';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { createSegmentKey, parseSegmentKey } from './volseg-api/utils';
import Markdown from 'react-markdown';
import { Asset } from '../../../mol-util/assets';
import { DescriptionData, SegmentAnnotationData } from './volseg-api/data';
import React from "react";
import JSONEditorComponent from './jsoneditor-component';
import { VolsegGeometricSegmentation } from './shape_primitives';
import { VolsegMeshSegmentation } from '../new-meshes/mesh-extension';
import { actionSelectSegment, actionToggleSegment, findNodesByRef } from '../common';
import { CVSXStateModel } from '../cvsx-extension/cvsx';
import { DescriptionsList, EntryDescriptionUI, MetadataTextFilter, SelectedSegmentDescription } from '../common-ui';
import { Script } from '../../../mol-script/script';

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

async function parseJSONwithAnnotationsOrDescriptions(v, entryData: VolsegEntryData) {
    console.log(v.target.files![0]);
    const file = Asset.File(v.target.files![0]);
    const asset = entryData.plugin.managers.asset.resolve(file, 'string');
    const data = (await asset.run()).data;
    const parsedData: DescriptionData[] | SegmentAnnotationData[] = JSON.parse(data);
    return parsedData;
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
    const metadata = useBehavior(entryData.metadata);
    const allDescriptions = entryData.metadata.value!.allDescriptions;
    const entryDescriptions = allDescriptions.filter(d => d.target_kind === 'entry');
    const parsedSelectedSegmentKey = parseSegmentKey(state.selectedSegment);
    const { segmentId, segmentationId, kind } = parsedSelectedSegmentKey;
    const selectedSegmentDescriptions = entryData.metadata.value!.getSegmentDescription(segmentId, segmentationId, kind);
    // NOTE: for now single description
    const selectedSegmentDescription = selectedSegmentDescriptions ? selectedSegmentDescriptions[0] : undefined;
    const visibleSegmentKeys = state.visibleSegments.map(seg => seg.segmentKey);
    console.log(visibleSegmentKeys);
    const visibleModels = state.visibleModels.map(model => model.pdbId);
    const allPdbs = entryData.pdbs;

    let text = 'TEXT';

    const currentTimeframe = useBehavior(entryData.currentTimeframe);
    console.log('Current timframe is: ', currentTimeframe);
    console.log('UI re-rendered');
    const annotationsJson = metadata!.raw.annotation;
    return <>
        {/* Title */}
        <div style={{ fontWeight: 'bold', padding: 8, paddingTop: 6, paddingBottom: 4, overflow: 'hidden' }}>
            {metadata!.raw.annotation?.name ?? 'Unnamed Annotation'}
        </div>
        {entryDescriptions.length > 0 && entryDescriptions.map(e =>
            <EntryDescriptionUI key={e.id} entryDescriptionData={e}></EntryDescriptionUI>)}
        {/* <JSONEditorComponent jsonData={annotationsJson} entryData={entryData}/> */}
        <Popup nested trigger={<Button>Open annotation JSON editor</Button>} modal>
            {/* <span> Modal content </span> */}
            {close => (
                <>
                    <button className="close" onClick={close}>
                        &times;
                    </button>
                    <JSONEditorComponent jsonData={annotationsJson} entryData={entryData} />
                </>

            )}
        </Popup>
        {/* Fitted models */}
        {allPdbs && allPdbs.length > 0 && <ExpandGroup header='Fitted models in PDB' initiallyExpanded>
            {allPdbs.map(pdb =>
                <WaitingButton key={pdb} onClick={() => entryData.actionShowFittedModel(visibleModels.includes(pdb) ? [] : [pdb])}
                    style={{ fontWeight: visibleModels.includes(pdb) ? 'bold' : undefined, textAlign: 'left', marginTop: 1 }}>
                    {pdb}
                </WaitingButton>
            )}
        </ExpandGroup>}

        {/* Volume */}
        <VolumeControls entryData={entryData} />
        <SegmentationControls model={entryData} />

        {/* Descriptions */}
        <SelectedSegmentDescription model={entryData} targetSegmentationId={segmentationId} targetKind={kind}></SelectedSegmentDescription>

        {<ExpandGroup header='Edit descriptions' initiallyExpanded>
            <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file' style={{ marginTop: '1px' }}>
                {'Load JSON with descriptions'} <input onChange={async v => {
                    const data = await parseJSONwithAnnotationsOrDescriptions(v, entryData);
                    await entryData.editDescriptions(data as DescriptionData[]);
                }} type='file' multiple={false} />
            </div>
        </ExpandGroup>}
        {<ExpandGroup header='Edit segment annotations' initiallyExpanded>
            <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file' style={{ marginTop: '1px' }}>
                {'Load JSON with segment annotations'} <input onChange={async v => {
                    const data = await parseJSONwithAnnotationsOrDescriptions(v, entryData);
                    await entryData.editSegmentAnnotations(data as SegmentAnnotationData[]);
                }} type='file' multiple={false} />
            </div>
        </ExpandGroup>}
    </>;
}

function TimeFrameSlider({ entryData }: { entryData: VolsegEntryData }) {
    // gets time info from volume
    const timeInfo = entryData.metadata.value!.raw.grid.volumes.time_info;
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


function _getVisualTransformFromProjectDataTransform(model: VolsegEntryData, projectDataTransform) {
    // TODO: if geometric segmentation - need child ref of child ref 
    // conform
    // if () {
    //     return 
    // } else {
    const childRef = model.plugin.state.data.tree.children.get(projectDataTransform.ref).toArray()[0];
    const segmentationRepresentation3DNode = findNodesByRef(model.plugin, childRef);
    // in case of CVSX segmentation.transform is already 3D representation
    // how to get segmentation Id from it?
    const transform = segmentationRepresentation3DNode.transform;
    if (transform.params.descriptions) {
        debugger;
        const childChildRef = model.plugin.state.data.tree.children.get(segmentationRepresentation3DNode.transform.ref).toArray()[0];
        const t = findNodesByRef(model.plugin, childChildRef);
        // TODO: get childRef of childRef
        debugger;
        return t.transform;
    } else {
        return transform;
    }
    // }

}
// TODO: TODO: TODO: exclude Opacity from state
function SegmentationSetControls({ model, segmentation, kind }: { model: VolsegEntryData, segmentation: StateObjectCell<PluginStateObject.Volume.Data> | StateObjectCell<VolsegGeometricSegmentation> | StateObjectCell<VolsegMeshSegmentation>, kind: 'lattice' | 'mesh' | 'primitive' }) {
    const projectDataTransform = segmentation.transform;
    if (!projectDataTransform) return null;
    const params: ProjectLatticeSegmentationDataParamsValues | ProjectGeometricSegmentationDataParamsValues | ProjectMeshSegmentationDataParamsValues = projectDataTransform.params;

    const segmentationId = params.segmentationId;


    // TODO: if geometric segmentation - need child ref of child ref 
    // const childRef = model.plugin.state.data.tree.children.get(projectDataTransform.ref).toArray()[0];
    // const segmentationRepresentation3DNode = findNodesByRef(model.plugin, childRef);
    // // in case of CVSX segmentation.transform is already 3D representation
    // // how to get segmentation Id from it?
    // const transform = segmentationRepresentation3DNode.transform;
    const transform = _getVisualTransformFromProjectDataTransform(model, projectDataTransform);
    if (!transform) return null;

    let opacity = undefined;
    if (transform.params?.type) {
        opacity = transform.params?.type.params.alpha;
    } else {
        // TODO: fix
        opacity = transform.params.alpha;
    }
    debugger;
    return <ExpandGroup header={`${segmentationId}`}>
        {/* TODO: use actual opacity */}
        {/* <div>Segmentation</div> */}
        <ControlRow label='Opacity' control={
            // TODO: problem is here
            <WaitingSlider min={0} max={1} value={opacity} step={0.05} onChange={async v => await model.actionSetOpacity(v, segmentationId, kind)} />
        } />
        <DescriptionsList
            model={model} targetSegmentationId={segmentationId} targetKind={kind}
        ></DescriptionsList>

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

export function SegmentationControls({ model }: { model: VolsegEntryData | CVSXStateModel }) {
    const h = useBehavior(model.state.hierarchy);
    if (!h) return null;
    return <>
        {/* <Button onClick={() => { console.log('volume cache, segmentation cache: ', entryData.cachedVolumeTimeframesData, entryData.cachedSegmentationTimeframesData); }}>Get volume and segmentation cache</Button> */}
        <ExpandGroup header='Segmentation data'>
            {/* TODO: just lattices, need geometric and mesh segmentations as well */}
            {h.segmentations.map((v) => {
                return <SegmentationSetControls key={v.transform.ref} model={model} segmentation={v} kind={'lattice'} />;
            })}
            {h.meshSegmentations.map((v) => {
                return <SegmentationSetControls key={v.transform.ref} model={model} segmentation={v} kind={'mesh'} />;
            })}
            {h.geometricSegmentations.map((v) => {
                return <SegmentationSetControls key={v.transform.ref} model={model} segmentation={v} kind={'primitive'} />;

            })}
        </ExpandGroup>
    </>;
}


type ComponentParams<T extends React.Component<any, any, any> | ((props: any) => JSX.Element)> =
    T extends React.Component<infer P, any, any> ? P : T extends (props: infer P) => JSX.Element ? P : never;

export function WaitingSlider({ value, onChange, ...etc }: { value: number, onChange: (value: number) => any } & ComponentParams<Slider>) {
    const [changing, sliderValue, execute] = useAsyncChange(value);

    return <Slider value={sliderValue} disabled={changing} onChange={newValue => execute(onChange, newValue)} {...etc} />;
}

export function WaitingButton({ onClick, ...etc }: { onClick: () => any } & ComponentParams<typeof Button>) {
    const [changing, _, execute] = useAsyncChange(undefined);

    return <Button disabled={changing} onClick={() => execute(onClick, undefined)} {...etc}>
        {etc.children}
    </Button>;
}

export function WaitingParameterControls<T extends PD.Params>({ values, onChangeValues, ...etc }: { values: PD.ValuesFor<T>, onChangeValues: (values: PD.ValuesFor<T>) => any } & ComponentParams<ParameterControls<T>>) {
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
