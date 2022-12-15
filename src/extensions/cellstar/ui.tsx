import { useEffect, useState } from 'react';

import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button, ControlRow, IconButton } from '../../mol-plugin-ui/controls/common';
import * as Icons from '../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import { Slider } from '../../mol-plugin-ui/controls/slider';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';
import { PluginContext } from '../../mol-plugin/context';
import { shallowEqualArrays } from '../../mol-util';
import { ParamDefinition } from '../../mol-util/param-definition';

import { CellstarEntry } from './entry-root';
import { VolumeTypeChoice } from './entry-state';
import { isDefined } from './helpers';


interface CellstarUIData {
    availableNodes: CellstarEntry[],
    activeNode?: CellstarEntry,
}
namespace CellstarUIData {
    export function changeAvailableNodes(data: CellstarUIData, newNodes: CellstarEntry[]): CellstarUIData {
        const newActiveNode = newNodes.length > data.availableNodes.length ?
            newNodes[newNodes.length - 1]
            : newNodes.find(node => node.data.ref === data.activeNode?.data.ref) ?? newNodes[0];
        return { availableNodes: newNodes, activeNode: newActiveNode };
    }
    export function changeActiveNode(data: CellstarUIData, newActiveRef: string): CellstarUIData {
        const newActiveNode = data.availableNodes.find(node => node.data.ref === newActiveRef) ?? data.availableNodes[0];
        return { availableNodes: data.availableNodes, activeNode: newActiveNode };
    }
    export function equals(data1: CellstarUIData, data2: CellstarUIData) {
        return shallowEqualArrays(data1.availableNodes, data2.availableNodes) && data1.activeNode === data2.activeNode;
    }
}

export class CellstarUI extends CollapsableControls<{}, { data: CellstarUIData }> {
    protected defaultState(): CollapsableState & { data: CellstarUIData } {
        return {
            header: 'Volume & Segmentation',
            isCollapsed: true,
            brand: { accent: 'orange', svg: Icons.ExtensionSvg },
            data: {
                availableNodes: [],
                activeNode: undefined,
            }
        };
    }
    protected renderControls(): JSX.Element | null {
        return <CellstarControls plugin={this.plugin} data={this.state.data} setData={d => this.setState({ data: d })} />;
    }
    componentDidMount(): void {
        this.setState({ isHidden: true, isCollapsed: false });
        // this.subscribe(this.plugin.state.data.events.cell.stateUpdated, e => console.log('cell.stateUpdated', e.ref, e.cell.obj?.label));
        // this.subscribe(this.plugin.state.data.events.cell.created, e => console.log('cell.created', e.ref, e.cell.obj?.label));
        // this.subscribe(this.plugin.state.data.events.cell.removed, e => console.log('cell.removed', e.ref));
        this.subscribe(this.plugin.state.data.events.changed, e => {
            // console.log('cell.changed', e);
            const nodes = e.state.selectQ(q => q.ofType(CellstarEntry)).map(cell => cell?.obj).filter(isDefined);
            // console.log('current nodes', ...this.state.data.availableNodes);
            // console.log('new nodes', ...nodes);
            const isHidden = nodes.length === 0;
            const newData = CellstarUIData.changeAvailableNodes(this.state.data, nodes);
            if (!CellstarUIData.equals(this.state.data, newData) || this.state.isHidden !== isHidden) {
                this.setState({ data: newData, isHidden: isHidden });
            }
        });
    }
}


function CellstarControls({ plugin, data, setData }: { plugin: PluginContext, data: CellstarUIData, setData: (d: CellstarUIData) => void }) {
    const entryData = data.activeNode?.data;
    if (!entryData) {
        return <p>No data!</p>;
    }

    const params = {
        /** Reference to the active CellstarEntry node */
        entry: ParamDefinition.Select(data.activeNode!.data.ref, data.availableNodes.map(entry => [entry.data.ref, entry.data.entryId]))
    };
    const values: ParamDefinition.Values<typeof params> = {
        entry: data.activeNode!.data.ref,
    };

    const state = useBehavior(entryData.currentState);

    const allSegments = entryData.metadata.allSegments;
    const selectedSegment = entryData.metadata.getSegment(state.selectedSegment);
    const visibleSegments = state.visibleSegments.map(seg => seg.segmentId);
    const visibleModels = state.visibleModels.map(model => model.pdbId);

    const allPdbs = entryData.pdbs;

    const volumeParams = {
        volumeType: VolumeTypeChoice.PDSelect(),
    };
    const volumeValues: ParamDefinition.Values<typeof volumeParams> = {
        volumeType: state.volumeType,
    };

    return <>
        {/* Entry select */}
        <ParameterControls params={params} values={values} onChangeValues={next => setData(CellstarUIData.changeActiveNode(data, next.entry))} />

        {/* Title */}
        <div style={{ padding: 8, overflow: 'hidden' }}>
            <p style={{ fontWeight: 'bold' }}>{entryData.metadata.raw.annotation?.name ?? 'Unnamed Annotation'}</p>
        </div>

        {/* Fitted models */}
        {allPdbs.length > 0 && <>
            <p style={{ margin: 5 }}><b>Fitted models in PDB:</b></p>
            {allPdbs.map(pdb =>
                <Button key={pdb} onClick={() => entryData.actionShowFittedModel(visibleModels.includes(pdb) ? [] : [pdb])}
                    style={{ fontWeight: visibleModels.includes(pdb) ? 'bold' : undefined, textAlign: 'left' }}>
                    {pdb}
                </Button>
            )}
        </>}

        {/* Volume */}
        <ParameterControls params={volumeParams} values={volumeValues} onChangeValues={next => entryData.actionSetVolumeVisual(next.volumeType)} />

        {/* Segment opacity slider */}
        <ControlRow label='Opacity' control={
            <WaitingSlider plugin={plugin} min={0} max={1} value={state.opacity} step={0.05} onChange={async v => await entryData.actionSetOpacity(v)} />
        } />

        {/* Segment toggles */}
        {allSegments.length > 0 && <>
            <Button onClick={() => entryData.actionToggleAllSegments()}
                style={{ marginTop: 1 }}>
                Toggle All segments
            </Button>
            <div style={{ maxHeight: 300, overflow: 'hidden', overflowY: 'auto', marginBlock: 1 }}>
                {allSegments.map(segment =>
                    <div style={{ display: 'flex', marginBottom: 1 }} key={segment.id}
                        onMouseEnter={() => entryData.actionHighlightSegment(segment)}
                        onMouseLeave={() => entryData.actionHighlightSegment()}>
                        <Button onClick={() => entryData.actionSelectSegment(segment !== selectedSegment ? segment.id : undefined)}
                            style={{ fontWeight: segment.id === selectedSegment?.id ? 'bold' : undefined, marginRight: 1, flexGrow: 1, textAlign: 'left' }}>
                            <div title={segment.biological_annotation.name ?? 'Unnamed segment'} style={{ maxWidth: 240, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
                                {segment.biological_annotation.name ?? 'Unnamed segment'}
                            </div>
                        </Button>
                        <IconButton svg={visibleSegments.includes(segment.id) ? Icons.VisibilityOutlinedSvg : Icons.VisibilityOffOutlinedSvg}
                            onClick={() => entryData.actionToggleSegment(segment.id)} />
                    </div>
                )}
            </div>
        </>}

        {/* Segment annotations */}
        <div style={{ padding: 8, maxHeight: 300, overflow: 'hidden', overflowY: 'auto' }}>
            {!selectedSegment && 'No segment selected'}
            {selectedSegment && <b>Segment {selectedSegment.id}:<br />{selectedSegment.biological_annotation.name ?? 'Unnamed segment'}</b>}
            {selectedSegment?.biological_annotation.external_references.map(ref =>
                <p key={ref.id} title={ref.description} style={{ marginTop: 4 }}>
                    {/* <b>{ref.resource}:{ref.accession}</b><br />
                    <i>{capitalize(ref.label)}:</i> {ref.description} */}
                    <small>{ref.resource}:{ref.accession}</small><br />
                    {capitalize(ref.label)}
                </p>)}
        </div>

    </>;
}

function WaitingSlider({ plugin, value, min, max, step, onChange }: { plugin: PluginContext, value: number, min: number, max: number, step: number, onChange: (value: number) => any }) {
    const [sliderValue, setSliderValue] = useState(value);
    const [changing, setChanging] = useState(false);
    useEffect(() => setSliderValue(value), [value]);

    return <Slider min={min} max={max} step={step} value={sliderValue} disabled={changing} onChange={async newValue => {
        setChanging(true);
        setSliderValue(newValue);
        try {
            await onChange(newValue);
        } catch (err) {
            setSliderValue(value); // reset original value
            throw err;
        } finally {
            setChanging(false);
        }
    }} />;
}


function capitalize(text: string) {
    const first = text.charAt(0);
    const rest = text.slice(1);
    return first.toUpperCase() + rest;
}