import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { PluginContext } from '../../mol-plugin/context';
import { Button, ControlRow, IconButton } from '../../mol-plugin-ui/controls/common';
import { ParamDefinition } from '../../mol-util/param-definition';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';

import { CellStarEntry } from './entry-root';
import { isDefined } from './helpers';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import * as Icons from '../../mol-plugin-ui/controls/icons';
import { Slider } from '../../mol-plugin-ui/controls/slider';


interface CellStarUIData {
    availableNodes: CellStarEntry[],
    activeNode?: CellStarEntry,
}
namespace CellStarUIData {
    export function changeAvailableNodes(data: CellStarUIData, newNodes: CellStarEntry[]): CellStarUIData {
        if (newNodes.length === data.availableNodes.length) {
            // No change
            return data;
        } else if (newNodes.length > data.availableNodes.length) {
            // Added
            return { availableNodes: newNodes, activeNode: newNodes[newNodes.length - 1] };
        } else {
            // Removed
            const newActiveNode = newNodes.find(node => node.id === data.activeNode?.id) ?? newNodes[0];
            return { availableNodes: newNodes, activeNode: newActiveNode };
        }
    }
    export function changeActiveNode(data: CellStarUIData, newActive: string): CellStarUIData {
        const newActiveNode = data.availableNodes.find(node => node.id === newActive) ?? data.availableNodes[0];
        return { availableNodes: data.availableNodes, activeNode: newActiveNode };
    }
}

export class CellStarUI extends CollapsableControls<{}, { data: CellStarUIData }> {
    protected defaultState(): CollapsableState & { data: CellStarUIData } {
        return {
            header: 'Cell* VolumeServer',
            isCollapsed: true,
            brand: { accent: 'orange', svg: Icons.ExtensionSvg },
            data: {
                availableNodes: [],
                activeNode: undefined,
            }
        };
    }
    protected renderControls(): JSX.Element | null {
        return <CellStarControls plugin={this.plugin} data={this.state.data} setData={d => this.setState({ data: d })} />;
    }
    componentDidMount(): void {
        this.setState({ isHidden: true, isCollapsed: false });
        this.subscribe(this.plugin.state.data.events.changed, e => {
            const nodes = e.state.selectQ(q => q.ofType(CellStarEntry)).map(cell => cell?.obj).filter(isDefined);
            const isHidden = nodes.length === 0;
            this.setState({ isHidden: isHidden });
            const newData = CellStarUIData.changeAvailableNodes(this.state.data, nodes);
            this.setState({ data: newData });
        });
    }
}


function CellStarControls({ plugin, data, setData }: { plugin: PluginContext, data: CellStarUIData, setData: (d: CellStarUIData) => void }) {
    const entryData = data.activeNode?.data;
    if (!entryData) {
        return <p>No data!</p>;
    }

    const params = {
        /** Reference to the active CellStarEntry node */
        entry: ParamDefinition.Select(data.activeNode!.id.toString(), data.availableNodes.map(entry => [entry.id.toString(), entry.data.entryId]))
    };
    const values: ParamDefinition.ValuesFor<typeof params> = {
        entry: data.activeNode!.id.toString(),
    };

    const allSegments = entryData.metadata.annotation?.segment_list ?? [];
    const currentSegment = useBehavior(entryData.currentSegment);
    const visibleSegments = useBehavior(entryData.visibleSegments);
    const opacity = useBehavior(entryData.opacity);

    const allPdbs = entryData.pdbs;
    const currentPdb = useBehavior(entryData.modelData.currentPdb);

    return <>
        <ParameterControls params={params} values={values} onChangeValues={next => setData(CellStarUIData.changeActiveNode(data, next.entry))} />

        {allPdbs.length > 0 && <>
            <p style={{ margin: 5 }}><b>Fitted models in PDB:</b></p>
            {allPdbs.map(pdb =>
                <Button key={pdb} onClick={() => entryData.modelData.showPdb(pdb === currentPdb ? undefined : pdb)}
                    style={{ fontWeight: pdb === currentPdb ? 'bold' : undefined, textAlign: 'left' }}>
                    {pdb}
                </Button>
            )}
        </>}


        <div style={{ padding: 8, maxHeight: 200, overflow: 'hidden', overflowY: 'auto' }}>
            <p style={{ fontWeight: 'bold' }}>{entryData.metadata.annotation?.name ?? 'Unnamed Annotation'}</p>
            {!currentSegment && 'No segment selected'}
            {currentSegment && `${currentSegment.biological_annotation.name} (${currentSegment.id})`}
            {currentSegment?.biological_annotation.external_references.map(ref =>
                <p key={ref.id} style={{ marginTop: 4 }}>
                    <b>{ref.resource}:{ref.accession}</b><br />
                    <i>{capitalize(ref.label)}:</i> {ref.description}
                </p>)}
        </div>

        <ControlRow label='Opacity' control={<Slider min={0} max={1} value={opacity} step={0.05} onChange={v => entryData.opacity.next(v)} />} />

        {allSegments.length > 0 && <>
            <Button onClick={() => entryData.toggleAllSegments()}
                style={{ marginTop: 1 }}>
                Toggle All segments
            </Button>
            {allSegments.map(segment =>
                <div style={{ display: 'flex', marginTop: 1 }} key={segment.id}
                    onMouseEnter={() => entryData.highlightSegment(segment)}
                    onMouseLeave={() => entryData.highlightSegment()}>
                    <Button onClick={() => entryData.showAnnotation(segment)}
                        style={{ fontWeight: segment.id === currentSegment?.id ? 'bold' : undefined, marginRight: 1, flexGrow: 1, textAlign: 'left' }}>
                        <div title={segment.biological_annotation.name ?? 'Unnamed segment'} style={{ maxWidth: 240, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
                            {segment.biological_annotation.name ?? 'Unnamed segment'}
                        </div>
                    </Button>
                    <IconButton svg={visibleSegments.includes(segment) ? Icons.VisibilityOutlinedSvg : Icons.VisibilityOffOutlinedSvg} onClick={() => entryData.toggleSegment(segment)} />
                </div>
            )}
        </>}
    </>;
}

function capitalize(text: string) {
    const first = text.charAt(0);
    const rest = text.slice(1);
    return first.toUpperCase() + rest;
}