import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { PluginContext } from '../../mol-plugin/context';
import { Button } from '../../mol-plugin-ui/controls/common';
import { ParamDefinition } from '../../mol-util/param-definition';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';

import { CellStarEntry } from './entry-root';
import { isDefined } from './helpers';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import * as Icons from '../../mol-plugin-ui/controls/icons';


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
            header: 'CellStar',
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
            const newData = CellStarUIData.changeAvailableNodes(this.state.data, nodes)
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
    }

    const allSegments = entryData.metadata.annotation.segment_list;
    const currentSegment = useBehavior(entryData.latticeSegmentationData.currentSegment);

    const allPdbs = entryData.pdbs;
    const currentPdb = useBehavior(entryData.modelData.currentPdb);

    return <>
        <ParameterControls params={params} values={values} onChangeValues={next => setData(CellStarUIData.changeActiveNode(data, next.entry))} />

        {allPdbs.length > 0 && <>
            <p style={{ margin: 5 }}><b>Fitted models in PDB:</b></p>
            {allPdbs.map(pdb =>
                <Button key={pdb} onClick={() => entryData.modelData.showPdb(pdb === currentPdb ? undefined : pdb)}
                    style={{ fontWeight: pdb === currentPdb ? 'bold' : undefined }}>
                    {pdb}
                </Button>
            )}
        </>}

        {allSegments.length > 0 && <>
            <p style={{ margin: 5 }}><b>Segmentation:</b></p>
            <Button onClick={() => entryData.showSegments(allSegments)}
                style={{ fontWeight: currentSegment?.id === undefined ? 'bold' : undefined }}>
                All segments
            </Button>
            {allSegments.map(segment =>
                <Button key={segment.id} onClick={() => entryData.showSegments([segment])}
                    style={{ fontWeight: segment.id === currentSegment?.id ? 'bold' : undefined }}>
                    {segment.biological_annotation.name ?? 'Unnamed segment'}
                </Button>
            )}
        </>}
        {currentSegment && <>
            {currentSegment.biological_annotation.external_references.map(ref => <p key={ref.id}>
                <b>{ref.resource}:{ref.accession}</b><br />
                {ref.description}
            </p>)}
        </>}
    </>;
}