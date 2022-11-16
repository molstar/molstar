import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { PluginContext } from '../../mol-plugin/context';
import { Button } from '../../mol-plugin-ui/controls/common';

import { CellStarEntry } from './entry-root';
import { isDefined } from './helpers';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';


interface CellStarUIData {
    entryNode?: CellStarEntry,
}

export class CellStarUI extends CollapsableControls<{}, { data: CellStarUIData }> {
    protected defaultState(): CollapsableState & { data: CellStarUIData } {
        return {
            header: 'CellStar',
            isCollapsed: true,
            // brand: { accent: 'cyan', svg: GetAppSvg }
            data: {
                entryNode: undefined
            }
        };
    }
    protected renderControls(): JSX.Element | null {
        return <CellStarControls plugin={this.plugin} data={this.state.data} />;
    }
    componentDidMount(): void {
        this.setState({ isHidden: true, isCollapsed: false });
        this.subscribe(this.plugin.state.data.events.changed, e => {
            const nodes = e.state.selectQ(q => q.ofType(CellStarEntry)).map(cell => cell?.obj).filter(isDefined);
            const isHidden = nodes.length === 0;
            this.setState({ isHidden: isHidden });
            this.setState({ data: { entryNode: nodes[0] } }); // TODO allow select entry if more entries
            console.log('event', e, nodes);
        });
    }
}


function CellStarControls({ plugin, data }: { plugin: PluginContext, data: CellStarUIData }) {
    // const [params, setParams] = useState(DefaultParams);
    // const [exporting, setExporting] = useState(false);
    // useBehavior(plugin.managers.structure.hierarchy.behaviors.selection); // triggers UI update
    // const isBusy = useBehavior(plugin.behaviors.state.isBusy);
    // const hierarchy = plugin.managers.structure.hierarchy.current;

    // let label: string = 'Nothing to Export';
    // if (hierarchy.structures.length === 1) {
    //     label = 'Export';
    // } if (hierarchy.structures.length > 1) {
    //     label = 'Export (as ZIP)';
    // }

    // const onExport = async () => {
    //     setExporting(true);
    //     try {
    //         await exportHierarchy(plugin, { format: params.format });
    //     } finally {
    //         setExporting(false);
    //     }
    // };
    const entryData = data.entryNode?.data;
    if (!entryData) {
        return <p>No data!</p>;
    }

    const allSegments = entryData.metadata.annotation.segment_list;
    const currentSegment = useBehavior(entryData.latticeSegmentationData.currentSegment);

    const allPdbs = entryData.pdbs;
    const currentPdb = useBehavior(entryData.modelData.currentPdb);

    return <>
        <p style={{margin: 5}}><b>Entry:</b> {entryData.entryId}</p>

        {allPdbs.length > 0 && <>
            <p style={{margin: 5}}><b>Fitted models in PDB:</b></p>
            {allPdbs.map(pdb =>
                <Button key={pdb} onClick={() => entryData.modelData.showPdb(pdb === currentPdb ? undefined : pdb)}
                    style={{ fontWeight: pdb === currentPdb ? 'bold' : undefined }}>
                    {pdb}
                </Button>
            )}
        </>}
        
        {allSegments.length > 0 && <>
            <p style={{margin: 5}}><b>Segmentation:</b></p>
            <Button onClick={() => entryData.latticeSegmentationData.showSegments(allSegments)}
                style={{ fontWeight: currentSegment?.id === undefined ? 'bold' : undefined }}>
                All segments
            </Button>
            {allSegments.map(segment =>
                <Button key={segment.id} onClick={() => entryData.latticeSegmentationData.showSegments([segment])}
                    style={{ fontWeight: segment.id === currentSegment?.id  ? 'bold' : undefined }}>
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
        
        {/* <ParameterControls params={Params} values={params} onChangeValues={setParams} isDisabled={isBusy || exporting} />
        <Button
            onClick={onExport}
            style={{ marginTop: 1 }}
            disabled={isBusy || hierarchy.structures.length === 0 || exporting}
            commit={hierarchy.structures.length ? 'on' : 'off'}
        >
            {label}
        </Button> */}
    </>;
}