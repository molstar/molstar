import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { PluginContext } from '../../mol-plugin/context';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';


export class CellStarUI extends CollapsableControls<{}, {}> {
    protected defaultState(): CollapsableState {
        return {
            header: 'CellStar',
            isCollapsed: true,
            // brand: { accent: 'cyan', svg: GetAppSvg }
        };
    }
    protected renderControls(): JSX.Element | null {
        return <CellStarControls plugin={this.plugin} />;
    }
}


function CellStarControls({ plugin }: { plugin: PluginContext }) {
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

    return <>
        <p>Blablabla</p>
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