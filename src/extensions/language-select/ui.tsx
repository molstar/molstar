/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

//import { DownloadFile } from '../../mol-plugin-state/actions/file';
//import { DownloadStructure, LoadTrajectory } from '../../mol-plugin-state/actions/structure';
//import { DownloadDensity } from '../../mol-plugin-state/actions/volume';
//import { CoordinatesFormatCategory } from '../../mol-plugin-state/formats/coordinates';
//import { TopologyFormatCategory } from '../../mol-plugin-state/formats/topology';
//import { TrajectoryFormatCategory } from '../../mol-plugin-state/formats/trajectory';
//import { VolumeFormatCategory } from '../../mol-plugin-state/formats/volume';
import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { SelectionModeSvg } from '../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
//import { PluginContext } from '../../mol-plugin/context';
//import { formatBytes } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';


type ScriptLanguage = undefined;

interface State {
    busy?: boolean
    languageValues: PD.Values<typeof ScriptImportParams>
    importParams?: ImportParams
    language?: ScriptLanguage
}

const ScriptImportParams = {
    language: PD.Text('', { description: 'Script Language' })
};

function createImportParams() {
    const params: PD.Params = {};
    let defaultType = '';
    return {
        type: PD.MappedStatic(defaultType, Object.keys(params).length ? params : { '': PD.EmptyGroup() })
    };
}
type ImportParams = ReturnType<typeof createImportParams>


export class ScriptImportUI extends CollapsableControls<{}, State> {
    protected defaultState(): State & CollapsableState {
        return {
            header: 'Scripting Language',
            isCollapsed: true,
            brand: { accent: 'cyan', svg: SelectionModeSvg },
            languageValues: PD.getDefaultValues(ScriptImportParams),
            importParams: undefined,
            language: undefined,
        };
    }
    private loadLanguage = async () => {
            this.plugin.log.message(`'${this.state.languageValues.language}'`);
    };

    private languageParamsOnChange = (values: any) => {
	this.setState({ languageValues: values });
    };

    private renderLanguageInfo(language: ScriptLanguage) {
	return <div style={{ marginBottom: 10 }}>
	    <div className='msp-help-text'>
	    <div>Language </div>
	    </div>
	    </div>;
    }
    
        
    private renderLoadLanguage() {
        return <div style={{ marginBottom: 10 }}>
            <ParameterControls params={ScriptImportParams} values={this.state.languageValues} onChangeValues={this.languageParamsOnChange} isDisabled={this.state.busy} />
            <Button onClick={this.loadLanguage} style={{ marginTop: 1 }} disabled={this.state.busy || !this.state.languageValues.language}>
                Set Language
            </Button>
        </div>;
    }

   
    protected renderControls(): JSX.Element | null {
	return <>
	{!this.state.language ? this.renderLoadLanguage() : null}
	{this.state.language ? this.renderLanguageInfo(this.state.language) : null}
	</>;
    }


}
