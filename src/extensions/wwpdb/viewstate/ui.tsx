/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useEffect, useState } from 'react';
import { PresetStructureRepresentations } from '../../../mol-plugin-state/builder/structure/representation-preset';
import { CollapsableControls, CollapsableState } from '../../../mol-plugin-ui/base';
import { Button } from '../../../mol-plugin-ui/controls/common';
import { BookmarksOutlinedSvg } from '../../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../../mol-plugin-ui/controls/parameters';
import { useBehavior } from '../../../mol-plugin-ui/hooks/use-behavior';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectCell } from '../../../mol-state';
import { download } from '../../../mol-util/download';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { wwPDBViewStateObject } from './model';

export class wwPDBViewStateEditUI extends CollapsableControls<{}, { cell?: StateObjectCell }> {
    protected defaultState(): CollapsableState {
        return {
            header: 'wwPDB View State',
            isCollapsed: true,
            isHidden: true,
            brand: { accent: 'cyan', svg: BookmarksOutlinedSvg }
        };
    }
    componentDidMount() {
        this.subscribe(this.plugin.state.data.events.changed, e => {
            if (e.inTransaction || this.plugin.behaviors.state.isAnimating.value) return;

            // TODO: support multiple cells and all that jazz
            const states = this.plugin.state.data.selectQ(q => q.ofType(wwPDBViewStateObject));
            if (states[0] !== this.state.cell) {
                this.setState(prev => ({
                    ...prev,
                    cell: states[0],
                    isHidden: states.length === 0
                }));
            }
        });
    }
    protected renderControls(): JSX.Element | null {
        if (!this.state.cell) return null;
        return <Controls plugin={this.plugin} cell={this.state.cell} />;
    }
}

const Params = {
    url: PD.Text('', { label: 'URL' }),
    presetName: PD.Select('atomic-detail', Object.keys(PresetStructureRepresentations).map(k => [k, k]))
};
const DefaultParams = PD.getDefaultValues(Params);
type Params = PD.ValuesFor<typeof Params>

function Controls({ plugin, cell }: { plugin: PluginContext, cell: StateObjectCell<wwPDBViewStateObject> }) {
    const [currentParams, setCurrentParams] = useState(DefaultParams);
    const isBusy = useBehavior(plugin.behaviors.state.isBusy);
    const hierarchy = plugin.managers.structure.hierarchy.current;

    const state = cell.obj?.data.state;

    useEffect(() => {
        if (!state) return;

        setCurrentParams({
            url: state.url,
            presetName: state.presetName as any,
        });
    }, [state]);

    if (!state) return null;

    const onExport = async () => {
        download(new Blob([JSON.stringify(state)]), `wwpdb-view-state-${Date.now()}.json`);
    };

    const applyParams = () => {
        plugin.build().to(cell).update({
            state: { ...state, ...currentParams }
        }).commit();
    };

    return <>
        <ParameterControls params={Params} values={currentParams} onChangeValues={setCurrentParams} isDisabled={isBusy} />
        <Button
            onClick={applyParams}
            style={{ marginTop: 1 }}
            commit={hierarchy.structures.length ? 'on' : 'off'}
        >
            Apply
        </Button>
        <Button
            onClick={onExport}
            style={{ marginTop: 1 }}
            commit={hierarchy.structures.length ? 'on' : 'off'}
        >
            Save
        </Button>
    </>;
}