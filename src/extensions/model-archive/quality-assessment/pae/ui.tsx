/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CollapsableControls, CollapsableState } from '../../../../mol-plugin-ui/base';
import { Button } from '../../../../mol-plugin-ui/controls/common';
import { ExtensionSvg } from '../../../../mol-plugin-ui/controls/icons';
import { drawPAEPng } from './plot';

interface State {
    plot?: string;
}

export class PAEPlotUI extends CollapsableControls<{}, State> {
    protected defaultState(): State & CollapsableState {
        return {
            header: 'PAE',
            isCollapsed: false,
            brand: { accent: 'cyan', svg: ExtensionSvg }
        };
    }

    private draw = () => {
        const model = this.plugin.managers.structure.hierarchy.current.models[0];
        const plot = drawPAEPng(model.cell.obj?.data!);
        this.setState({ plot });
    };

    protected renderControls(): JSX.Element | null {
        const plot = this.state.plot;

        return <>
            <Button onClick={this.draw}>Draw</Button>
            <div style={{ margin: '8px 8px 0 8px', position: 'relative' }}>
                {plot && <img src={plot} style={{ width: '100%', height: 'auto', border: '1px solid black' }} />}
            </div>
        </>;
    }
}