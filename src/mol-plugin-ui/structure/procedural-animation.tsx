/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clearStructureWiggle, setStructureWiggleFromUncertainty } from '../../mol-plugin-state/helpers/structure-wiggle';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { Button } from '../controls/common';
import { AnimationSvg } from '../controls/icons';

export class StructureProceduralAnimationControls extends CollapsableControls {
    defaultState() {
        return {
            isCollapsed: false,
            header: 'Procedural Animation',
            brand: { accent: 'gray' as const, svg: AnimationSvg }
        };
    }

    renderControls() {
        return <StructureProceduralAnimation />;
    }
}

interface StructureProceduralAnimationState {
    busy: boolean;
}

class StructureProceduralAnimation extends PurePluginUIComponent<{}, StructureProceduralAnimationState> {
    state: StructureProceduralAnimationState = { busy: false };

    private get components() {
        return this.plugin.managers.structure.hierarchy.selection.structures.flatMap(s => s.components);
    }

    async applyUncertaintyWiggle() {
        this.setState({ busy: true });
        try {
            const options = this.plugin.managers.structure.component.state.options;
            await this.plugin.managers.structure.component.setOptions({
                ...options,
                animation: {
                    ...options.animation,
                    wiggleAmplitude: 0,
                    tumbleAmplitude: 0
                }
            });
            await setStructureWiggleFromUncertainty(this.plugin, this.components);
        } finally {
            this.setState({ busy: false });
        }
    }

    async applyDynamics() {
        this.setState({ busy: true });
        try {
            const options = this.plugin.managers.structure.component.state.options;
            await this.plugin.managers.structure.component.setOptions({
                ...options,
                animation: {
                    ...options.animation,
                    wiggleSpeed: 7,
                    wiggleAmplitude: 1,
                    wiggleFrequency: 0.2,
                }
            });
            await clearStructureWiggle(this.plugin, this.components);
        } finally {
            this.setState({ busy: false });
        }
    }

    async clearWiggle() {
        this.setState({ busy: true });
        try {
            const options = this.plugin.managers.structure.component.state.options;
            await this.plugin.managers.structure.component.setOptions({
                ...options,
                animation: {
                    ...options.animation,
                    wiggleAmplitude: 0,
                    tumbleAmplitude: 0
                }
            });
            await clearStructureWiggle(this.plugin, this.components);
        } finally {
            this.setState({ busy: false });
        }
    }

    render() {
        return <>
            <div className='msp-control-group-wrapper'>
                <div className='msp-control-group-header'><div><b>Apply Wiggle</b></div></div>
                <div className='msp-flex-row'>
                    <Button title='Set wiggle speed to 8 and amplitude to 1'
                        onClick={() => this.applyDynamics()} disabled={this.state.busy} >
                        Dynamics
                    </Button>
                    <Button title='Set per-group wiggle amplitude based on B-factor / RMSF uncertainty'
                        onClick={() => this.applyUncertaintyWiggle()} disabled={this.state.busy} >
                        Uncertainty
                    </Button>
                    <Button title='Set wiggle/tumble amplitude to zero and remove per-group wiggle layers'
                        onClick={() => this.clearWiggle()} disabled={this.state.busy} >
                        Clear
                    </Button>
                </div>
            </div>
        </>;
    }
}
