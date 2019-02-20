/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { ParameterControls, ParamOnChange } from '../controls/parameters';

export class AnimationControls extends PluginUIComponent<{ }> {
    componentDidMount() {
        this.subscribe(this.plugin.state.animation.events.updated, () => this.forceUpdate());
    }

    updateParams: ParamOnChange = p => {
        this.plugin.state.animation.updateParams({ [p.name]: p.value });
    }

    updateCurrentParams: ParamOnChange = p => {
        this.plugin.state.animation.updateCurrentParams({ [p.name]: p.value });
    }

    startOrStop = () => {
        const anim = this.plugin.state.animation;
        if (anim.state.animationState === 'playing') anim.stop();
        else anim.start();
    }

    render() {
        const anim = this.plugin.state.animation;
        if (anim.isEmpty) return null;

        const isDisabled = anim.state.animationState === 'playing';

        return <div className='msp-animation-section'>
            <div className='msp-section-header'>Animations</div>

            <ParameterControls params={anim.getParams()} values={anim.state.params} onChange={this.updateParams} isDisabled={isDisabled} />
            <ParameterControls params={anim.current.params} values={anim.current.paramValues} onChange={this.updateCurrentParams} isDisabled={isDisabled} />

            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.startOrStop}>
                    {anim.state.animationState === 'playing' ? 'Stop' : 'Start'}
                </button>
            </div>
        </div>
    }
}