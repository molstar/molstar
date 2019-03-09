/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { Icon } from '../controls/common';

export class AnimationControlsWrapper extends PluginUIComponent<{ }> {
    render() {
        const anim = this.plugin.state.animation;
        if (anim.isEmpty) return null;
        return <div className='msp-contols-section'>
            <div className='msp-section-header'>Animations</div>
            <AnimationControls />
        </div>
    }
}

export class AnimationControls extends PluginUIComponent<{ onStart?: () => void }> {
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
        else {
            if (this.props.onStart) this.props.onStart();
            anim.start();
        }
    }

    render() {
        const anim = this.plugin.state.animation;
        if (anim.isEmpty) return null;

        const isDisabled = anim.state.animationState === 'playing';

        return <>
            <ParameterControls params={anim.getParams()} values={anim.state.params} onChange={this.updateParams} isDisabled={isDisabled} />
            <ParameterControls params={anim.current.params} values={anim.current.paramValues} onChange={this.updateCurrentParams} isDisabled={isDisabled} />

            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.startOrStop}>
                    {anim.state.animationState !== 'playing' && <Icon name='play' />}
                    {anim.state.animationState === 'playing' ? 'Stop' : 'Start'}
                </button>
            </div>
        </>;
    }
}