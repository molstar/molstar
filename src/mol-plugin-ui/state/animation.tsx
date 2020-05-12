/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { Button } from '../controls/common';
import { PlayArrowSvg } from '../controls/icons';

export class AnimationControls extends PluginUIComponent<{ onStart?: () => void }> {
    componentDidMount() {
        this.subscribe(this.plugin.managers.animation.events.updated, () => this.forceUpdate());
    }

    updateParams: ParamOnChange = p => {
        this.plugin.managers.animation.updateParams({ [p.name]: p.value });
    }

    updateCurrentParams: ParamOnChange = p => {
        this.plugin.managers.animation.updateCurrentParams({ [p.name]: p.value });
    }

    startOrStop = () => {
        const anim = this.plugin.managers.animation;
        if (anim.state.animationState === 'playing') anim.stop();
        else {
            if (this.props.onStart) this.props.onStart();
            anim.start();
        }
    }

    render() {
        const anim = this.plugin.managers.animation;
        if (anim.isEmpty) return null;

        const isDisabled = anim.state.animationState === 'playing';
        const canApply = anim.current.anim.canApply?.(this.plugin);

        return <>
            <ParameterControls params={anim.getParams()} values={anim.state.params} onChange={this.updateParams} isDisabled={isDisabled} />
            <ParameterControls params={anim.current.params} values={anim.current.paramValues} onChange={this.updateCurrentParams} isDisabled={isDisabled} />

            <div className='msp-flex-row'>
                <Button icon={anim.state.animationState !== 'playing' ? void 0 : PlayArrowSvg} onClick={this.startOrStop} disabled={canApply !== void 0 && canApply.canApply === false}>
                    {anim.state.animationState === 'playing' ? 'Stop' : canApply === void 0 || canApply.canApply ? 'Start' : canApply.reason || 'Start'}
                </Button>
            </div>
        </>;
    }
}