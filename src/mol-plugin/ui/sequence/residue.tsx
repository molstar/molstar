/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { PurePluginUIComponent } from '../base';
import { getButtons, getModifiers } from '../../../mol-util/input/input-observer';
import { BaseSequence } from './base';

export class Residue extends PurePluginUIComponent<{ seqId: number, letter: string, parent: BaseSequence, marker: number }> {

    mouseEnter = (e: React.MouseEvent) => {
        const modifiers = getModifiers(e.nativeEvent)
        this.props.parent.highlight(this.props.seqId, modifiers);
    }

    mouseLeave = () => {
        this.props.parent.highlight();
    }

    mouseDown = (e: React.MouseEvent) => {
        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)
        this.props.parent.click(this.props.seqId, buttons, modifiers);
        e.stopPropagation() // so that `parent.mouseDown` is not called
    }

    getBackgroundColor() {
        // TODO make marker color configurable
        if (this.props.marker === 0) return ''
        if (this.props.marker % 2 === 0) return 'rgb(51, 255, 25)' // selected
        if (this.props.marker === undefined) console.error('unexpected marker value')
        return 'rgb(255, 102, 153)' // highlighted
    }

    render() {
        return <span
            onMouseEnter={this.mouseEnter}
            onMouseLeave={this.mouseLeave}
            onMouseDown={this.mouseDown}
            style={{ backgroundColor: this.getBackgroundColor() }}>
            {this.props.letter}
        </span>;
    }
}