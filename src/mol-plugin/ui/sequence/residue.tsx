/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { PurePluginUIComponent } from '../base';
import { getButtons, getModifiers } from '../../../mol-util/input/input-observer';
import { Sequence } from './sequence';
import { Color } from '../../../mol-util/color';

export class Residue extends PurePluginUIComponent<{ seqIdx: number, label: string, parent: Sequence<any>, marker: number, color: Color }> {

    mouseEnter = (e: React.MouseEvent) => {
        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)
        this.props.parent.hover(this.props.seqIdx, buttons, modifiers);
    }

    mouseLeave = (e: React.MouseEvent) => {
        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)
        this.props.parent.hover(undefined, buttons, modifiers);
    }

    mouseDown = (e: React.MouseEvent) => {
        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)
        this.props.parent.click(this.props.seqIdx, buttons, modifiers);
        e.stopPropagation() // so that `parent.mouseDown` is not called
    }

    get backgroundColor() {
        // TODO make marker color configurable
        if (this.props.marker === 0) return ''
        if (this.props.marker % 2 === 0) return 'rgb(51, 255, 25)' // selected
        if (this.props.marker === undefined) console.error('unexpected marker value')
        return 'rgb(255, 102, 153)' // highlighted
    }

    get margin() {
        return this.props.label.length > 1
            ? (this.props.seqIdx === 0 ? `0px 2px 0px 0px` : `0px 2px 0px 2px`)
            : undefined
    }

    render() {
        return <span
            onMouseEnter={this.mouseEnter}
            onMouseLeave={this.mouseLeave}
            onMouseDown={this.mouseDown}
            style={{
                color: Color.toStyle(this.props.color),
                backgroundColor: this.backgroundColor,
                margin: this.margin
            }}>
            {this.props.label}
        </span>;
    }
}