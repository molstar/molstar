/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { camelCaseToWords, stringToWords } from '../../mol-util/string';
import * as React from 'react';
import { _Props, _State } from '../base';
import { ParamProps } from './parameters';
import { TextInput, Button } from './common';
import { DefaultColorSwatch } from '../../mol-util/color/swatches';

export class CombinedColorControl extends React.PureComponent<ParamProps<PD.Color>, { isExpanded: boolean }> {
    state = { isExpanded: !!this.props.param.isExpanded }

    protected update(value: Color) {
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    onClickSwatch = (e: React.MouseEvent<HTMLButtonElement>) => {
        const value = Color(+(e.currentTarget.getAttribute('data-color') || '0'));
        if (value !== this.props.value) {
            if (!this.props.param.isExpanded) this.setState({ isExpanded: false });
            this.update(value);
        }
    }

    onChangeText = (value: Color) => {
        if (value !== this.props.value) {
            this.update(value);
        }
    }

    swatch() {
        // const def = this.props.param.defaultValue;
        return <div className='msp-combined-color-swatch'>
            {DefaultColorSwatch.map(c => <Button key={c[1]} inline data-color={c[1]} onClick={this.onClickSwatch} style={{ background: Color.toStyle(c[1]) }} />)}
        </div>;
    }

    render() {
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <>
            <div className='msp-control-row'>
                <span title={this.props.param.description}>{label}</span>
                <div>
                    <Button onClick={this.toggleExpanded} inline className='msp-combined-color-button' style={{ background: Color.toStyle(this.props.value) }} />
                </div>
            </div>
            {this.state.isExpanded && <div className='msp-control-offset'>
                {this.swatch()}
                <div className='msp-control-row'>
                    <span>RGB</span>
                    <div>
                        <TextInput onChange={this.onChangeText} value={this.props.value}
                            fromValue={formatColorRGB} toValue={getColorFromString} isValid={isValidColorString}
                            className='msp-form-control' onEnter={this.props.onEnter} blurOnEnter={true} blurOnEscape={true}
                            placeholder='e.g. 127 127 127' delayMs={250} />
                    </div>
                </div>
            </div>}
        </>;
    }
}

function formatColorRGB(c: Color) {
    const [r, g, b] = Color.toRgb(c);
    return `${r} ${g} ${b}`;
}

function getColorFromString(s: string) {
    const cs = s.split(/\s+/g);
    return Color.fromRgb(+cs[0], +cs[1], +cs[2]);
}

function isValidColorString(s: string) {
    const cs = s.split(/\s+/g);
    if (cs.length !== 3 && !(cs.length === 4 && cs[3] === '')) return false;
    for (const c of cs) {
        if (c === '') continue;
        const n = +c;
        if ('' + n !== c) return false;
        if (n < 0 || n > 255) return false;
    }
    return true;
}

let _colors: React.ReactFragment | undefined = void 0;
export function ColorOptions() {
    if (_colors) return _colors;
    _colors = <>{DefaultColorSwatch.map(v =>
        <option key={v[1]} value={v[1]} style={{ background: `${Color.toStyle(v[1])}` }} >
            {stringToWords(v[0])}
        </option>
    )}</>;
    return _colors;
}

const DefaultColorSwatchMap = (function () {
    const map = new Map<Color, string>()
    for (const v of DefaultColorSwatch) map.set(v[1], v[0])
    return map
})()
export function ColorValueOption(color: Color) {
    return !DefaultColorSwatchMap.has(color) ? <option key={Color.toHexString(color)} value={color} style={{ background: `${Color.toStyle(color)}` }} >
        {Color.toRgbString(color)}
    </option> : null
}