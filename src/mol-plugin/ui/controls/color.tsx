/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '../../../mol-util/color';
import { ColorNames, ColorNamesValueMap } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { camelCaseToWords } from '../../../mol-util/string';
import * as React from 'react';
import { _Props, _State } from '../base';
import { ParamProps } from './parameters';
import { TextInput } from './common';

export class CombinedColorControl extends React.PureComponent<ParamProps<PD.Color>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    protected update(value: Color) {
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    onChangeSelect = (e: React.ChangeEvent<HTMLSelectElement>) => {
        const value = Color(parseInt(e.target.value));
        if (value !== this.props.value) {
            this.update(value);
        }
    }

    onClickSwatch = (e: React.MouseEvent<HTMLButtonElement>) => {
        const value = Color(+(e.currentTarget.getAttribute('data-color') || '0'));
        if (value !== this.props.value) {
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
            {/* <button title='Default Color' key={def} className='msp-form-control msp-btn' data-color={def} onClick={this.onClickSwatch} style={{ background: Color.toStyle(def) }}></button> */}
            {SwatchColors.map(c => <button key={c} className='msp-form-control msp-btn' data-color={c} onClick={this.onClickSwatch} style={{ background: Color.toStyle(c) }}></button>)}
        </div>;
    }

    stripStyle(): React.CSSProperties {
        return {
            background: Color.toStyle(this.props.value),
            position: 'absolute',
            bottom: '0',
            height: '4px',
            right: '0',
            left: '0'
        };
    }

    render() {
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <>
            <div className='msp-control-row'>
                <span>{label}</span>
                <div>
                    <button onClick={this.toggleExpanded} className='msp-combined-color-button' style={{ background: Color.toStyle(this.props.value) }}></button>
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
                <div className='msp-control-row'>
                    <span>Color List</span>
                    <div>
                        <select value={this.props.value} onChange={this.onChangeSelect}>
                            {ColorValueOption(this.props.value)}
                            {ColorOptions()}
                        </select>
                        <div style={this.stripStyle()} />
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

// the 1st color is the default value.
const SwatchColors = [
    0x000000, 0x808080, 0xFFFFFF, 0xD33115, 0xE27300, 0xFCC400,
    0x68BC00, 0x16A5A5, 0x009CE0, 0x7B64FF, 0xFA28FF, 0x7D2187
].map(Color);

let _colors: React.ReactFragment | undefined = void 0;
function ColorOptions() {
    if (_colors) return _colors;
    _colors = <>{Object.keys(ColorNames).map(name =>
        <option key={name} value={(ColorNames as { [k: string]: Color })[name]} style={{ background: `${Color.toStyle((ColorNames as { [k: string]: Color })[name])}` }} >
            {name}
        </option>
    )}</>;
    return _colors;
}

function ColorValueOption(color: Color) {
    return !ColorNamesValueMap.has(color) ? <option key={Color.toHexString(color)} value={color} style={{ background: `${Color.toStyle(color)}` }} >
        {Color.toRgbString(color)}
    </option> : null
}