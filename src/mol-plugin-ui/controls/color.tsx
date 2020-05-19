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
import { TextInput, Button, ControlRow } from './common';
import { DefaultColorSwatch } from '../../mol-util/color/swatches';

export class CombinedColorControl extends React.PureComponent<ParamProps<PD.Color>, { isExpanded: boolean, lightness: number }> {
    state = {
        isExpanded: !!this.props.param.isExpanded,
        lightness: 0
    }

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

    onR = (v: number) => {
        const [, g, b] = Color.toRgb(this.props.value);
        const value = Color.fromRgb(v, g, b);
        if (value !== this.props.value) this.update(value);
    }

    onG = (v: number) => {
        const [r, , b] = Color.toRgb(this.props.value);
        const value = Color.fromRgb(r, v, b);
        if (value !== this.props.value) this.update(value);
    }

    onB = (v: number) => {
        const [r, g, ] = Color.toRgb(this.props.value);
        const value = Color.fromRgb(r, g, v);
        if (value !== this.props.value) this.update(value);
    }

    onLighten = () => {
        this.update(Color.lighten(this.props.value, 0.1));
    }

    onDarken = () => {
        this.update(Color.darken(this.props.value, 0.1));
    }

    swatch() {
        return <div className='msp-combined-color-swatch'>
            {DefaultColorSwatch.map(c => <Button key={c[1]} inline data-color={c[1]} onClick={this.onClickSwatch} style={{ background: Color.toStyle(c[1]) }} />)}
        </div>;
    }

    render() {
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const [r, g, b] = Color.toRgb(this.props.value);
        return <>
            <ControlRow title={this.props.param.description}
                label={label}
                control={<Button onClick={this.toggleExpanded} inline className='msp-combined-color-button' style={{ background: Color.toStyle(this.props.value) }} />} />
            {this.state.isExpanded && <div className='msp-control-offset'>
                {this.swatch()}
                <ControlRow label='RGB' className='msp-control-label-short' control={<div style={{ display: 'flex', textAlignLast: 'center', left: '80px' }}>
                    <TextInput onChange={this.onR} numeric value={r} delayMs={250} style={{ order: 1, flex: '1 1 auto', minWidth: 0 }} className='msp-form-control' onEnter={this.props.onEnter} blurOnEnter={true} blurOnEscape={true} />
                    <TextInput onChange={this.onG} numeric value={g} delayMs={250} style={{ order: 2, flex: '1 1 auto', minWidth: 0 }} className='msp-form-control' onEnter={this.props.onEnter} blurOnEnter={true} blurOnEscape={true} />
                    <TextInput onChange={this.onB} numeric value={b} delayMs={250} style={{ order: 3, flex: '1 1 auto', minWidth: 0 }} className='msp-form-control' onEnter={this.props.onEnter} blurOnEnter={true} blurOnEscape={true} />
                </div>}/>
                <div style={{ display: 'flex', textAlignLast: 'center' }}>
                    <Button  onClick={this.onLighten} style={{ order: 1, flex: '1 1 auto', minWidth: 0 }} className='msp-form-control'>Lighten</Button>
                    <Button onClick={this.onDarken} style={{ order: 1, flex: '1 1 auto', minWidth: 0 }} className='msp-form-control'>Darken</Button>
                </div>
            </div>}
        </>;
    }
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
    const map = new Map<Color, string>();
    for (const v of DefaultColorSwatch) map.set(v[1], v[0]);
    return map;
})();
export function ColorValueOption(color: Color) {
    return !DefaultColorSwatchMap.has(color) ? <option key={Color.toHexString(color)} value={color} style={{ background: `${Color.toStyle(color)}` }} >
        {Color.toRgbString(color)}
    </option> : null;
}