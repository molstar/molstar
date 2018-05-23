/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'

import { View } from '../view';
import { Controller } from '../../controller/controller';
import { Toggle } from '../controls/common';
import { SpacefillEntity } from 'mol-view/state/entity';
import { SpacefillUpdate } from 'mol-view/state/transform'
import { StateContext } from 'mol-view/state/context';
import { ColorTheme } from 'mol-geo/theme';
import { Color, ColorNames } from 'mol-util/color';
import { Slider } from '../controls/slider';

export const ColorThemeInfo = {
    'atom-index': {},
    'chain-id': {},
    'element-symbol': {},
    'instance-index': {},
    'uniform': {}
}
export type ColorThemeInfo = keyof typeof ColorThemeInfo

interface SpacefillState {
    doubleSided: boolean
    detail: number
    colorTheme: ColorTheme
    colorValue: Color
    visible: boolean
    alpha: number
    depthMask: boolean
}

export class Spacefill extends View<Controller<any>, SpacefillState, { transform: SpacefillUpdate, entity: SpacefillEntity, ctx: StateContext }> {
    state = {
        doubleSided: true,
        detail: 2,
        colorTheme: { name: 'element-symbol' } as ColorTheme,
        colorValue: 0x000000,
        visible: true,
        alpha: 1,
        depthMask: true
    }

    update(state?: Partial<SpacefillState>) {
        const { transform, entity, ctx } = this.props
        console.log('update spacefill', transform, entity)
        const newState = { ...this.state, ...state }
        this.setState(newState)
        transform.apply(ctx, entity, newState)
    }

    render() {
        const { transform } = this.props

        const sphereDetailOptions = [0, 1, 2, 3].map((value, idx) => {
            return <option key={value} value={value}>{value.toString()}</option>
        })

        const colorThemeOptions = Object.keys(ColorThemeInfo).map((name, idx) => {
            return <option key={name} value={name}>{name}</option>
        })

        const colorValueOptions = Object.keys(ColorNames).map((name, idx) => {
            return <option key={name} value={(ColorNames as any)[name]}>{name}</option>
        })

        return <div className='molstar-transformer-wrapper'>
            <div className='molstar-panel molstar-control molstar-transformer molstar-panel-expanded'>
                <div className='molstar-panel-header'>
                    <button
                        className='molstar-btn molstar-btn-link molstar-panel-expander'
                        onClick={() => this.update()}
                    >
                        <span>[{transform.kind}] {transform.inputKind} -> {transform.outputKind}</span>
                    </button>
                </div>
                <div className='molstar-panel-body'>
                    <div>
                        <div className='molstar-control-row molstar-options-group'>
                            <span>Sphere detail</span>
                            <div>
                                <select
                                    className='molstar-form-control'
                                    value={this.state.detail}
                                    onChange={(e) => this.update({ detail: parseInt(e.target.value) })}
                                >
                                    {sphereDetailOptions}
                                </select>
                            </div>
                        </div>
                        <div className='molstar-control-row molstar-options-group'>
                            <span>Color theme</span>
                            <div>
                                <select
                                    className='molstar-form-control'
                                    value={this.state.colorTheme.name}
                                    onChange={(e) => {
                                        const colorThemeName = e.target.value as ColorThemeInfo
                                        if (colorThemeName === 'uniform') {
                                            this.update({
                                                colorTheme: {
                                                    name: colorThemeName,
                                                    value: this.state.colorValue
                                                }
                                            })
                                        } else {
                                            this.update({
                                                colorTheme: { name: colorThemeName }
                                            })
                                        }
                                    }}
                                >
                                    {colorThemeOptions}
                                </select>
                            </div>
                        </div>
                        <div className='molstar-control-row molstar-options-group'>
                            <span>Color value</span>
                            <div>
                                <select
                                    className='molstar-form-control'
                                    value={this.state.colorValue}
                                    onChange={(e) => {
                                        const colorValue = parseInt(e.target.value)
                                        this.update({
                                            colorTheme: {
                                                name: 'uniform',
                                                value: colorValue
                                            },
                                            colorValue
                                        })
                                    }}
                                >
                                    {colorValueOptions}
                                </select>
                            </div>
                        </div>
                        <div className='molstar-control-row molstar-options-group'>
                            <div>
                                <Toggle
                                    value={this.state.visible}
                                    label='Visibility'
                                    onChange={value => this.update({ visible: value })}
                                />
                            </div>
                        </div>
                        <div className='molstar-control-row molstar-options-group'>
                            <div>
                                <Toggle
                                    value={this.state.depthMask}
                                    label='Depth write'
                                    onChange={value => this.update({ depthMask: value })}
                                />
                            </div>
                        </div>
                        <div className='molstar-control-row molstar-options-group'>
                            <div>
                                <Slider
                                    value={this.state.alpha}
                                    label='Opacity'
                                    min={0}
                                    max={1}
                                    step={0.01}
                                    callOnChangeWhileSliding={true}
                                    onChange={value => this.update({ alpha: value })}
                                />
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>;
    }
}