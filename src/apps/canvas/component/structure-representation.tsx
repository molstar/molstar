/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureRepresentation, StructureProps } from 'mol-geo/representation/structure';
import Viewer from 'mol-view/viewer';
import { ColorThemeProps, ColorThemeName, ColorThemeNames, ColorTheme } from 'mol-view/theme/color';
import { Color } from 'mol-util/color';
import { Progress } from 'mol-task';
import { VisualQuality, VisualQualityNames } from 'mol-geo/geometry/geometry';

export interface StructureRepresentationComponentProps {
    viewer: Viewer
    representation: StructureRepresentation<StructureProps>
}

export interface StructureRepresentationComponentState {
    label: string
    visible: boolean
    alpha: number
    quality: VisualQuality
    colorTheme: ColorThemeProps

    flatShaded?: boolean
    resolutionFactor?: number
    radiusOffset?: number
    smoothness?: number
}

export class StructureRepresentationComponent extends React.Component<StructureRepresentationComponentProps, StructureRepresentationComponentState> {
    state = {
        label: this.props.representation.label,
        visible: this.props.representation.props.visible,
        alpha: this.props.representation.props.alpha,
        quality: this.props.representation.props.quality,
        colorTheme: this.props.representation.props.colorTheme,

        flatShaded: (this.props.representation.props as any).flatShaded,
        resolutionFactor: (this.props.representation.props as any).resolutionFactor,
        radiusOffset: (this.props.representation.props as any).radiusOffset,
        smoothness: (this.props.representation.props as any).smoothness,
    }

    componentWillMount() {
        const repr = this.props.representation

        this.setState({
            ...this.state,
            label: repr.label,
            visible: repr.props.visible,
            alpha: repr.props.alpha,
            quality: repr.props.quality,
            colorTheme: repr.props.colorTheme,

            flatShaded: (repr.props as any).flatShaded,
            resolutionFactor: (repr.props as any).resolutionFactor,
            radiusOffset: (repr.props as any).probeRadius,
            smoothness: (repr.props as any).smoothness,
        })
    }

    async update(state: Partial<StructureRepresentationComponentState>) {
        const repr = this.props.representation
        const props: Partial<StructureProps> = {}

        if (state.visible !== undefined) props.visible = state.visible
        if (state.quality !== undefined) props.quality = state.quality
        if (state.alpha !== undefined) props.alpha = state.alpha
        if (state.colorTheme !== undefined) props.colorTheme = state.colorTheme

        if (state.flatShaded !== undefined) (props as any).flatShaded = state.flatShaded
        if (state.resolutionFactor !== undefined) (props as any).resolutionFactor = state.resolutionFactor
        if (state.radiusOffset !== undefined) (props as any).probeRadius = state.radiusOffset
        if (state.smoothness !== undefined) (props as any).smoothness = state.smoothness

        await repr.createOrUpdate(props).run(
            progress => console.log(Progress.format(progress))
        )
        this.props.viewer.add(repr)
        this.props.viewer.draw(true)
        console.log(this.props.viewer.stats)

        console.log(
            'drawCount',
            repr.renderObjects[0].values.drawCount.ref.version,
            repr.renderObjects[0].values.drawCount.ref.value,
            'dColorType',
            repr.renderObjects[0].values.dColorType.ref.version,
            repr.renderObjects[0].values.dColorType.ref.value
        )

        const newState = {
            ...this.state,
            visible: repr.props.visible,
            quality: repr.props.quality,
            alpha: repr.props.alpha,
            colorTheme: repr.props.colorTheme,

            flatShaded: (repr.props as any).flatShaded,
            resolutionFactor: (repr.props as any).resolutionFactor,
            probeRadius: (repr.props as any).probeRadius,
            isoValue: (repr.props as any).isoValue,
        }
        this.setState(newState)
    }

    render() {
        const { label, visible, quality, alpha, colorTheme } = this.state

        const ct = ColorTheme(colorTheme)

        if (ct.legend && ct.legend.kind === 'scale-legend') {
            // console.log(`linear-gradient(to right, ${ct.legend.colors.map(c => Color.toStyle(c)).join(', ')})`)
        }

        return <div>
            <div>
                <h4>{label}</h4>
            </div>
            <div>
                <div>
                    <span>Visible </span>
                    <button onClick={(e) => this.update({ visible: !visible }) }>
                        {visible ? 'Hide' : 'Show'}
                    </button>
                </div>
                { this.state.flatShaded !== undefined ? <div>
                    <span>Flat Shaded </span>
                    <button onClick={(e) => this.update({ flatShaded: !this.state.flatShaded }) }>
                        {this.state.flatShaded ? 'Deactivate' : 'Activate'}
                    </button>
                </div> : '' }
                <div>
                    <span>Quality </span>
                    <select value={quality} onChange={(e) => this.update({ quality: e.target.value as VisualQuality }) }>
                        {VisualQualityNames.map(name => <option key={name} value={name}>{name}</option>)}
                    </select>
                </div>
                <div>
                    <span>Opacity </span>
                    <input type='range'
                        defaultValue={alpha.toString()}
                        min='0'
                        max='1'
                        step='0.05'
                        onInput={(e) => this.update({ alpha: parseFloat(e.currentTarget.value) })}
                    >
                    </input>
                </div>
                { this.state.resolutionFactor !== undefined ? <div>
                    <span>Resolution Factor </span>
                    <input type='range'
                        defaultValue={this.state.resolutionFactor.toString()}
                        min='4'
                        max='9'
                        step='1'
                        onInput={(e) => this.update({ resolutionFactor: parseInt(e.currentTarget.value) })}
                    >
                    </input>
                </div> : '' }
                { this.state.smoothness !== undefined ? <div>
                    <span>Smoothness </span>
                    <input type='range'
                        defaultValue={this.state.smoothness.toString()}
                        min='1'
                        max='3'
                        step='0.1'
                        onInput={(e) => this.update({ smoothness: parseFloat(e.currentTarget.value) })}
                    >
                    </input>
                </div> : '' }
                { this.state.radiusOffset !== undefined ? <div>
                    <span>Radius Offset </span>
                    <input type='range'
                        defaultValue={this.state.radiusOffset.toString()}
                        min='0'
                        max='10'
                        step='0.1'
                        onInput={(e) => this.update({ radiusOffset: parseFloat(e.currentTarget.value) })}
                    >
                    </input>
                </div> : '' }
                <div>
                    <span>Color Theme </span>
                    <select value={colorTheme.name} onChange={(e) => this.update({ colorTheme: { name: e.target.value as ColorThemeName } }) }>
                        {ColorThemeNames.map(name => <option key={name} value={name}>{name}</option>)}
                    </select>
                    {ct.description ? <div><i>{ct.description}</i></div> : ''}
                    {
                        ct.legend && ct.legend.kind === 'scale-legend'
                            ? <div
                                style={{
                                    width: '100%',
                                    height: '30px',
                                    background: `linear-gradient(to right, ${ct.legend.colors.map(c => Color.toStyle(c)).join(', ')})`
                                }}
                            >
                                <span style={{float: 'left', padding: '6px', color: 'white', fontWeight: 'bold', backgroundColor: 'rgba(0, 0, 0, 0.2)'}}>{ct.legend.min}</span>
                                <span style={{float: 'right', padding: '6px', color: 'white', fontWeight: 'bold', backgroundColor: 'rgba(0, 0, 0, 0.2)'}}>{ct.legend.max}</span>
                            </div>
                        : ct.legend && ct.legend.kind === 'table-legend'
                            ? <div>
                                {ct.legend.table.map((value, i) => {
                                    const [name, color] = value
                                    return <div key={i} style={{minWidth: '60px', marginRight: '5px', display: 'inline-block'}}>
                                        <div style={{width: '30px', height: '20px', backgroundColor: Color.toStyle(color), display: 'inline-block'}}></div>
                                        {name}
                                    </div>
                                })}
                            </div>
                        : ''
                    }
                </div>
            </div>
        </div>;
    }
}