/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureRepresentation, StructureProps } from 'mol-geo/representation/structure';
import Viewer from 'mol-view/viewer';
import { VisualQuality, VisualQualityNames } from 'mol-geo/representation/util';
import { ColorThemeProps, ColorThemeName, ColorThemeNames, ColorTheme } from 'mol-view/theme/color';
import { Color } from 'mol-util/color';

export interface StructureRepresentationComponentProps {
    viewer: Viewer
    representation: StructureRepresentation<StructureProps>
}

export interface StructureRepresentationComponentState {
    label: string
    visible: boolean
    quality: VisualQuality
    colorTheme: ColorThemeProps
}

export class StructureRepresentationComponent extends React.Component<StructureRepresentationComponentProps, StructureRepresentationComponentState> {
    state = {
        label: this.props.representation.label,
        visible: this.props.representation.props.visible,
        quality: this.props.representation.props.quality,
        colorTheme: this.props.representation.props.colorTheme,
    }

    componentWillMount() {
        const repr = this.props.representation

        this.setState({
            ...this.state,
            label: repr.label,
            visible: repr.props.visible,
            quality: repr.props.quality,
            colorTheme: repr.props.colorTheme,
        })
    }

    async update(state: Partial<StructureRepresentationComponentState>) {
        const repr = this.props.representation
        const props: Partial<StructureProps> = {}

        if (state.visible !== undefined) props.visible = state.visible
        if (state.quality !== undefined) props.quality = state.quality
        if (state.colorTheme !== undefined) props.colorTheme = state.colorTheme

        await repr.createOrUpdate(props).run()
        this.props.viewer.add(repr)
        this.props.viewer.requestDraw(true)
        console.log(this.props.viewer.stats)

        const newState = {
            ...this.state,
            visible: repr.props.visible,
            quality: repr.props.quality,
            colorTheme: repr.props.colorTheme,
        }
        this.setState(newState)
    }

    render() {
        const { label, visible, quality, colorTheme } = this.state

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
                <div>
                    <span>Quality </span>
                    <select value={quality} onChange={(e) => this.update({ quality: e.target.value as VisualQuality }) }>
                        {VisualQualityNames.map(name => <option key={name} value={name}>{name}</option>)}
                    </select>
                </div>
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