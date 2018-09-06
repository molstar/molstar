/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { StructureRepresentation, StructureProps } from 'mol-geo/representation/structure';
import Viewer from 'mol-view/viewer';
import { VisualQuality, VisualQualityNames } from 'mol-geo/representation/util';
import { ColorThemeProps, ColorThemeName, ColorThemeNames } from 'mol-view/theme/color';

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
        this.props.viewer.requestDraw()

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

        return <div>
            <div>
                <h4>{label}</h4>
            </div>
            <div>
                <div>
                    <span>Visible</span>
                    <button onClick={(e) => this.update({ visible: !visible }) }>
                        {visible ? 'Hide' : 'Show'}
                    </button>
                </div>
                <div>
                    <span>Quality</span>
                    <select value={quality} onChange={(e) => this.update({ quality: e.target.value as VisualQuality }) }>
                        {VisualQualityNames.map(name => <option key={name} value={name}>{name}</option>)}
                    </select>
                </div>
                <div>
                    <span>Color Theme</span>
                    <select value={colorTheme.name} onChange={(e) => this.update({ colorTheme: { name: e.target.value as ColorThemeName } }) }>
                        {ColorThemeNames.map(name => <option key={name} value={name}>{name}</option>)}
                    </select>
                </div>
            </div>
        </div>;
    }
}