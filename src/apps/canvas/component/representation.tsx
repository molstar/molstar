/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import Canvas3D from 'mol-canvas3d/canvas3d';
import { App } from '../app';
import { Params } from 'mol-util/parameter';
import { Representation } from 'mol-repr';
import { ParametersComponent } from 'mol-app/component/parameters';
import { ColorTheme } from 'mol-theme/color';
import { getColorThemeProps } from 'mol-geo/geometry/color-data';
import { ColorThemeComponent } from 'mol-app/component/color-theme';

export interface RepresentationComponentProps {
    app: App
    canvas3d: Canvas3D
    repr: Representation<Params>
}

export interface RepresentationComponentState {
    label: string
    reprParams: Params
    reprProps: Readonly<{}>
}

export class RepresentationComponent extends React.Component<RepresentationComponentProps, RepresentationComponentState> {

    private stateFromRepr(repr: Representation<Params>) {
        return {
            label: this.props.repr.label,
            reprParams: this.props.repr.params,
            reprProps: this.props.repr.props
        }
    }

    componentWillMount() {
        this.setState(this.stateFromRepr(this.props.repr))
    }

    async onChange(k: string, v: any) {
        const ctx = { webgl: this.props.canvas3d.webgl }
        await this.props.app.runTask(this.props.repr.createOrUpdate(ctx, { [k]: v }).run(
            progress => this.props.app.log(progress)
        ), 'Representation Update')
        this.props.canvas3d.add(this.props.repr)
        this.props.canvas3d.requestDraw(true)
        this.setState(this.stateFromRepr(this.props.repr))
    }

    render() {
        const { label, reprParams, reprProps } = this.state
        let colorTheme: ColorTheme | undefined = undefined
        if ('colorTheme' in reprProps) {
            colorTheme = ColorTheme(getColorThemeProps(reprProps))
        }

        return <div>
            <div>
                <h4>{label}</h4>
            </div>
            <div>
                <ParametersComponent
                    params={reprParams}
                    values={reprProps}
                    onChange={(k, v) => this.onChange(k as string, v)}
                />
            </div>
            { colorTheme !== undefined ? <ColorThemeComponent colorTheme={colorTheme} /> : '' }
        </div>;
    }
}