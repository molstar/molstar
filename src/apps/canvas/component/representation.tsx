/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { App } from '../app';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Representation } from 'mol-repr/representation';
import { ParametersComponent } from 'mol-app/component/parameters';

export interface RepresentationComponentProps<P extends PD.Params> {
    app: App
    canvas3d: Canvas3D
    repr: Representation<P>
}

export interface RepresentationComponentState {
    label: string
    reprParams: PD.Params
    reprProps: Readonly<{}>
}

export class RepresentationComponent<P extends PD.Params> extends React.Component<RepresentationComponentProps<P>, RepresentationComponentState> {

    private stateFromRepr(repr: Representation<P>) {
        return {
            label: repr.label,
            reprParams: repr.params,
            reprProps: repr.props
        }
    }

    componentWillMount() {
        this.setState(this.stateFromRepr(this.props.repr))
    }

    async onChange(k: string, v: any) {
        await this.props.app.runTask(this.props.repr.createOrUpdate(this.props.app.reprCtx, { [k]: v }).run(
            progress => this.props.app.log(progress)
        ), 'Representation Update')
        this.setState(this.stateFromRepr(this.props.repr))
    }

    render() {
        const { label, reprParams, reprProps } = this.state
        // let colorTheme: ColorTheme | undefined = undefined
        // if ('colorTheme' in reprProps) {
        //     colorTheme = ColorTheme(getColorThemeProps(reprProps))
        // }

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
            {/* { colorTheme !== undefined ? <ColorThemeComponent colorTheme={colorTheme} /> : '' } */}
        </div>;
    }
}