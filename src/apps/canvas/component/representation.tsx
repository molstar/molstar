/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import Viewer from 'mol-view/viewer';
import { App } from '../app';
import { Params } from 'mol-view/parameter';
import { Representation } from 'mol-geo/representation';
import { ParametersComponent } from 'mol-app/component/parameters';
import { Progress } from 'mol-task';

export interface RepresentationComponentProps {
    app: App
    viewer: Viewer
    repr: Representation<Params>
}

export interface RepresentationComponentState {

}

export class RepresentationComponent extends React.Component<RepresentationComponentProps, RepresentationComponentState> {

    async onChange(k: string, v: any) {
        await this.props.app.runTask(this.props.repr.createOrUpdate({ [k]: v }).run(
            progress => console.log(Progress.format(progress))
        ), 'Representation Update')
        this.props.viewer.add(this.props.repr)
        this.props.viewer.requestDraw(true)
    }

    render() {
        const { repr } = this.props
        // const ct = ColorTheme(colorTheme)

        return <div>
            <div>
                <h4>{repr.label}</h4>
            </div>
            <div>
                <ParametersComponent
                    params={repr.params}
                    values={repr.props}
                    onChange={(k, v) => this.onChange(k as string, v)}
                />
            </div>
        </div>;
    }
}