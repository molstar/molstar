/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { merge } from 'rxjs';
import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { GetAppSvg, CubeSendSvg } from '../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import { download } from '../../mol-util/download';
import { GeometryParams, GeometryControls } from './controls';

interface State {
    busy?: boolean
}

export class GeometryExporterUI extends CollapsableControls<{}, State> {
    private _controls: GeometryControls | undefined;

    get controls() {
        return this._controls || (this._controls = new GeometryControls(this.plugin));
    }

    protected defaultState(): State & CollapsableState {
        return {
            header: 'Export Geometry',
            isCollapsed: true,
            brand: { accent: 'cyan', svg: CubeSendSvg }
        };
    }

    protected renderControls(): JSX.Element {
        const ctrl = this.controls;
        return <>
            <ParameterControls
                params={GeometryParams}
                values={ctrl.behaviors.params.value}
                onChangeValues={xs => ctrl.behaviors.params.next(xs)}
                isDisabled={this.state.busy}
            />
            <Button icon={GetAppSvg}
                onClick={this.save} style={{ marginTop: 1 }}
                disabled={this.state.busy || !this.plugin.canvas3d?.reprCount.value}>
                Save
            </Button>
        </>;
    }

    componentDidMount() {
        const merged = merge(
            this.controls.behaviors.params,
            this.plugin.canvas3d!.reprCount
        );

        this.subscribe(merged, () => {
            if (!this.state.isCollapsed) this.forceUpdate();
        });
    }

    componentWillUnmount() {
        this._controls?.dispose();
        this._controls = void 0;
    }

    save = async () => {
        try {
            this.setState({ busy: true });
            const data = await this.controls.exportGeometry();
            this.setState({ busy: false });

            download(data.blob, data.filename);
        } catch {
            this.setState({ busy: false });
        }
    }
}