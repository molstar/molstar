/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { merge } from 'rxjs';
import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { GetAppSvg, CubeScanSvg, CubeSendSvg } from '../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import { download } from '../../mol-util/download';
import { GeometryParams, GeometryControls } from './controls';

interface State {
    busy?: boolean
}

export class GeometryExporterUI extends CollapsableControls<{}, State> {
    private _controls: GeometryControls | undefined;
    private isARSupported: boolean | undefined;

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
        if (this.isARSupported === undefined) {
            this.isARSupported = !!document.createElement('a').relList?.supports?.('ar');
        }
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
            {this.isARSupported && ctrl.behaviors.params.value.format === 'usdz' &&
                <Button icon={CubeScanSvg}
                    onClick={this.viewInAR} style={{ marginTop: 1 }}
                    disabled={this.state.busy || !this.plugin.canvas3d?.reprCount.value}>
                    View in AR
                </Button>
            }
        </>;
    }

    componentDidMount() {
        if (!this.plugin.canvas3d) return;

        const merged = merge(
            this.controls.behaviors.params,
            this.plugin.canvas3d!.reprCount
        );

        this.subscribe(merged, () => {
            if (!this.state.isCollapsed) this.forceUpdate();
        });
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        this._controls?.dispose();
        this._controls = void 0;
    }

    save = async () => {
        try {
            this.setState({ busy: true });
            const data = await this.controls.exportGeometry();
            download(data.blob, data.filename);
        } catch (e) {
            console.error(e);
        } finally {
            this.setState({ busy: false });
        }
    };

    viewInAR = async () => {
        try {
            this.setState({ busy: true });
            const data = await this.controls.exportGeometry();
            const a = document.createElement('a');
            a.rel = 'ar';
            a.href = URL.createObjectURL(data.blob);
            // For in-place viewing of USDZ on iOS, the link must contain a single child that is either an img or picture.
            // https://webkit.org/blog/8421/viewing-augmented-reality-assets-in-safari-for-ios/
            a.appendChild(document.createElement('img'));
            setTimeout(() => URL.revokeObjectURL(a.href), 4E4); // 40s
            setTimeout(() => a.dispatchEvent(new MouseEvent('click')));
        } catch (e) {
            console.error(e);
        } finally {
            this.setState({ busy: false });
        }
    };
}