/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { GetAppSvg, CubeSendSvg } from '../../mol-plugin-ui/controls/icons';
import { download } from '../../mol-util/download';
import { GeometryControls } from './controls';

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
            header: 'Export Geometries',
            isCollapsed: true,
            brand: { accent: 'cyan', svg: CubeSendSvg }
        };
    }

    protected renderControls(): JSX.Element {
        return <>
            <Button icon={GetAppSvg}
                onClick={this.saveObj} style={{ marginTop: 1 }}
                disabled={this.state.busy || !this.plugin.canvas3d?.reprCount.value}>
                Save OBJ + MTL
            </Button>
        </>;
    }

    componentDidMount() {
        this.subscribe(this.plugin.canvas3d!.reprCount, () => {
            if (!this.state.isCollapsed) this.forceUpdate();
        });
    }

    componentWillUnmount() {
        this._controls?.dispose();
        this._controls = void 0;
    }

    saveObj = async () => {
        try {
            this.setState({ busy: true });
            const data = await this.controls.exportObj();
            this.setState({ busy: false });

            download(new Blob([data.zipData]), data.filename);
        } catch {
            this.setState({ busy: false });
        }
    }
}