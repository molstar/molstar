/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { ParameterControls } from '../controls/parameters';
import { PluginUIComponent } from '../base';
import { debounceTime } from 'rxjs/operators';
import { Subject } from 'rxjs';
import { ViewportScreenshotHelper, ViewportScreenshotHelperParams } from '../../mol-plugin/util/viewport-screenshot';
import { Button, ExpandGroup } from '../controls/common';
import { CameraHelperProps } from '../../mol-canvas3d/helper/camera-helper';
import { PluginCommands } from '../../mol-plugin/commands';
import { StateExportImportControls, LocalStateSnapshotParams } from '../state/snapshots';
import { CopySvg, GetAppSvg } from '../controls/icons';

interface ImageControlsState {
    showPreview: boolean
    isDisabled: boolean

    resolution?: ViewportScreenshotHelper.ResolutionSettings,
    transparent?: boolean,
    axes?: CameraHelperProps['axes']
}

export class DownloadScreenshotControls extends PluginUIComponent<{ close: () => void }, ImageControlsState> {
    state: ImageControlsState = {
        showPreview: true,
        isDisabled: false
    } as ImageControlsState

    private canvasRef = React.createRef<HTMLCanvasElement>()
    private updateQueue = new Subject();

    private preview = () => {
        if (!this.state.showPreview || !this.canvasRef.current || this.state.isDisabled) return;

        const { canvas, width, height } = this.plugin.helpers.viewportScreenshot?.getPreview()!;
        const ctx = this.canvasRef.current.getContext('2d');
        if (!ctx) return;

        const a0 = width / height;

        const target = this.canvasRef.current;
        const w = this.canvasRef.current.clientWidth;
        const h = this.canvasRef.current.clientHeight;
        target.width = w;
        target.height = h;
        const a1 = w / h;

        ctx.clearRect(0, 0, w, h);
        if (a0 <= a1) {
            const t = h * a0;
            ctx.drawImage(canvas, Math.round((w - t) / 2), 0, Math.round(t), h);
        } else {
            const t = w / a0;
            ctx.drawImage(canvas, 0, Math.round((h - t) / 2), w, Math.round(t));
        }
    }

    private download = () => {
        this.plugin.helpers.viewportScreenshot?.download();
        this.props.close();
    }

    private copy = async () => {
        await this.plugin.helpers.viewportScreenshot?.copyToClipboard();
        PluginCommands.Toast.Show(this.plugin, {
            message: 'Copied to clipboard.',
            title: 'Screenshot',
            timeoutMs: 1500
        });
    }

    componentDidMount() {
        if (!this.plugin.canvas3d) return;

        this.subscribe(this.updateQueue.pipe(debounceTime(33)), () => this.preview());

        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => {
            this.updateQueue.next();
        });

        this.subscribe(this.plugin.canvas3d.didDraw.pipe(debounceTime(150)), () => {
            if (this.state.isDisabled) return;
            this.updateQueue.next();
        });

        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
            if (!v) this.updateQueue.next();
        });

        this.subscribe(this.plugin.helpers.viewportScreenshot!.behaviors.values, () => {
            this.forceUpdate();
            this.updateQueue.next();
        });

        this.preview();
    }

    private setValues = (p: ViewportScreenshotHelperParams) => {
        this.plugin.helpers.viewportScreenshot!.behaviors.values.next(p);
    }

    downloadToFileJson = () => {
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { type: 'json' });
    }

    downloadToFileZip = () => {
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { type: 'zip' });
    }

    open = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (!e.target.files || !e.target.files![0]) return;
        PluginCommands.State.Snapshots.OpenFile(this.plugin, { file: e.target.files![0] });
    }

    onCanvasClick = (e: React.MouseEvent) => {
        e.preventDefault();
    };

    render() {
        const values = this.plugin.helpers.viewportScreenshot!.values;

        return <div>
            <div className='msp-image-preview'>
                <canvas ref={this.canvasRef} onClick={this.onCanvasClick} onContextMenu={this.onCanvasClick} style={{ width: '100%', height: '180px' }} className={values.transparent ? 'msp-transparent-screenshot' : void 0}></canvas><br />
            </div>
            <div className='msp-flex-row'>
                {!!(navigator.clipboard as any).write && <Button icon={CopySvg} onClick={this.copy} disabled={this.state.isDisabled}>Copy</Button>}
                <Button icon={GetAppSvg} onClick={this.download} disabled={this.state.isDisabled}>Download</Button>
            </div>
            <ParameterControls params={this.plugin.helpers.viewportScreenshot!.params} values={values} onChangeValues={this.setValues} isDisabled={this.state.isDisabled} />
            <ExpandGroup header='State'>
                <StateExportImportControls onAction={this.props.close} />
                <ExpandGroup header='Save Options' initiallyExpanded={false} noOffset>
                    <LocalStateSnapshotParams />
                </ExpandGroup>
            </ExpandGroup>
        </div>;
    }
}