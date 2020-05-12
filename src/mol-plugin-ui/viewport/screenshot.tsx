/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParameterControls } from '../controls/parameters';
import { PluginUIComponent } from '../base';
import { debounceTime } from 'rxjs/operators';
import { Subject } from 'rxjs';
import { ViewportScreenshotHelper } from '../../mol-plugin/util/viewport-screenshot';
import { Button, ExpandGroup } from '../controls/common';
import { CameraHelperProps } from '../../mol-canvas3d/helper/camera-helper';
import { PluginCommands } from '../../mol-plugin/commands';
import { StateExportImportControls, LocalStateSnapshotParams } from '../state/snapshots';
import { GetAppSvg, LaunchSvg } from '../controls/icons';

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
        isDisabled: false,
        resolution: this.plugin.helpers.viewportScreenshot?.currentResolution,
        transparent: this.plugin.helpers.viewportScreenshot?.transparent,
        axes: this.plugin.helpers.viewportScreenshot?.axes
    } as ImageControlsState

    private imgRef = React.createRef<HTMLImageElement>()
    private updateQueue = new Subject();

    private preview = async () => {
        if (!this.imgRef.current) return;
        this.imgRef.current!.src = await this.plugin.helpers.viewportScreenshot!.imageData();
    }

    private download = () => {
        this.plugin.helpers.viewportScreenshot?.download();
        this.props.close();
    }

    private openTab = () => {
        // modified from https://stackoverflow.com/questions/16245767/creating-a-blob-from-a-base64-string-in-javascript/16245768#16245768

        const base64 = this.imgRef.current!.src;
        const byteCharacters = atob(base64.substr(`data:image/png;base64,`.length));
        const byteArrays = [];

        const sliceSize = Math.min(byteCharacters.length, 1024 * 1024);
        for (let offset = 0; offset < byteCharacters.length; offset += sliceSize) {
            const byteNumbers = new Uint8Array(Math.min(sliceSize, byteCharacters.length - offset));
            for (let i = 0, _i = byteNumbers.length; i < _i; i++) {
                byteNumbers[i] = byteCharacters.charCodeAt(offset + i);
            }
            byteArrays.push(byteNumbers);
        }
        const blob = new Blob(byteArrays, { type: 'image/png' });
        const blobUrl = URL.createObjectURL(blob);

        window.open(blobUrl, '_blank');
        this.props.close();
    }

    private handlePreview() {
        if (this.state.showPreview) {
            this.preview();
        }
    }

    componentDidUpdate() {
        this.updateQueue.next();
    }

    componentDidMount() {
        if (!this.plugin.canvas3d) return;

        this.subscribe(debounceTime(250)(this.updateQueue), () => this.handlePreview());

        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => {
            this.plugin.helpers.viewportScreenshot!.imagePass.setProps({
                multiSample: { mode: 'on', sampleLevel: 2 },
                postprocessing: this.plugin.canvas3d?.props.postprocessing
            });
            this.updateQueue.next();
        });

        this.subscribe(debounceTime(250)(this.plugin.canvas3d.didDraw), () => {
            if (this.state.isDisabled) return;
            this.updateQueue.next();
        });

        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
            if (!v) this.updateQueue.next();
        });

        this.handlePreview();
    }

    private setProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'resolution') {
            this.plugin.helpers.viewportScreenshot!.currentResolution = p.value;
            this.setState({ resolution: p.value });
        } else if (p.name === 'transparent') {
            this.plugin.helpers.viewportScreenshot!.transparent = p.value;
            this.setState({ transparent: p.value });
        } else if (p.name === 'axes') {
            this.plugin.helpers.viewportScreenshot!.axes = p.value;
            this.setState({ axes: p.value });
        }
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

    render() {
        return <div>
            <div className='msp-image-preview'>
                <img ref={this.imgRef} /><br />
                <span>Right-click the image to Copy.</span>
            </div>
            <div className='msp-flex-row'>
                <Button icon={GetAppSvg} onClick={this.download} disabled={this.state.isDisabled}>Download</Button>
                <Button icon={LaunchSvg} onClick={this.openTab} disabled={this.state.isDisabled}>Open in new Tab</Button>
            </div>
            <ParameterControls params={this.plugin.helpers.viewportScreenshot!.params} values={this.plugin.helpers.viewportScreenshot!.values} onChange={this.setProps} isDisabled={this.state.isDisabled} />
            <ExpandGroup header='State'>
                <StateExportImportControls onAction={this.props.close} />
                <ExpandGroup header='Save Options' initiallyExpanded={false} noOffset>
                    <LocalStateSnapshotParams />
                </ExpandGroup>
            </ExpandGroup>
        </div>;
    }
}