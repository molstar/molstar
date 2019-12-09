/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParameterControls } from '../controls/parameters';
import { PluginUIComponent } from '../base';
import { Icon } from '../controls/common';
import { debounceTime } from 'rxjs/operators';
import { Subject } from 'rxjs';
import { ViewportScreenshotHelper } from '../../mol-plugin/util/viewport-screenshot';

interface ImageControlsState {
    showPreview: boolean

    resolution?: ViewportScreenshotHelper.ResolutionSettings,
    isDisabled: boolean
}

export class DownloadScreenshotControls extends PluginUIComponent<{ close: () => void }, ImageControlsState> {
    state: ImageControlsState = {
        showPreview: true,
        resolution: this.plugin.helpers.viewportScreenshot?.currentResolution,
        isDisabled: false
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
            this.preview()
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
            })
            this.updateQueue.next();
        })

        this.subscribe(debounceTime(250)(this.plugin.canvas3d.didDraw), () => {
            if (this.state.isDisabled) return;
            this.updateQueue.next();
        })

        this.subscribe(this.plugin.state.dataState.events.isUpdating, v => {
            this.setState({ isDisabled: v })
            if (!v) this.updateQueue.next();
        })

        this.handlePreview();
    }

    private setProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'resolution') {
            const resolution = p.value as ViewportScreenshotHelper.ResolutionSettings
            if (resolution.name === 'custom') {
                this.plugin.helpers.viewportScreenshot!.currentResolution.type = 'custom';
                this.plugin.helpers.viewportScreenshot!.currentResolution.width = resolution.params.width;
                this.plugin.helpers.viewportScreenshot!.currentResolution.height = resolution.params.height;
            } else {
                this.plugin.helpers.viewportScreenshot!.currentResolution.type = resolution.name;
            }
            this.setState({ resolution });
        }
    }

    render() {
        return <div>
            <div className='msp-image-preview'>
                <img ref={this.imgRef} /><br />
                <span>Right-click the image to Copy.</span>
            </div>
            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.download} disabled={this.state.isDisabled}><Icon name='download' /> Download</button>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.openTab} disabled={this.state.isDisabled}><Icon name='export' /> Open in new Tab</button>
            </div>
            <ParameterControls params={this.plugin.helpers.viewportScreenshot!.params} values={this.plugin.helpers.viewportScreenshot!.values} onChange={this.setProps} isDisabled={this.state.isDisabled} />
        </div>
    }
}