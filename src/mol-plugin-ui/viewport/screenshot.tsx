/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { useRef, useState } from 'react';
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/operators';
import { Viewport } from '../../mol-canvas3d/camera/util';
import { CameraHelperProps } from '../../mol-canvas3d/helper/camera-helper';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { ViewportScreenshotHelper, ViewportScreenshotHelperParams } from '../../mol-plugin/util/viewport-screenshot';
import { PluginUIComponent } from '../base';
import { Button, ExpandGroup, ToggleButton } from '../controls/common';
import { CopySvg, CropFreeSvg, CropOrginalSvg, CropSvg, GetAppSvg } from '../controls/icons';
import { ParameterControls } from '../controls/parameters';
import { useBehavior } from '../hooks/use-behavior';
import { LocalStateSnapshotParams, StateExportImportControls } from '../state/snapshots';

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

        const target = this.canvasRef.current;
        const w = target.clientWidth;
        const h = target.clientHeight;
        target.width = w;
        target.height = h;

        ctx.clearRect(0, 0, w, h);
        const frame = getViewportFrame(width, height, w, h);
        if (this.plugin.helpers.viewportScreenshot?.values.transparent) {
            this.drawCheckerboard(ctx, frame);
        }
        ctx.drawImage(canvas, frame.x, frame.y, frame.width, frame.height);
    }

    private drawCheckerboard(ctx: CanvasRenderingContext2D, frame: Viewport) {
        // must be odd number!
        const s = 13;
        for (let i = 0; i < frame.width; i += s) {
            for (let j = 0; j < frame.height; j += s) {
                ctx.fillStyle = (i + j) % 2 ? '#ffffff' : '#bfbfbf';

                const x = frame.x + i, y = frame.y + j;
                const w = i + s > frame.width ? frame.width - i : s;
                const h = j + s > frame.height ? frame.height - j : s;
                ctx.fillRect(x, y, w, h);
            }
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

    onCanvasContextMenu = (e: React.MouseEvent) => {
        e.preventDefault();
        e.stopPropagation();
    };

    render() {
        const values = this.plugin.helpers.viewportScreenshot!.values;

        return <div>
            <div className='msp-image-preview'>
                <div style={{ position: 'relative', width: '100%', height: '100%' }}>
                    <canvas ref={this.canvasRef} onClick={this.onCanvasContextMenu} onContextMenu={this.onCanvasContextMenu} style={{ width: '100%', height: '180px' }}></canvas>
                    <ViewportFrame plugin={this.plugin} canvasRef={this.canvasRef} />
                </div>
                <CropControls plugin={this.plugin} />
            </div>
            <div className='msp-flex-row'>
                {/* TODO: figure out how to do copy/paste in Firefox */}
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

function isFullFrame(crop: Viewport) {
    return equalEps(crop.x, 0, 1e-5) && equalEps(crop.y, 0, 1e-5) && equalEps(crop.width, 1, 1e-5) && equalEps(crop.height, 1, 1e-5);
}

function CropControls({ plugin }: { plugin: PluginContext }) {
    const helper = plugin.helpers.viewportScreenshot;
    const params = useBehavior(helper?.behaviors.values);
    const crop = useBehavior(helper?.behaviors.relativeCrop);

    if (!params || !crop || !helper) return null;

    return <div style={{ width: '100%', height: '24px' }}>
        <ToggleButton icon={CropOrginalSvg} title='Auto-crop' inline isSelected={params.autoCrop}
            style={{ background: 'transparent', float: 'left', width: 'auto',  height: '24px', lineHeight: '24px' }}
            toggle={() => helper.behaviors.values.next({ ...params, autoCrop: !params.autoCrop })} label={'Auto-crop ' + (params.autoCrop ? 'On' : 'Off') } />

        {!params.autoCrop && <Button icon={CropSvg} title='Crop'
            style={{ background: 'transparent', float: 'right', height: '24px', lineHeight: '24px', width: '24px', padding: '0' }}
            onClick={() => helper.autocrop()} />}
        {!isFullFrame(crop) && <Button icon={CropFreeSvg} title='Reset Crop'
            style={{ background: 'transparent', float: 'right', height: '24px', lineHeight: '24px', width: '24px', padding: '0' }}
            onClick={() => helper.resetCrop()} />}
    </div>;
}

function ViewportFrame({ plugin, canvasRef }: { plugin: PluginContext, canvasRef: React.RefObject<HTMLCanvasElement> }) {
    const helper = plugin.helpers.viewportScreenshot;
    const params = useBehavior(helper?.behaviors.values);
    const crop = useBehavior(helper?.behaviors.relativeCrop);
    const cropFrameRef = useRef<Viewport>({ x: 0, y: 0, width: 0, height: 0 });
    useBehavior(params?.resolution.name === 'viewport' ? plugin.canvas3d?.resized : void 0);

    const [drag, setDrag] = React.useState<string>('');
    const [start, setStart] = useState([0, 0]);
    const [current, setCurrent] = useState([0, 0]);

    if (!helper || !crop || !canvasRef.current) return null;

    const { width, height } = helper.getSizeAndViewport();
    const canvas = canvasRef.current;

    const frame = getViewportFrame(width, height, canvas.clientWidth, canvas.clientHeight);

    const cropFrame: Viewport = {
        x: frame.x + Math.floor(frame.width * crop.x),
        y: frame.y + Math.floor(frame.height * crop.y),
        width: Math.ceil(frame.width * crop.width),
        height: Math.ceil(frame.height * crop.height)
    };

    const rectCrop = toRect(cropFrame);
    const rectFrame = toRect(frame);

    if (drag === 'move') {
        rectCrop.l += current[0] - start[0];
        rectCrop.r += current[0] - start[0];
        rectCrop.t += current[1] - start[1];
        rectCrop.b += current[1] - start[1];
    } else if (drag) {
        if (drag.indexOf('left') >= 0) {
            rectCrop.l += current[0] - start[0];
        } else if (drag.indexOf('right') >= 0) {
            rectCrop.r += current[0] - start[0];
        }

        if (drag.indexOf('top') >= 0) {
            rectCrop.t += current[1] - start[1];
        } else if (drag.indexOf('bottom') >= 0) {
            rectCrop.b += current[1] - start[1];
        }
    }

    if (rectCrop.l > rectCrop.r) {
        const t = rectCrop.l;
        rectCrop.l = rectCrop.r;
        rectCrop.r = t;
    }

    if (rectCrop.t > rectCrop.b) {
        const t = rectCrop.t;
        rectCrop.t = rectCrop.b;
        rectCrop.b = t;
    }

    const pad = 40;
    rectCrop.l = Math.min(rectFrame.r - pad, Math.max(rectFrame.l, rectCrop.l));
    rectCrop.r = Math.max(rectFrame.l + pad, Math.min(rectFrame.r, rectCrop.r));
    rectCrop.t = Math.min(rectFrame.b - pad, Math.max(rectFrame.t, rectCrop.t));
    rectCrop.b = Math.max(rectFrame.t + pad, Math.min(rectFrame.b, rectCrop.b));

    cropFrame.x = rectCrop.l;
    cropFrame.y = rectCrop.t;
    cropFrame.width = rectCrop.r - rectCrop.l + 1;
    cropFrame.height = rectCrop.b - rectCrop.t + 1;

    cropFrameRef.current = cropFrame;

    const onMove = (e: MouseEvent) => {
        e.preventDefault();
        setCurrent([e.pageX, e.pageY]);
    };

    const onStart = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        setDrag(e.currentTarget.getAttribute('data-drag')! as any);
        const p = [e.pageX, e.pageY];
        setStart(p);
        setCurrent(p);
        window.addEventListener('mouseup', onEnd);
        window.addEventListener('mousemove', onMove);
    };

    const onEnd = () => {
        window.removeEventListener('mouseup', onEnd);
        window.removeEventListener('mousemove', onMove);

        const cropFrame = cropFrameRef.current;
        if (params?.autoCrop) {
            helper.behaviors.values.next({ ...params, autoCrop: false });
        }
        helper?.behaviors.relativeCrop.next({
            x: (cropFrame.x - frame.x) / frame.width,
            y: (cropFrame.y - frame.y) / frame.height,
            width: cropFrame.width / frame.width,
            height: cropFrame.height / frame.height
        });
        setDrag('');
        const p = [0, 0];
        setStart(p);
        setCurrent(p);
    };

    const contextMenu = (e: React.MouseEvent) => {
        e.preventDefault();
        e.stopPropagation();
    };

    const d = 4;
    const border = `3px solid rgba(255, 87, 45, 0.75)`;
    const transparent = 'transparent';

    return <>
        <div data-drag='move' style={{ position: 'absolute', left: cropFrame.x, top: cropFrame.y, width: cropFrame.width, height: cropFrame.height, border, cursor: 'move' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />

        <div data-drag='left' style={{ position: 'absolute', left: cropFrame.x - d, top: cropFrame.y + d, width: 4 * d, height: cropFrame.height - d, background: transparent, cursor: 'w-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />
        <div data-drag='right' style={{ position: 'absolute', left: rectCrop.r - 2 * d, top: cropFrame.y, width: 4 * d, height: cropFrame.height - d, background: transparent, cursor: 'w-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />
        <div data-drag='top' style={{ position: 'absolute', left: cropFrame.x - d, top: cropFrame.y - d, width: cropFrame.width + 2 * d, height: 4 * d, background: transparent, cursor: 'n-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />
        <div data-drag='bottom' style={{ position: 'absolute', left: cropFrame.x - d, top: rectCrop.b - 2 * d, width: cropFrame.width + 2 * d, height: 4 * d, background: transparent, cursor: 'n-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />

        <div data-drag='top, left' style={{ position: 'absolute', left: rectCrop.l - d, top: rectCrop.t - d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'nw-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />
        <div data-drag='bottom, right' style={{ position: 'absolute', left: rectCrop.r - 2 * d, top: rectCrop.b - 2 * d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'nw-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />
        <div data-drag='top, right' style={{ position: 'absolute', left: rectCrop.r - 2 * d, top: rectCrop.t - d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'ne-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />
        <div data-drag='bottom, left' style={{ position: 'absolute', left: rectCrop.l - d, top: rectCrop.b - 2 * d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'ne-resize' }} onMouseDown={onStart} draggable={false} onContextMenu={contextMenu} />
    </>;
}

function toRect(viewport: Viewport) {
    return { l: viewport.x, t: viewport.y, r: viewport.x + viewport.width - 1, b: viewport.y + viewport.height - 1 };
}

function getViewportFrame(srcWidth: number, srcHeight: number, w: number, h: number): Viewport {
    const a0 = srcWidth / srcHeight;
    const a1 = w / h;

    if (a0 <= a1) {
        const t = h * a0;
        return { x: Math.round((w - t) / 2), y: 0, width: Math.round(t), height: h };
    } else {
        const t = w / a0;
        return { x: 0, y: Math.round((h - t) / 2), width: w, height: Math.round(t) };
    }
}