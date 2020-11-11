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
import { CopySvg, GetAppSvg, RefreshSvg } from '../controls/icons';
import { PluginContext } from '../../mol-plugin/context';
import { useBehavior } from '../hooks/use-behavior';
import { Viewport } from '../../mol-canvas3d/camera/util';
import { useEffect, useRef, useState } from 'react';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';

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
        ctx.drawImage(canvas, frame.x, frame.y, frame.width, frame.height);
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

        this.subscribe(this.plugin.helpers.viewportScreenshot!.behaviors.relativeCrop, () => {
            this.forceUpdate();
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
        const crop = this.plugin.helpers.viewportScreenshot!.relativeCrop;

        return <div>
            <div className='msp-image-preview'>
                <div style={{ position: 'relative', width: '100%', height: '100%' }}>
                    <canvas ref={this.canvasRef} onClick={this.onCanvasClick} onContextMenu={this.onCanvasClick} style={{ width: '100%', height: '180px' }} className={values.transparent ? 'msp-transparent-screenshot' : void 0}></canvas>
                    <ViewportFrame plugin={this.plugin} canvasRef={this.canvasRef} />
                </div>
                <span>
                    Drag the frame to crop the image.
                    {!isFullFrame(crop) && <Button icon={RefreshSvg} title='Reset Crop' inline
                        style={{ width: 'auto', background: 'transparent', height: '15px', lineHeight: '15px', padding: '0 4px', marginLeft: '8px', border: 'none' }}
                        onClick={() => this.plugin.helpers.viewportScreenshot!.behaviors.relativeCrop.next({ x: 0, y: 0, width: 1, height: 1 })} />}
                </span>
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

function isFullFrame(crop: Viewport) {
    return equalEps(crop.x, 0, 1e-5) && equalEps(crop.y, 0, 1e-5) && equalEps(crop.width, 1, 1e-5) && equalEps(crop.height, 1, 1e-5);
}

export function ViewportFrame({ plugin, canvasRef }: { plugin: PluginContext, canvasRef: React.RefObject<HTMLCanvasElement> }) {
    const helper = plugin.helpers.viewportScreenshot;
    const params = useBehavior(helper?.behaviors.values);
    const crop = useBehavior(helper?.behaviors.relativeCrop);
    const cropFrameRef = useRef<Viewport>({ x: 0, y: 0, width: 0, height: 0 });

    useBehavior(params?.resolution.name === 'viewport' ? plugin.canvas3d?.resized : void 0);
    useEffect(() => {
        return () => {
            if (onMove) document.removeEventListener('mousemove', onMove);
            if (onEnd) document.removeEventListener('mouseup', onEnd);
        };
    }, []);

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
    cropFrame.width = rectCrop.r - rectCrop.l;
    cropFrame.height = rectCrop.b - rectCrop.t;

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
        document.addEventListener('mouseup', onEnd);
        document.addEventListener('mousemove', onMove);
    };

    const onEnd = () => {
        document.removeEventListener('mouseup', onEnd);
        document.removeEventListener('mousemove', onMove);

        const cropFrame = cropFrameRef.current;
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

    const d = 4;
    const border = `6px solid rgba(255, 87, 45, 0.75)`;
    const transparent = 'transparent';

    return <>
        <div style={{ position: 'absolute', left: frame.x, top: frame.y, width: frame.width, height: frame.height, border: '1px solid rgba(0, 0, 0, 0.25)' }} />

        <div data-drag='move' style={{ position: 'absolute', left: cropFrame.x, top: cropFrame.y, width: cropFrame.width, height: cropFrame.height, border, cursor: 'move' }} onMouseDown={onStart} draggable={false} />

        <div data-drag='left' style={{ position: 'absolute', left: cropFrame.x - d, top: cropFrame.y + d, width: 4 * d, height: cropFrame.height - d, background: transparent, cursor: 'w-resize' }} onMouseDown={onStart} draggable={false} />
        <div data-drag='right' style={{ position: 'absolute', left: rectCrop.r - 2 * d, top: cropFrame.y, width: 4 * d, height: cropFrame.height - d, background: transparent, cursor: 'w-resize' }} onMouseDown={onStart} draggable={false} />
        <div data-drag='top' style={{ position: 'absolute', left: cropFrame.x - d, top: cropFrame.y - d, width: cropFrame.width + 2 * d, height: 4 * d, background: transparent, cursor: 'n-resize' }} onMouseDown={onStart} draggable={false} />
        <div data-drag='bottom' style={{ position: 'absolute', left: cropFrame.x - d, top: rectCrop.b - 2 * d, width: cropFrame.width + 2 * d, height: 4 * d, background: transparent, cursor: 'n-resize' }} onMouseDown={onStart} draggable={false} />

        <div data-drag='top, left' style={{ position: 'absolute', left: rectCrop.l - d, top: rectCrop.t - d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'nw-resize' }} onMouseDown={onStart} draggable={false} />
        <div data-drag='bottom, right' style={{ position: 'absolute', left: rectCrop.r - 2 * d, top: rectCrop.b - 2 * d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'nw-resize' }} onMouseDown={onStart} draggable={false} />
        <div data-drag='top, right' style={{ position: 'absolute', left: rectCrop.r - 2 * d, top: rectCrop.t - d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'ne-resize' }} onMouseDown={onStart} draggable={false} />
        <div data-drag='bottom, left' style={{ position: 'absolute', left: rectCrop.l - d, top: rectCrop.b - 2 * d, width: 4 * d, height: 4 * d, background: transparent, cursor: 'ne-resize' }} onMouseDown={onStart} draggable={false} />
    </>;
}

function toRect(viewport: Viewport) {
    return { l: viewport.x, t: viewport.y, r: viewport.x + viewport.width, b: viewport.y + viewport.height };
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