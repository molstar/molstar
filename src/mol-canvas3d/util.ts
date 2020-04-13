/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Set canvas size taking `devicePixelRatio` into account */
export function setCanvasSize(canvas: HTMLCanvasElement, width: number, height: number) {
    canvas.width = Math.round(window.devicePixelRatio * width);
    canvas.height = Math.round(window.devicePixelRatio * height);
    Object.assign(canvas.style, { width: `${width}px`, height: `${height}px` });
}

/** Resize canvas to container element taking `devicePixelRatio` into account */
export function resizeCanvas (canvas: HTMLCanvasElement, container: Element) {
    let width = window.innerWidth;
    let height = window.innerHeight;
    if (container !== document.body) {
        let bounds = container.getBoundingClientRect();
        width = bounds.right - bounds.left;
        height = bounds.bottom - bounds.top;
    }
    setCanvasSize(canvas, width, height);
}

function _canvasToBlob(canvas: HTMLCanvasElement, callback: BlobCallback, type?: string, quality?: any) {
    const bin = atob(canvas.toDataURL(type, quality).split(',')[1]);
    const len = bin.length;
    const len32 = len >> 2;
    const a8 = new Uint8Array(len);
    const a32 = new Uint32Array(a8.buffer, 0, len32);

    let j = 0;
    for (let i = 0; i < len32; ++i) {
        a32[i] = bin.charCodeAt(j++) |
            bin.charCodeAt(j++) << 8 |
            bin.charCodeAt(j++) << 16 |
            bin.charCodeAt(j++) << 24;
    }

    let tailLength = len & 3;
    while (tailLength--) a8[j] = bin.charCodeAt(j++);

    callback(new Blob([a8], { type: type || 'image/png' }));
}

export async function canvasToBlob(canvas: HTMLCanvasElement, type?: string, quality?: any): Promise<Blob> {
    return new Promise((resolve, reject) => {
        const callback = (blob: Blob | null) => {
            if (blob) resolve(blob);
            else reject('no blob returned');
        };

        if (!HTMLCanvasElement.prototype.toBlob) {
            _canvasToBlob(canvas, callback, type, quality);
        } else {
            canvas.toBlob(callback, type, quality);
        }
    });
}