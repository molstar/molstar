/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export class PluginContainer {
    readonly parent: HTMLDivElement;
    readonly canvas: HTMLCanvasElement;

    mount(target: HTMLElement) {
        if (this.parent.parentElement !== target) {
            this.parent.parentElement?.removeChild(this.parent);
        }
        target.appendChild(this.parent);
    }

    unmount() {
        this.parent.parentElement?.removeChild(this.parent);
    }

    constructor(public options?: { checkeredCanvasBackground?: boolean, canvas?: HTMLCanvasElement }) {
        const parent = document.createElement('div');
        Object.assign(parent.style, {
            position: 'absolute',
            left: 0,
            top: 0,
            right: 0,
            bottom: 0,
            '-webkit-user-select': 'none',
            'user-select': 'none',
            '-webkit-tap-highlight-color': 'rgba(0,0,0,0)',
            '-webkit-touch-callout': 'none',
            'touch-action': 'manipulation',
        });
        let canvas = options?.canvas;
        if (!canvas) {
            canvas = document.createElement('canvas');
            if (options?.checkeredCanvasBackground) {
                Object.assign(canvas.style, {
                    'background-image': 'linear-gradient(45deg, lightgrey 25%, transparent 25%, transparent 75%, lightgrey 75%, lightgrey), linear-gradient(45deg, lightgrey 25%, transparent 25%, transparent 75%, lightgrey 75%, lightgrey)',
                    'background-size': '60px 60px',
                    'background-position': '0 0, 30px 30px'
                });
            }
            parent.appendChild(canvas);
        }

        this.canvas = canvas;
        this.parent = parent;
    }
}