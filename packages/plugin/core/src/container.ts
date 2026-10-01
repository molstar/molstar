/**
 * Copyright (c) 2024-25 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

    /**
     * options.checkeredCanvasBackground has no effect. Use canvas3d.checkeredTransparentBackground instead.
     * TODO: remove in v6
     */
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
            parent.appendChild(canvas);
        }

        this.canvas = canvas;
        this.parent = parent;
    }
}