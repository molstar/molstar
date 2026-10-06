/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Scene } from '../../mol-gl/scene';
import { Mat4 } from '../../mol-math/linear-algebra';

/** CPU helper geometry ownership; native rendering uploads these objects directly. */
export class WebGPUHelperScene {
    readonly view = Mat4.identity();
    readonly renderObjects: GraphicsRenderObject[] = [];
    add(object: GraphicsRenderObject) { if (!this.renderObjects.includes(object)) this.renderObjects.push(object); }
    remove(object: GraphicsRenderObject) { const index = this.renderObjects.indexOf(object); if (index >= 0) this.renderObjects.splice(index, 1); }
    clear() { this.renderObjects.length = 0; }
    commit() { return true; }
    update: Scene['update'] = () => {};
}
