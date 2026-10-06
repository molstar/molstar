/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { isDebugMode } from '../../mol-util/debug';
import { GraphicsRenderable } from '../../mol-gl/renderable';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { WebGPUContext } from '../../mol-gl/webgpu/context';
import { WebGPUHelperScene } from './webgpu-scene';

export type DebugHelperScene = Pick<Scene, 'add' | 'remove' | 'clear' | 'update' | 'commit'>;
export interface DebugHelperParent extends Pick<Scene, 'boundingSphere' | 'boundingSphereVisible' | 'has'> {
    forEach(callback: (value: Pick<GraphicsRenderable, 'values'>, object: GraphicsRenderObject) => void): void;
}

export interface DebugHelper<T extends {} = {}, S extends DebugHelperScene = Scene> {
    readonly scene: S;
    update(): void;
    syncVisibility(): void;
    clear(): void;
    readonly isEnabled: boolean;
    readonly props: T;
    setProps(props: Partial<T>): void;
}

export class DebugRegistry<S extends DebugHelperScene = Scene, C = WebGLContext, P extends DebugHelperParent = Scene> {
    readonly ctx: C;
    readonly parent: P;

    private readonly entries = new Map<string, DebugHelper<{}, S>>();

    constructor(ctx: C, parent: P) {
        this.ctx = ctx;
        this.parent = parent;
    }

    register<T extends {}>(name: string, entry: DebugHelper<T, S>) {
        if (this.entries.has(name)) {
            if (isDebugMode) {
                console.warn(`Debug helper with name '${name}' already exists, replacing.`);
            }
            this.entries.get(name)!.clear();
        }
        this.entries.set(name, entry);
    }

    unregister(name: string) {
        const entry = this.entries.get(name);
        if (entry) {
            entry.clear();
            this.entries.delete(name);
        }
    }

    get scenes(): S[] {
        return Array.from(this.entries.values()).map(e => e.scene);
    }

    update() {
        this.entries.forEach(entry => {
            if (entry.isEnabled) entry.update();
        });
    }

    syncVisibility() {
        this.entries.forEach(entry => {
            entry.syncVisibility();
        });
    }

    clear() {
        this.entries.forEach(entry => {
            entry.clear();
        });
    }

    get isEnabled() {
        let enabled = false;
        this.entries.forEach(entry => {
            if (entry.isEnabled) enabled = true;
        });
        return enabled;
    }

    setProps<T extends {}>(props: Partial<T>) {
        this.entries.forEach(entry => {
            entry.setProps(props);
        });
    }
}

/** Debug helpers use CPU scene ownership with the same geometry builders as WebGL. */
export class WebGPUDebugRegistry extends DebugRegistry<WebGPUHelperScene, WebGPUContext, DebugHelperParent> {
    getRenderObjects() {
        this.update();
        this.syncVisibility();
        return this.scenes.flatMap(scene => scene.renderObjects).filter(object => object.state.visible);
    }
}
