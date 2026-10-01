/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Scene } from '@molstar/graphics/gl/scene';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import { isDebugMode } from '@molstar/core/util/debug';

export interface DebugHelper<T extends {} = {}> {
    readonly scene: Scene;
    update(): void;
    syncVisibility(): void;
    clear(): void;
    readonly isEnabled: boolean;
    readonly props: T;
    setProps(props: Partial<T>): void;
}

export class DebugRegistry {
    readonly ctx: WebGLContext;
    readonly parent: Scene;

    private readonly entries = new Map<string, DebugHelper>();

    constructor(ctx: WebGLContext, parent: Scene) {
        this.ctx = ctx;
        this.parent = parent;
    }

    register<T extends {}>(name: string, entry: DebugHelper<T>) {
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

    get scenes(): Scene[] {
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
