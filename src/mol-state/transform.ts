/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer } from './transformer';
import { UUID } from '../mol-util';

export { Transform as StateTransform };

interface Transform<T extends StateTransformer = StateTransformer> {
    readonly parent: Transform.Ref,
    readonly transformer: T,
    readonly state: Transform.State,
    readonly tags?: string[],
    readonly ref: Transform.Ref,
    /**
     * Sibling-like dependency
     * Do NOT make a cell dependent on its ancestor.
     */
    readonly dependsOn?: Transform.Ref[],
    readonly params?: StateTransformer.Params<T>,
    readonly version: string
}

namespace Transform {
    export type Ref = string
    export type Transformer<T extends Transform> = T extends Transform<infer S> ? S : never

    export const RootRef = '-=root=-' as Ref;

    export interface State {
        // is the cell shown in the UI
        isGhost?: boolean,
        // can the corresponding be deleted by the user.
        isLocked?: boolean,
        // is the representation associated with the cell hidden
        isHidden?: boolean,
        // is the tree node collapsed?
        isCollapsed?: boolean
    }

    export function areStatesEqual(a: State, b: State) {
        return !!a.isHidden !== !!b.isHidden || !!a.isCollapsed !== !!b.isCollapsed
            || !!a.isGhost !== !!b.isGhost || !!a.isLocked !== !!b.isLocked;
    }

    export function isStateChange(a: State, b?: Partial<State>) {
        if (!b) return false;
        if (typeof b.isCollapsed !== 'undefined' && a.isCollapsed !== b.isCollapsed) return true;
        if (typeof b.isHidden !== 'undefined' && a.isHidden !== b.isHidden) return true;
        if (typeof b.isGhost !== 'undefined' && a.isGhost !== b.isGhost) return true;
        if (typeof b.isLocked !== 'undefined' && a.isLocked !== b.isLocked) return true;
        return false;
    }

    export function assignState(a: State, b?: Partial<State>): boolean {
        if (!b) return false;

        let changed = false;
        for (const k of Object.keys(b)) {
            const s = (b as any)[k], t = (a as any)[k];
            if (!!s === !!t) continue;
            changed = true;
            (a as any)[k] = s;
        }
        return changed;
    }

    export function syncState(a: State, b?: Partial<State>): boolean {
        if (!b) return false;

        let changed = false;
        for (const k of Object.keys(b)) {
            const s = (b as any)[k], t = (a as any)[k];
            if (!!s === !!t) continue;
            changed = true;
            (a as any)[k] = s;
        }
        for (const k of Object.keys(a)) {
            const s = (b as any)[k], t = (a as any)[k];
            if (!!s === !!t) continue;
            changed = true;
            (a as any)[k] = s;
        }
        return changed;
    }

    export interface Options {
        ref?: string,
        tags?: string | string[],
        state?: State,
        dependsOn?: Ref[]
    }

    export function create<T extends StateTransformer>(parent: Ref, transformer: T, params?: StateTransformer.Params<T>, options?: Options): Transform<T> {
        const ref = options && options.ref ? options.ref : UUID.create22() as string as Ref;
        let tags: string[] | undefined = void 0;
        if (options && options.tags) {
            tags = typeof options.tags === 'string' ? [options.tags] : options.tags;
            if (tags.length === 0) tags = void 0;
            else tags.sort();
        }
        return {
            parent,
            transformer,
            state: options?.state || { },
            tags,
            ref,
            dependsOn: options && options.dependsOn,
            params,
            version: UUID.create22()
        };
    }

    export function withParams(t: Transform, params: any): Transform {
        return { ...t, params, version: UUID.create22() };
    }

    export function withState(t: Transform, state?: Partial<State>): Transform {
        if (!state) return t;
        return { ...t, state: { ...t.state, ...state } };
    }

    export function withTags(t: Transform, newTags?: string | string[]): Transform {
        let tags: string[] | undefined = void 0;
        if (newTags) {
            tags = typeof newTags === 'string' ? [newTags] : newTags;
            if (tags.length === 0) tags = void 0;
            else tags.sort();
        }
        return { ...t, tags, version: UUID.create22() };
    }

    export function withParent(t: Transform, parent: Ref): Transform {
        return { ...t, parent, version: UUID.create22() };
    }

    export function createRoot(state?: State): Transform {
        return create(RootRef, StateTransformer.ROOT, {}, { ref: RootRef, state });
    }

    export function hasTag(t: Transform, tag: string) {
        if (!t.tags) return false;
        return t.tags.indexOf(tag) >= 0;
    }

    export function hasTags(t: Transform, tags: string | string[]) {
        if (!t.tags) return typeof tags !== 'string' && tags.length === 0;
        if (typeof tags === 'string') return hasTag(t, tags);
        for (const tag of tags) {
            if (t.tags.indexOf(tag) < 0) return false;
        }
        return true;
    }

    export interface Serialized {
        parent: string,
        transformer: string,
        params: any,
        state?: State,
        tags?: string[],
        isDecorator?: boolean,
        ref: string,
        dependsOn?: string[]
        version: string
    }

    function _id(x: any) { return x; }
    export function toJSON(t: Transform): Serialized {
        const pToJson = t.transformer.definition.customSerialization
            ? t.transformer.definition.customSerialization.toJSON
            : _id;
        let state: any = void 0;
        for (const k of Object.keys(t.state)) {
            const s = (t.state as any)[k];
            if (!s) continue;
            if (!state) state = { };
            state[k] = true;
        }
        return {
            parent: t.parent,
            transformer: t.transformer.id,
            params: t.params ? pToJson(t.params) : void 0,
            state,
            tags: t.tags,
            ref: t.ref,
            dependsOn: t.dependsOn,
            version: t.version
        };
    }

    export function fromJSON(t: Serialized): Transform {
        const transformer = StateTransformer.get(t.transformer);
        const pFromJson = transformer.definition.customSerialization
            ? transformer.definition.customSerialization.toJSON
            : _id;
        return {
            parent: t.parent as Ref,
            transformer,
            params: t.params ? pFromJson(t.params) : void 0,
            state: t.state || { },
            tags: t.tags,
            ref: t.ref as Ref,
            dependsOn: t.dependsOn,
            version: t.version
        };
    }
}