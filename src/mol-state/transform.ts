/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransformer } from './transformer';
import { UUID } from '../mol-util';
import { hashMurmur128o } from '../mol-data/util/hash-functions';
import { ParamDefinition as PD } from '../mol-util/param-definition';

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
            if (s !== void 0) {
                (a as any)[k] = s;
            } else {
                delete (a as any)[k];
            }
        }
        for (const k of Object.keys(a)) {
            const s = (b as any)[k], t = (a as any)[k];
            if (!!s === !!t) continue;
            changed = true;
            if (s !== void 0) {
                (a as any)[k] = s;
            } else {
                delete (a as any)[k];
            }
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

    export function withDependsOn(t: Transform, newDependsOn?: string | string[]): Transform {
        let dependsOn: string[] | undefined = void 0;
        if (newDependsOn) {
            dependsOn = typeof newDependsOn === 'string' ? [newDependsOn] : newDependsOn;
            if (dependsOn.length === 0) dependsOn = void 0;
            else dependsOn.sort();
        }
        return { ...t, dependsOn, version: UUID.create22() };
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

    /**
     * Compute the effective set of sibling-like dependencies for a transform.
     *
     * Combines (in order, de-duplicated):
     *   1. Explicit `t.dependsOn` (back-compat / non-param refs).
     *   2. Refs from `transformer.definition.getDependencies(params)` if defined.
     *   3. Refs collected from `PD.ValueRef` / `PD.DataRef` parameter values.
     *
     * Self-references and the root ref are filtered out. `globalCtx` is forwarded
     * to `params(undefined, globalCtx)` for schema acquisition. If the schema
     * can't be obtained (no `params` function or it throws), auto-derivation
     * falls back to a structural scan of parameter values for `{ ref, getValue }`
     * shaped objects.
     */
    export function getEffectiveDependsOn(t: Transform, globalCtx?: unknown): Ref[] {
        const out: Ref[] = [];
        const seen = new Set<string>();
        const add = (ref: string | undefined) => {
            if (!ref || ref === t.ref || ref === RootRef) return;
            if (seen.has(ref)) return;
            seen.add(ref);
            out.push(ref as Ref);
        };

        if (t.dependsOn) {
            for (const r of t.dependsOn) add(r);
        }

        const def = t.transformer.definition;
        const params = t.params as any;

        if (def.getDependencies && params) {
            try {
                const extra = def.getDependencies(params);
                if (extra) for (const r of extra) add(r);
            } catch {
                // Keep reconciliation robust if a user hook misbehaves.
            }
        }

        if (params) {
            let schema: PD.Params | undefined = void 0;
            if (def.params) {
                try {
                    schema = def.params(undefined as any, globalCtx) as PD.Params;
                } catch {
                    schema = void 0;
                }
            }
            if (schema) {
                const refs = PD.collectRefs(schema, params);
                refs.forEach(r => add(r));
            } else {
                collectStructuralRefs(params, add);
            }
        }

        return out;
    }

    function collectStructuralRefs(value: any, add: (ref: string) => void, depth = 0) {
        if (!value || typeof value !== 'object' || depth > 6) return;
        if (Array.isArray(value)) {
            for (const v of value) collectStructuralRefs(v, add, depth + 1);
            return;
        }
        if (typeof (value as any).ref === 'string' && typeof (value as any).getValue === 'function') {
            add((value as any).ref);
            return;
        }
        for (const k of Object.keys(value)) {
            collectStructuralRefs(value[k], add, depth + 1);
        }
    }

    const _emptyParams = {};
    /** Updates the version of the transform to be computed as hash of the parameters */
    export function setParamsHashVersion(t: Transform) {
        let version: string;
        try {
            version = hashMurmur128o(t.params ?? _emptyParams);
        } catch {
            const pToJson = t.transformer.definition.customSerialization
                ? t.transformer.definition.customSerialization.toJSON
                : _id;
            version = hashMurmur128o(pToJson(t.params ?? _emptyParams));
        }
        (t as { version: string }).version = version;
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
            ? transformer.definition.customSerialization.fromJSON
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