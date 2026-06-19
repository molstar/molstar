/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { State, StateObject, StateObjectCell, StateTransform, StateTransformer, StateTreeCycleError } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

interface TypeInfo { name: string; typeClass: 'Root' | 'Data' }
const Create = StateObject.factory<TypeInfo>();

class Root extends Create({ name: 'Root', typeClass: 'Root' }) { }
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) { }

const NS = 'state-deps-spec';
let counter = 0;
const uniq = (s: string) => `${s}-${counter++}`;

function newState() {
    return State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
}

/** Plain leaf created from Root with a number param. */
function constLeaf() {
    return StateTransformer.create<Root, Leaf, { value: number }>(NS, {
        name: uniq('const-leaf'),
        from: [Root],
        to: [Leaf],
        display: { name: 'Const Leaf' },
        params: () => ({ value: PD.Numeric(0) }) as any,
        apply({ params }) { return new Leaf({ value: params.value }); },
        update({ oldParams, newParams }) {
            return oldParams.value === newParams.value
                ? StateTransformer.UpdateResult.Unchanged
                : StateTransformer.UpdateResult.Recreate;
        }
    });
}

/** Leaf whose value is read from a single explicit dependsOn ref. */
function deriveFromDep(depRef: string) {
    return StateTransformer.create<Root, Leaf, {}>(NS, {
        name: uniq('derive-from-dep'),
        from: [Root],
        to: [Leaf],
        display: { name: 'Derive From Dep' },
        params: () => ({}) as any,
        apply({ dependencies }) {
            const dep = dependencies?.[depRef] as Leaf;
            if (!dep) throw new Error('missing dep');
            return new Leaf({ value: dep.data.value + 100 });
        },
        update({ b, dependencies }) {
            const dep = dependencies?.[depRef] as Leaf;
            if (!dep) throw new Error('missing dep');
            (b.data as { value: number }).value = dep.data.value + 100;
            return StateTransformer.UpdateResult.Updated;
        }
    });
}

describe('State dependencies - linking', () => {
    it('explicit dependsOn establishes an edge and passes the dep object to apply', async () => {
        const state = newState();
        const A = constLeaf();
        const B = deriveFromDep('leaf-a');

        const builder = state.build();
        builder.toRoot<Root>().apply(A as any, { value: 7 }, { ref: 'leaf-a' });
        builder.toRoot<Root>().apply(B as any, {}, { ref: 'leaf-b', dependsOn: ['leaf-a'] });
        await state.runTask(state.updateTree(builder));

        const b = state.cells.get('leaf-b')!;
        expect(b.dependencies.dependsOn.map(c => c.transform.ref)).toEqual(['leaf-a']);
        expect((b.obj as Leaf).data.value).toBe(107);

        const a = state.cells.get('leaf-a')!;
        expect(a.dependencies.dependentBy.map(c => c.transform.ref)).toEqual(['leaf-b']);
    });

    it('re-evaluates dependents when the source updates', async () => {
        const state = newState();
        const A = constLeaf();
        const B = deriveFromDep('leaf-a');

        const builder1 = state.build();
        builder1.toRoot<Root>().apply(A as any, { value: 1 }, { ref: 'leaf-a' });
        builder1.toRoot<Root>().apply(B as any, {}, { ref: 'leaf-b', dependsOn: ['leaf-a'] });
        await state.runTask(state.updateTree(builder1));
        expect((state.cells.get('leaf-b')!.obj as Leaf).data.value).toBe(101);

        const builder2 = state.build();
        builder2.to('leaf-a').update({ value: 5 });
        await state.runTask(state.updateTree(builder2));
        expect((state.cells.get('leaf-b')!.obj as Leaf).data.value).toBe(105);
    });

    it('throws when an explicit dependsOn references a non-existent transform', async () => {
        const state = newState();
        const B = deriveFromDep('missing-ref');

        const builder = state.build();
        builder.toRoot<Root>().apply(B as any, {}, { ref: 'leaf-b', dependsOn: ['missing-ref'] });
        await expect(state.runTask(state.updateTree(builder))).rejects.toThrow(/non-existent transform/);
    });

    it('honors getDependencies(params) and relinks when params change', async () => {
        const state = newState();
        const A = constLeaf();
        const A2 = constLeaf();

        const PickViaParams = StateTransformer.create<Root, Leaf, { which: string }>(NS, {
            name: uniq('pick-via-params'),
            from: [Root],
            to: [Leaf],
            display: { name: 'Pick' },
            params: () => ({ which: PD.Text('leaf-a') }) as any,
            getDependencies(params) { return params.which ? [params.which as StateTransform.Ref] : []; },
            apply({ params, dependencies }) {
                const dep = dependencies?.[params.which] as Leaf;
                return new Leaf({ value: dep ? dep.data.value : -1 });
            },
            update({ b, newParams, dependencies }) {
                const dep = dependencies?.[newParams.which] as Leaf;
                (b.data as { value: number }).value = dep ? dep.data.value : -1;
                return StateTransformer.UpdateResult.Updated;
            }
        });

        const builder = state.build();
        builder.toRoot<Root>().apply(A as any, { value: 11 }, { ref: 'leaf-a' });
        builder.toRoot<Root>().apply(A2 as any, { value: 22 }, { ref: 'leaf-a2' });
        builder.toRoot<Root>().apply(PickViaParams as any, { which: 'leaf-a' }, { ref: 'pick' });
        await state.runTask(state.updateTree(builder));

        const pick = state.cells.get('pick')!;
        expect(pick.dependencies.dependsOn.map(c => c.transform.ref)).toEqual(['leaf-a']);
        expect((pick.obj as Leaf).data.value).toBe(11);

        const update = state.build();
        update.to('pick').update({ which: 'leaf-a2' });
        await state.runTask(state.updateTree(update));

        const pick2 = state.cells.get('pick')!;
        expect(pick2.dependencies.dependsOn.map(c => c.transform.ref)).toEqual(['leaf-a2']);
        expect((pick2.obj as Leaf).data.value).toBe(22);
        // Old source no longer reverse-linked.
        expect(state.cells.get('leaf-a')!.dependencies.dependentBy.length).toBe(0);
        expect(state.cells.get('leaf-a2')!.dependencies.dependentBy.map(c => c.transform.ref)).toEqual(['pick']);
    });

    it('auto-collects refs from PD.ValueRef parameter values', async () => {
        const state = newState();
        const A = constLeaf();

        const ViaValueRef = StateTransformer.create<Root, Leaf, { target: { ref: string, getValue: () => Leaf } }>(NS, {
            name: uniq('via-value-ref'),
            from: [Root],
            to: [Leaf],
            display: { name: 'Via ValueRef' },
            params: () => ({
                target: PD.ValueRef<Leaf>(() => [], (ref, getData) => getData(ref))
            }) as any,
            apply({ params, dependencies }) {
                const dep = dependencies?.[params.target.ref] as Leaf;
                return new Leaf({ value: dep ? dep.data.value * 2 : -1 });
            }
        });

        const builder = state.build();
        builder.toRoot<Root>().apply(A as any, { value: 9 }, { ref: 'leaf-a' });
        builder.toRoot<Root>().apply(ViaValueRef as any, {
            target: { ref: 'leaf-a', getValue: () => null as any }
        }, { ref: 'vr' });
        await state.runTask(state.updateTree(builder));

        const vr = state.cells.get('vr')!;
        expect(vr.dependencies.dependsOn.map(c => c.transform.ref)).toEqual(['leaf-a']);
        expect((vr.obj as Leaf).data.value).toBe(18);
    });

    it('falls back to a structural scan when the schema is unavailable', async () => {
        const state = newState();
        const A = constLeaf();

        // No `def.params` - params normalization will drop unknown fields at
        // evaluation time, but link-time collection (via the structural
        // fallback) still happens against the original transform.params.
        const Structural = StateTransformer.create<Root, Leaf, any>(NS, {
            name: uniq('structural'),
            from: [Root],
            to: [Leaf],
            display: { name: 'Structural' },
            apply({ dependencies }) {
                const ref = dependencies ? Object.keys(dependencies)[0] : undefined;
                const dep = ref ? dependencies![ref] as Leaf : undefined;
                return new Leaf({ value: dep ? dep.data.value + 1000 : -1 });
            }
        });

        const builder = state.build();
        builder.toRoot<Root>().apply(A as any, { value: 3 }, { ref: 'leaf-a' });
        builder.toRoot<Root>().apply(Structural as any, {
            link: { ref: 'leaf-a', getValue: () => null }
        }, { ref: 'struct' });
        await state.runTask(state.updateTree(builder));

        const s = state.cells.get('struct')!;
        expect(s.dependencies.dependsOn.map(c => c.transform.ref)).toEqual(['leaf-a']);
        expect((s.obj as Leaf).data.value).toBe(1003);
    });

    it('filters out self and root refs from getDependencies', async () => {
        const state = newState();

        const SelfRef = StateTransformer.create<Root, Leaf, {}>(NS, {
            name: uniq('self-ref'),
            from: [Root],
            to: [Leaf],
            display: { name: 'Self Ref' },
            params: () => ({}) as any,
            getDependencies() { return ['self', StateTransform.RootRef as any]; },
            apply() { return new Leaf({ value: 42 }); }
        });

        const builder = state.build();
        builder.toRoot<Root>().apply(SelfRef as any, {}, { ref: 'self' });
        await state.runTask(state.updateTree(builder));

        const cell = state.cells.get('self')!;
        expect(cell.dependencies.dependsOn.length).toBe(0);
        expect((cell.obj as Leaf).data.value).toBe(42);
    });
});

describe('State dependencies - cycle detection', () => {
    it('throws StateTreeCycleError for a direct A → B → A cycle', async () => {
        const state = newState();

        // Two transformers, each declaring a getDependencies pointing at the other.
        const A = StateTransformer.create<Root, Leaf, {}>(NS, {
            name: uniq('cycle-a'),
            from: [Root],
            to: [Leaf],
            display: { name: 'Cycle A' },
            params: () => ({}) as any,
            getDependencies() { return ['cyc-b' as any]; },
            apply() { return new Leaf({ value: 0 }); }
        });
        const B = StateTransformer.create<Root, Leaf, {}>(NS, {
            name: uniq('cycle-b'),
            from: [Root],
            to: [Leaf],
            display: { name: 'Cycle B' },
            params: () => ({}) as any,
            getDependencies() { return ['cyc-a' as any]; },
            apply() { return new Leaf({ value: 0 }); }
        });

        const builder = state.build();
        builder.toRoot<Root>().apply(A as any, {}, { ref: 'cyc-a' });
        builder.toRoot<Root>().apply(B as any, {}, { ref: 'cyc-b' });

        let caught: unknown;
        try {
            await state.runTask(state.updateTree(builder));
        } catch (e) { caught = e; }
        expect(caught).toBeInstanceOf(StateTreeCycleError);
        const cycle = (caught as StateTreeCycleError).cycle;
        expect(cycle[0]).toBe(cycle[cycle.length - 1]);
        expect(cycle).toEqual(expect.arrayContaining(['cyc-a', 'cyc-b']));
    });
});

describe('State dependencies - deferred resolution', () => {
    /** Force evaluation order: place dependent subtree first under root so
     *  tree pre-order visits it before its dependency. */
    it('resolves cross-subtree deps even when the dependent is scheduled first', async () => {
        const state = newState();
        const A = constLeaf();
        const B = deriveFromDep('leaf-a');

        const builder = state.build();
        // B added FIRST - so its subtree comes before A's in tree pre-order.
        builder.toRoot<Root>().apply(B as any, {}, { ref: 'leaf-b', dependsOn: ['leaf-a'] });
        builder.toRoot<Root>().apply(A as any, { value: 4 }, { ref: 'leaf-a' });

        await state.runTask(state.updateTree(builder));

        expect((state.cells.get('leaf-a')!.obj as Leaf).data.value).toBe(4);
        expect((state.cells.get('leaf-b')!.obj as Leaf).data.value).toBe(104);
    });

    it('propagates a clear error when a dep has errored and cannot resolve', async () => {
        const state = newState();

        const Boom = StateTransformer.create<Root, Leaf, {}>(NS, {
            name: uniq('boom'),
            from: [Root],
            to: [Leaf],
            display: { name: 'Boom' },
            params: () => ({}) as any,
            apply() { throw new Error('intentional'); }
        });
        const B = deriveFromDep('boom');

        const builder = state.build();
        builder.toRoot<Root>().apply(Boom as any, {}, { ref: 'boom' });
        builder.toRoot<Root>().apply(B as any, {}, { ref: 'leaf-b', dependsOn: ['boom'] });
        // The state surfaces transform errors via console.error; suppress the noise.
        const err = jest.spyOn(console, 'error').mockImplementation(() => {});
        try {
            await state.runTask(state.updateTree(builder));
        } finally {
            err.mockRestore();
        }

        const b: StateObjectCell = state.cells.get('leaf-b')!;
        expect(b.status).toBe('error');
        expect(b.errorText).toMatch(/Unresolved dependency|missing dep|intentional/);
    });
});
