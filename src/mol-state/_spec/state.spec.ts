/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { State, StateObject, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';

interface TypeInfo { name: string; typeClass: 'Root' | 'Data' }
const Create = StateObject.factory<TypeInfo>();

class Root extends Create({ name: 'Root', typeClass: 'Root' }) { }
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) { }

const NS = 'state-dispose-spec';
let counter = 0;

function leafTransformer(spy: () => void) {
    return StateTransformer.create<Root, Leaf, { value: number }>(NS, {
        name: `create-leaf-${counter++}`,
        from: [Root],
        to: [Leaf],
        display: { name: 'Create Leaf' },
        params: () => ({} as any),
        apply({ params }) { return new Leaf({ value: params.value }); },
        dispose() { spy(); }
    });
}

function chainedTransformer(spy: () => void) {
    return StateTransformer.create<Leaf, Leaf, {}>(NS, {
        name: `chained-leaf-${counter++}`,
        from: [Leaf],
        to: [Leaf],
        display: { name: 'Chained Leaf' },
        apply({ a }) { return new Leaf({ value: a.data.value + 1 }); },
        dispose() { spy(); }
    });
}

function newState() {
    return State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
}

describe('State.dispose', () => {
    it('calls transformer.dispose for every live cell', async () => {
        const leafSpy = jest.fn();
        const chainSpy = jest.fn();
        const A = leafTransformer(leafSpy);
        const B = chainedTransformer(chainSpy);

        const state = newState();
        const builder = state.build();
        builder.toRoot<Root>().apply(A as any, { value: 1 }).apply(B as any, {});
        await state.runTask(state.updateTree(builder));

        // root + 2 transformer outputs.
        expect(state.cells.size).toBe(3);

        state.dispose();

        expect(leafSpy).toHaveBeenCalledTimes(1);
        expect(chainSpy).toHaveBeenCalledTimes(1);
    });

    it('disposes all sibling subtrees', async () => {
        const spyA = jest.fn();
        const spyB = jest.fn();
        const A = leafTransformer(spyA);
        const B = leafTransformer(spyB);

        const state = newState();
        const builder = state.build();
        builder.toRoot<Root>().apply(A as any, { value: 1 });
        builder.toRoot<Root>().apply(B as any, { value: 2 });
        await state.runTask(state.updateTree(builder));

        state.dispose();

        expect(spyA).toHaveBeenCalledTimes(1);
        expect(spyB).toHaveBeenCalledTimes(1);
    });

    it('does not throw when a transformer dispose throws', async () => {
        const goodSpy = jest.fn();
        const Throwing = StateTransformer.create<Root, Leaf, { value: number }>(NS, {
            name: `throwing-leaf-${counter++}`,
            from: [Root],
            to: [Leaf],
            display: { name: 'Throwing Leaf' },
            apply({ params }) { return new Leaf({ value: params.value }); },
            dispose() { throw new Error('boom'); }
        });
        const Good = leafTransformer(goodSpy);

        const state = newState();
        const builder = state.build();
        builder.toRoot<Root>().apply(Throwing as any, { value: 1 });
        builder.toRoot<Root>().apply(Good as any, { value: 2 });
        await state.runTask(state.updateTree(builder));

        const warn = jest.spyOn(console, 'warn').mockImplementation(() => {});
        try {
            expect(() => state.dispose()).not.toThrow();
        } finally {
            warn.mockRestore();
        }
        expect(goodSpy).toHaveBeenCalledTimes(1);
    });

    it('is a no-op for transformers without a dispose definition', async () => {
        const NoDispose = StateTransformer.create<Root, Leaf, { value: number }>(NS, {
            name: `no-dispose-${counter++}`,
            from: [Root],
            to: [Leaf],
            display: { name: 'No-dispose Leaf' },
            apply({ params }) { return new Leaf({ value: params.value }); }
        });

        const state = newState();
        const builder = state.build();
        builder.toRoot<Root>().apply(NoDispose as any, { value: 1 });
        await state.runTask(state.updateTree(builder));

        expect(() => state.dispose()).not.toThrow();
    });
});
