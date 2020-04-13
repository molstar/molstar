/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, Progress, Scheduler, MultistepTask, chunkedSubtask } from '../mol-task';
import { now } from '../mol-util/now';

export async function test1() {
    const t = Task.create('test', async () => 1);
    const r = await t.run();
    console.log(r);
}

function messageTree(root: Progress.Node, prefix = ''): string {
    const p = root.progress;
    if (!root.children.length) {
        if (p.isIndeterminate) return `${prefix}${p.taskName}: ${p.message}`;
        return `${prefix}${p.taskName}: [${p.current}/${p.max}] ${p.message}`;
    }

    const newPrefix = prefix + '  |_ ';
    const subTree = root.children.map(c => messageTree(c, newPrefix));
    if (p.isIndeterminate) return `${prefix}${p.taskName}: ${p.message}\n${subTree.join('\n')}`;
    return `${prefix}${p.taskName}: [${p.current}/${p.max}] ${p.message}\n${subTree.join('\n')}`;
}

function createTask<T>(delayMs: number, r: T): Task<T> {
    return Task.create('delayed value ' + r, async ctx => {
        ctx.update(`Processing delayed ${r} after ${delayMs}ms`, true);
        await Scheduler.delay(delayMs);
        if (ctx.shouldUpdate) await ctx.update({ message: `hello from delayed ${r} ${delayMs}` });
        return r;
    }, () => console.log('On abort called ' + r));
}

export function abortAfter(delay: number) {
    return Task.create('abort after ' + delay, async ctx => {
        await Scheduler.delay(delay);
        throw Task.Aborted('test');
        // if (ctx.shouldUpdate) await ctx.update({ message: 'hello from delayed... ' });
        // return r;
    });
}

export function testTree() {
    return Task.create('test o', async ctx => {
        await Scheduler.delay(250);
        if (ctx.shouldUpdate) await ctx.update({ message: 'hi! 1' });
        await Scheduler.delay(125);
        if (ctx.shouldUpdate) await ctx.update({ message: 'hi! 2' });
        await Scheduler.delay(250);
        if (ctx.shouldUpdate) await ctx.update('hi! 3');

        // ctx.update('Running children...', true);
        const c1 = createTask(250, 1).runAsChild(ctx);
        const c2 = createTask(500, 2).runAsChild(ctx);
        const c3 = createTask(750, 3).runAsChild(ctx);

        // await ctx.runChild(abortAfter(350));

        const r = await c1 + await c2 + await c3;
        if (ctx.shouldUpdate) await ctx.update({ message: 'Almost done...' });
        return r + 1;
    }, () => console.log('On abort O'));
}

export type ChunkedState = { i: number, current: number, total: number }

function processChunk(n: number, state: ChunkedState): number {
    const toProcess = Math.min(state.current + n, state.total);
    const start = state.current;
    for (let i = start; i < toProcess; i++) {
        for (let j = 0; j < 1000000; j++) {
            state.i += (i * j + 1 + state.i) % 1023;
            state.i = state.i % 1000;
        }
    }
    state.current = toProcess;
    return toProcess - start;
}

export const ms = MultistepTask('ms-task', ['step 1', 'step 2', 'step 3'], async (p: { i: number }, step, ctx) => {
    await step(0);

    const child = Task.create('chunked', async ctx => {
        const s = await chunkedSubtask(ctx, 25, { i: 0, current: 0, total: 125 }, processChunk, (ctx, s, p) => ctx.update('chunk test ' + p));
        return s.i;
    });

    await child.runAsChild(ctx);
    await Scheduler.delay(250);
    await step(1);
    await chunkedSubtask(ctx, 25, { i: 0, current: 0, total: 80 }, processChunk, (ctx, s, p) => ctx.update('chunk test ' + p));
    await Scheduler.delay(250);
    await step(2);
    await Scheduler.delay(250);
    return p.i + 3;
});


export function abortingObserver(p: Progress) {
    console.log(messageTree(p.root));
    if (now() - p.root.progress.startedTime > 1000) {
        p.requestAbort('test');
    }
}

export function logP(p: Progress) { console.log(messageTree(p.root)); }

async function test() {
    try {
        // const r = await Run(testTree(), p => console.log(messageTree(p.root)), 250);
        // const r = await Run(testTree(), abortingObserver, 250);
        // console.log(r);

        const m = await testTree().run(abortingObserver, 50);
        console.log(m);
    } catch (e) {
        console.error(e);
    }
}

test();
// testObs();