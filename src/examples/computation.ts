/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, Run, Progress } from 'mol-task'

async function test() {
    const t = Task.create('test', async () => 1);
    const r = await Run(t);
    console.log(r);
}

function delay(ms: number) {
    return new Promise(r => setTimeout(r, ms));
}

function messageTree(root: Progress.Node, prefix = ''): string {
    if (!root.children.length) return `${prefix}${root.progress.message}`;

    const newPrefix = prefix + '  |_ ';
    const subTree = root.children.map(c => messageTree(c, newPrefix));
    return `${prefix}${root.progress.message}\n${subTree.join('\n')}`;
}

function createTask<T>(delayMs: number, r: T): Task<T> {
    return Task.create('delayed', async ctx => {
        ctx.updateProgress('Processing delayed... ' + r);
        await delay(delayMs);
        if (ctx.needsYield) await ctx.yield({ message: 'hello from delayed... ' + r });
        return r;
    });
}

async function testObs() {
    const t = Task.create('test o', async ctx => {
        await delay(250);
        if (ctx.needsYield) await ctx.yield({ message: 'hi! 1' });
        await delay(125);
        if (ctx.needsYield) await ctx.yield({ message: 'hi! 2' });
        await delay(250);
        if (ctx.needsYield) await ctx.yield('hi! 3');

        ctx.updateProgress('Running children...');
        const c1 = ctx.runChild(createTask(250, 1));
        const c2 = ctx.runChild(createTask(500, 2));
        const c3 = ctx.runChild(createTask(750, 3));
        const r = await c1 + await c2 + await c3;
        if (ctx.needsYield) await ctx.yield({ message: 'Almost done...' });
        return r + 1;
    });
    const r = await Run(t, p => console.log(messageTree(p.root)), 250);
    console.log(r);
}

test();
testObs();