/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, Run, Progress, Scheduler } from 'mol-task'

async function test() {
    const t = Task.create('test', async () => 1);
    const r = await Run(t);
    console.log(r);
}

function messageTree(root: Progress.Node, prefix = ''): string {
    if (!root.children.length) return `${prefix}${root.progress.taskName}: ${root.progress.message}`;

    const newPrefix = prefix + '  |_ ';
    const subTree = root.children.map(c => messageTree(c, newPrefix));
    return `${prefix}${root.progress.taskName}: ${root.progress.message}\n${subTree.join('\n')}`;
}

function createTask<T>(delayMs: number, r: T): Task<T> {
    return Task.create('delayed value ' + r, async ctx => {
        ctx.update('Processing delayed... ' + r, true);
        await Scheduler.delay(delayMs);
        if (ctx.shouldUpdate) await ctx.update({ message: 'hello from delayed... ' });
        return r;
    });
}

async function testObs() {
    const t = Task.create('test o', async ctx => {
        await Scheduler.delay(250);
        if (ctx.shouldUpdate) await ctx.update({ message: 'hi! 1' });
        await Scheduler.delay(125);
        if (ctx.shouldUpdate) await ctx.update({ message: 'hi! 2' });
        await Scheduler.delay(250);
        if (ctx.shouldUpdate) await ctx.update('hi! 3');

        ctx.update('Running children...', true);
        const c1 = ctx.runChild(createTask(250, 1));
        const c2 = ctx.runChild(createTask(500, 2));
        const c3 = ctx.runChild(createTask(750, 3));
        const r = await c1 + await c2 + await c3;
        if (ctx.shouldUpdate) await ctx.update({ message: 'Almost done...' });
        return r + 1;
    });
    const r = await Run(t, p => console.log(messageTree(p.root)), 250);
    console.log(r);
}

test();
testObs();