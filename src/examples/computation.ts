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

function messageTree(root: Progress.Node, prefix = ''): string {
    if (!root.children.length) return `${prefix}${root.progress.message}`;

    const newPrefix = prefix + '  |_ ';
    const subTree = root.children.map(c => messageTree(c, newPrefix));
    return `${prefix}${root.progress.message}\n${subTree.join('\n')}`;
}

function createTask<T>(delay: number, r: T): Task<T> {
    return Task.create('delayed', async ctx => {
        await new Promise(r => setTimeout(r, delay));
        if (ctx.requiresUpdate) await ctx.update({ message: 'hello from delayed...' });
        return r;
    });
}

async function testObs() {
    const t = Task.create('test o', async ctx => {
        await new Promise(r => setTimeout(r, 250));
        if (ctx.requiresUpdate) await ctx.update({ message: 'hi! 1' });
        await new Promise(r => setTimeout(r, 125));
        if (ctx.requiresUpdate) await ctx.update({ message: 'hi! 2' });
        await new Promise(r => setTimeout(r, 250));
        if (ctx.requiresUpdate) await ctx.update({ message: 'hi! 3' });

        const r = await ctx.runChild({ message: 'Running child!' }, createTask(250, 100));
        if (ctx.requiresUpdate) await ctx.update({ message: 'Almost done...' });
        return r + 1;
    });
    const r = await Run(t, p => console.log(messageTree(p.root)), 250);
    console.log(r);
}

test();
testObs();