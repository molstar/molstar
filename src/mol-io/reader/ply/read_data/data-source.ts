/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from LiteMol
 */

import { Task, RuntimeContext } from 'mol-task';

export function readFromFile(file: File) {
    return <Task<number | string>>readFromFileInternal(file);
}


async function processFile(ctx: RuntimeContext, e: any) {
    const data = (e.target as FileReader).result;
    return  data as string;
}

function readData(ctx: RuntimeContext, action: string, data: XMLHttpRequest | FileReader): Promise<any> {
    return new Promise<any>((resolve, reject) => {
        data.onerror = (e: any) => {
            const error = (<FileReader>e.target).error;
            reject(error ? error : 'Failed.');
        };

        data.onabort = () => reject(Task.Aborted(''));

        data.onprogress = (e: ProgressEvent) => {
            if (e.lengthComputable) {
                ctx.update({ message: action, isIndeterminate: false, current: e.loaded, max: e.total });
            } else {
                ctx.update({ message: `${action} ${(e.loaded / 1024 / 1024).toFixed(2)} MB`, isIndeterminate: true });
            }
        }
        data.onload = (e: any) => resolve(e);
    });
}

function readFromFileInternal(file: File): Task<string | number> {
    let reader: FileReader | undefined = void 0;
    return Task.create('Read File', async ctx => {
        try {
            reader = new FileReader();
            reader.readAsBinaryString(file);

            ctx.update({ message: 'Opening file...', canAbort: true });
            const e = await readData(ctx, 'Reading...', reader);
            const result = processFile(ctx, e);
            return result;
        } finally {
            reader = void 0;
        }
    }, () => {
        if (reader) reader.abort();
    });
}
