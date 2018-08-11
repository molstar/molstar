/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export async function retryIf<T>(promiseProvider: () => Promise<T>, params: {
    retryThenIf?: (result: T) => boolean,
    retryCatchIf?: (error: any) => boolean,
    retryCount: number
}) {
    let count = 0;
    while (count <= params.retryCount) {
        try {
            const result = await promiseProvider();
            if (params.retryThenIf && params.retryThenIf(result)) {
                count++;
                continue;
            }
            return result;
        } catch (e) {
            if (!params.retryCatchIf || params.retryCatchIf(e)) {
                count++;
                continue;
            }
            throw e;
        }
    }

    throw new Error('Maximum retry count exceeded.');
}