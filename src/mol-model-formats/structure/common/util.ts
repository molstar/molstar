/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { memoize1 } from '../../../mol-util/memoize';

const ChainIdAlphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

function _getChainId (index: number) {
    const n = ChainIdAlphabet.length;
    let j = index;
    let k = 0;
    let chainId = ChainIdAlphabet[j % n];
    while (j >= n) {
        j = Math.floor(j / n);
        chainId += ChainIdAlphabet[j % n];
        k += 1;
    }
    if (k >= 5) {
        console.warn('getChainId overflow');
    }
    return chainId;
}

export const getChainId = memoize1(_getChainId);