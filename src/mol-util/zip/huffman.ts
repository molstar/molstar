/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 */

import { NumberArray } from '../type-helpers';

export type HufTree = {
    lit: number,
    f: number,
    l?: HufTree,
    r?: HufTree,
    d: number
}
export function _hufTree(hst: NumberArray, tree: number[], MAXL: number) {
    const list: HufTree[] = [];
    const hl = hst.length, tl = tree.length;
    for(let i = 0; i < tl; i += 2) {
        tree[i] = 0;
        tree[i + 1] = 0;
    }
    for(let i = 0; i < hl; i++) if(hst[i] !== 0) list.push({ lit: i, f: hst[i], d: undefined as any });
    const end = list.length, l2 = list.slice(0);
    if(end === 0) return 0;  // empty histogram (usually for dist)
    if(end === 1) {
        const lit = list[0].lit, l2 = lit === 0 ? 1 : 0;
        tree[(lit << 1) + 1] = 1;
        tree[(l2 << 1) + 1] = 1;
        return 1;
    }
    list.sort(function(a, b){return a.f - b.f;});
    let a = list[0], b = list[1], i0 = 0, i1 = 1, i2 = 2;
    list[0] = {
        lit: -1,
        f: a.f + b.f,
        l: a,
        r: b,
        d: 0
    };
    while(i1 !== end - 1) {
        if(i0 !== i1 && (i2 === end || list[i0].f < list[i2].f)) {
            a = list[i0++];
        } else {
            a = list[i2++];
        }
        if(i0 !== i1 && (i2 === end || list[i0].f < list[i2].f)) {
            b = list[i0++];
        } else {
            b = list[i2++];
        }
        list[i1++] = {
            lit: -1,
            f: a.f + b.f,
            l: a,
            r: b,
            d: undefined as any
        };
    }
    let maxl = setDepth(list[i1 - 1], 0);
    if(maxl > MAXL) {
        restrictDepth(l2, MAXL, maxl);
        maxl = MAXL;
    }
    for(let i = 0; i < end; i++) tree[(l2[i].lit << 1) + 1] = l2[i].d;
    return maxl;
}

function setDepth(t: HufTree, d: number): number {
    if(t.lit !== -1) {
        t.d = d;
        return d;
    }
    return Math.max(setDepth(t.l!, d + 1), setDepth(t.r!, d + 1));
}

function restrictDepth(dps: HufTree[], MD: number, maxl: number) {
    let i = 0, bCost = 1 << (maxl - MD), dbt = 0;
    dps.sort(function(a: HufTree, b: HufTree){return b.d === a.d ? a.f - b.f : b.d - a.d;});

    for(i = 0; i < dps.length; i++) {
        if(dps[i].d > MD) {
            const od = dps[i].d;
            dps[i].d = MD;
            dbt += bCost - (1 << (maxl - od));
        } else {
            break;
        }
    }
    dbt = dbt >>> (maxl - MD);
    while(dbt > 0) {
        const od = dps[i].d;
        if(od < MD) {
            dps[i].d++;
            dbt -= (1 << (MD - od - 1));
        } else {
            i++;
        }
    }
    for(; i >= 0; i--) {
        if(dps[i].d === MD && dbt < 0) {
            dps[i].d--;
            dbt++;
        }
    }
    if(dbt !== 0) console.log('debt left');
}