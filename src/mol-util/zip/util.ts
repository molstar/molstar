/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 */

export const U = (function(){
    const u16 = Uint16Array, u32 = Uint32Array;
    return {
        next_code : new u16(16),
        bl_count  : new u16(16),
        ordr : [ 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 ],
        of0  : [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 999, 999, 999],
        exb  : [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,  4,  5,  5,  5,  5,  0,  0,  0,  0],
        ldef : new u16(32),
        df0  : [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577, 65535, 65535],
        dxb  : [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5,  6,  6,  7,  7,  8,  8,   9,   9,  10,  10,  11,  11,  12,   12,   13,   13,     0,     0],
        ddef : new u32(32),
        flmap: new u16(  512),  fltree: [] as number[],
        fdmap: new u16(   32),  fdtree: [] as number[],
        lmap : new u16(32768),  ltree : [] as number[],  ttree:[] as number[],
        dmap : new u16(32768),  dtree : [] as number[],
        imap : new u16(  512),  itree : [] as number[],
        // rev9 : new u16(  512)
        rev15: new u16(1 << 15),
        lhst : new u32(286), dhst : new u32( 30), ihst : new u32(19),
        lits : new u32(15000),
        strt : new u16(1 << 16),
        prev : new u16(1 << 15)
    };
})();

(function(){
    const len = 1 << 15;
    for(let i = 0; i < len; i++) {
        let x = i;
        x = (((x & 0xaaaaaaaa) >>> 1) | ((x & 0x55555555) << 1));
        x = (((x & 0xcccccccc) >>> 2) | ((x & 0x33333333) << 2));
        x = (((x & 0xf0f0f0f0) >>> 4) | ((x & 0x0f0f0f0f) << 4));
        x = (((x & 0xff00ff00) >>> 8) | ((x & 0x00ff00ff) << 8));
        U.rev15[i] = (((x >>> 16) | (x << 16))) >>> 17;
    }

    function pushV(tgt: number[], n: number, sv: number) {
        while(n-- !== 0) tgt.push(0, sv);
    }

    for(let i = 0; i < 32; i++) {
        U.ldef[i] = (U.of0[i] << 3) | U.exb[i];
        U.ddef[i] = (U.df0[i] << 4) | U.dxb[i];
    }

    pushV(U.fltree, 144, 8);
    pushV(U.fltree, 255 - 143, 9);
    pushV(U.fltree, 279 - 255, 7);
    pushV(U.fltree, 287 - 279, 8);
    /*
    var i = 0;
    for(; i<=143; i++) U.fltree.push(0,8);
    for(; i<=255; i++) U.fltree.push(0,9);
    for(; i<=279; i++) U.fltree.push(0,7);
    for(; i<=287; i++) U.fltree.push(0,8);
    */
    makeCodes(U.fltree, 9);
    codes2map(U.fltree, 9, U.flmap);
    revCodes (U.fltree, 9);

    pushV(U.fdtree, 32, 5);
    // for(i=0;i<32; i++) U.fdtree.push(0,5);
    makeCodes(U.fdtree, 5);
    codes2map(U.fdtree, 5, U.fdmap);
    revCodes (U.fdtree, 5);

    pushV(U.itree, 19, 0);  pushV(U.ltree, 286, 0);  pushV(U.dtree, 30, 0);  pushV(U.ttree, 320, 0);
    /*
    for(var i=0; i< 19; i++) U.itree.push(0,0);
    for(var i=0; i<286; i++) U.ltree.push(0,0);
    for(var i=0; i< 30; i++) U.dtree.push(0,0);
    for(var i=0; i<320; i++) U.ttree.push(0,0);
    */
})();

export function codes2map(tree: number[], MAX_BITS: number, map: Uint16Array) {
    const max_code = tree.length;
    const r15 = U.rev15;
    for(let i = 0; i < max_code; i += 2) {
        if(tree[i + 1] !== 0)  {
            const lit = i >> 1;
            const cl = tree[i + 1], val = (lit << 4) | cl; // :  (0x8000 | (U.of0[lit-257]<<7) | (U.exb[lit-257]<<4) | cl);
            const rest = (MAX_BITS - cl);
            let i0 = tree[i] << rest;
            const i1 = i0 + (1 << rest);
            // tree[i]=r15[i0]>>>(15-MAX_BITS);
            while(i0 !== i1) {
                const p0 = r15[i0] >>> (15 - MAX_BITS);
                map[p0] = val;  i0++;
            }
        }
    }
}

export function makeCodes(tree: number[], MAX_BITS: number) {  // code, length
    const max_code = tree.length;

    const bl_count = U.bl_count;
    for(let i = 0; i <= MAX_BITS; i++) bl_count[i] = 0;
    for(let i = 1; i < max_code; i += 2) bl_count[tree[i]]++;

    const next_code = U.next_code;	// smallest code for each length

    let code = 0;
    bl_count[0] = 0;
    for (let bits = 1; bits <= MAX_BITS; bits++) {
        code = (code + bl_count[bits - 1]) << 1;
        next_code[bits] = code;
    }

    for (let n = 0; n < max_code; n += 2) {
        const len = tree[n + 1];
        if (len !== 0) {
            tree[n] = next_code[len];
            next_code[len]++;
        }
    }
}

export function revCodes(tree: number[], MAX_BITS: number) {
    const r15 = U.rev15, imb = 15 - MAX_BITS;
    for(let i = 0; i < tree.length; i += 2) {
        const i0 = (tree[i] << (MAX_BITS - tree[i + 1]));
        tree[i] = r15[i0] >>> imb;
    }
}