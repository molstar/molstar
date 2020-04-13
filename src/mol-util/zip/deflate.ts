/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 */

import { NumberArray } from '../type-helpers';
import { _hufTree } from './huffman';
import { U, revCodes, makeCodes } from './util';

export function _deflateRaw(data: Uint8Array, out: Uint8Array, opos: number, lvl: number) {
    const opts = [
        /*
            ush good_length; /* reduce lazy search above this match length
            ush max_lazy;    /* do not perform lazy search above this match length
            ush nice_length; /* quit search above this match length
        */
        /*      good lazy nice chain */
        /* 0 */ [ 0,   0,   0,    0, 0], /* store only */
        /* 1 */ [ 4,   4,   8,    4, 0], /* max speed, no lazy matches */
        /* 2 */ [ 4,   5,  16,    8, 0],
        /* 3 */ [ 4,   6,  16,   16, 0],

        /* 4 */ [ 4,  10,  16,   32, 0], /* lazy matches */
        /* 5 */ [ 8,  16,  32,   32, 0],
        /* 6 */ [ 8,  16, 128,  128, 0],
        /* 7 */ [ 8,  32, 128,  256, 0],
        /* 8 */ [32, 128, 258, 1024, 1],
        /* 9 */ [32, 258, 258, 4096, 1] /* max compression */
    ];

    const opt = opts[lvl];

    let i = 0, pos = opos << 3, cvrd = 0;
    const dlen = data.length;

    if(lvl === 0) {
        while(i < dlen) {
            const len = Math.min(0xffff, dlen - i);
            _putsE(out, pos, (i + len === dlen ? 1 : 0));
            pos = _copyExact(data, i, len, out, pos + 8);
            i += len;
        }
        return pos >>> 3;
    }

    const { lits, strt, prev } = U;
    let li = 0, lc = 0, bs = 0, ebits = 0, c = 0, nc = 0;  // last_item, literal_count, block_start
    if(dlen > 2) {
        nc = _hash(data, 0);
        strt[nc] = 0;
    }

    // let nmch = 0
    // let nmci = 0

    for(i = 0; i < dlen; i++)  {
        c = nc;
        //*
        if(i + 1 < dlen - 2) {
            nc = _hash(data, i + 1);
            const ii = ((i + 1) & 0x7fff);
            prev[ii] = strt[nc];
            strt[nc] = ii;
        } // */
        if(cvrd <= i) {
            if((li > 14000 || lc > 26697) && (dlen - i) > 100) {
                if(cvrd < i) {
                    lits[li] = i - cvrd;
                    li += 2;
                    cvrd = i;
                }
                pos = _writeBlock(((i === dlen - 1) || (cvrd === dlen)) ? 1 : 0, lits, li, ebits, data, bs, i - bs, out, pos);
                li = lc = ebits = 0;
                bs = i;
            }

            let mch = 0;
            // if(nmci==i) mch= nmch;  else
            if(i < dlen - 2) {
                mch = _bestMatch(data, i, prev, c, Math.min(opt[2], dlen - i), opt[3]);
            }
            /*
            if(mch!=0 && opt[4]==1 && (mch>>>16)<opt[1] && i+1<dlen-2) {
                nmch = UZIP.F._bestMatch(data, i+1, prev, nc, opt[2], opt[3]);  nmci=i+1;
                //var mch2 = UZIP.F._bestMatch(data, i+2, prev, nnc);  //nmci=i+1;
                if((nmch>>>16)>(mch>>>16)) mch=0;
            }//*/
            // const len = mch>>>16, dst = mch & 0xffff;  // if(i-dst<0) throw "e";
            if(mch !== 0) {
                const len = mch >>> 16, dst = mch & 0xffff;  // if(i-dst<0) throw "e";
                const lgi = _goodIndex(len, U.of0);  U.lhst[257 + lgi]++;
                const dgi = _goodIndex(dst, U.df0);  U.dhst[    dgi]++;  ebits += U.exb[lgi] + U.dxb[dgi];
                lits[li] = (len << 23) | (i - cvrd);  lits[li + 1] = (dst << 16) | (lgi << 8) | dgi;  li += 2;
                cvrd = i + len;
            } else {
                U.lhst[data[i]]++;
            }
            lc++;
        }
    }
    if(bs !== i || data.length === 0) {
        if(cvrd < i) {
            lits[li] = i - cvrd;
            li += 2;
            cvrd = i;
        }
        pos = _writeBlock(1, lits, li, ebits, data, bs, i - bs, out, pos);
        li = 0;
        lc = 0;
        li = lc = ebits = 0;
        bs = i;
    }
    while((pos & 7) !== 0) pos++;
    return pos >>> 3;
}

function _bestMatch(data: Uint8Array, i: number, prev: Uint16Array, c: number, nice: number, chain: number) {
    let ci = (i & 0x7fff), pi = prev[ci];
    // console.log("----", i);
    let dif = ((ci - pi + (1 << 15)) & 0x7fff);
    if(pi === ci || c !== _hash(data, i - dif)) return 0;
    let tl = 0, td = 0;  // top length, top distance
    const dlim = Math.min(0x7fff, i);
    while(dif <= dlim && --chain !== 0 && pi !== ci /* && c==UZIP.F._hash(data,i-dif)*/) {
        if(tl === 0 || (data[i + tl] === data[i + tl - dif])) {
            let cl = _howLong(data, i, dif);
            if(cl > tl) {
                tl = cl;  td = dif;  if(tl >= nice) break;    //*
                if(dif + 2 < cl) cl = dif + 2;
                let maxd = 0; // pi does not point to the start of the word
                for(let j = 0; j < cl - 2; j++) {
                    const ei =  (i - dif + j + (1 << 15)) & 0x7fff;
                    const li = prev[ei];
                    const curd = (ei - li + (1 << 15)) & 0x7fff;
                    if(curd > maxd) {  maxd = curd;  pi = ei; }
                }
            }
        }

        ci = pi;  pi = prev[ci];
        dif += ((ci - pi + (1 << 15)) & 0x7fff);
    }
    return (tl << 16) | td;
}

function _howLong(data: Uint8Array, i: number, dif: number) {
    if(data[i] !== data[i - dif] || data[i + 1] !== data[i + 1 - dif] || data[i + 2] !== data[i + 2 - dif]) return 0;
    const oi = i, l = Math.min(data.length, i + 258);
    i += 3;
    // while(i+4<l && data[i]==data[i-dif] && data[i+1]==data[i+1-dif] && data[i+2]==data[i+2-dif] && data[i+3]==data[i+3-dif]) i+=4;
    while(i < l && data[i] === data[i - dif]) i++;
    return i - oi;
}

function _hash(data: Uint8Array, i: number) {
    return (((data[i] << 8) | data[i + 1]) + (data[i + 2] << 4)) & 0xffff;
    // var hash_shift = 0, hash_mask = 255;
    // var h = data[i+1] % 251;
    // h = (((h << 8) + data[i+2]) % 251);
    // h = (((h << 8) + data[i+2]) % 251);
    // h = ((h<<hash_shift) ^ (c) ) & hash_mask;
    // return h | (data[i]<<8);
    // return (data[i] | (data[i+1]<<8));
}

function _writeBlock(BFINAL: number, lits: Uint32Array, li: number, ebits: number, data: Uint8Array, o0: number, l0: number, out: Uint8Array, pos: number) {
    U.lhst[256]++;
    const [ ML, MD, MH, numl, numd, numh, lset, dset ] = getTrees();

    const cstSize = (((pos + 3) & 7) === 0 ? 0 : 8 - ((pos + 3) & 7)) + 32 + (l0 << 3);
    const fxdSize = ebits + contSize(U.fltree, U.lhst) + contSize(U.fdtree, U.dhst);
    let dynSize = ebits + contSize(U.ltree, U.lhst) + contSize(U.dtree, U.dhst);
    dynSize += 14 + 3 * numh + contSize(U.itree, U.ihst) + (U.ihst[16] * 2 + U.ihst[17] * 3 + U.ihst[18] * 7);

    for(let j = 0; j < 286; j++) U.lhst[j] = 0;
    for(let j = 0; j < 30; j++) U.dhst[j] = 0;
    for(let j = 0; j < 19; j++) U.ihst[j] = 0;

    const BTYPE = (cstSize < fxdSize && cstSize < dynSize) ? 0 : ( fxdSize < dynSize ? 1 : 2 );
    _putsF(out, pos, BFINAL);
    _putsF(out, pos + 1, BTYPE);
    pos += 3;

    // let opos = pos;
    if(BTYPE === 0) {
        while((pos & 7) !== 0) pos++;
        pos = _copyExact(data, o0, l0, out, pos);
    } else {
        let ltree: number[], dtree: number[];
        if(BTYPE === 1) {
            ltree = U.fltree;  dtree = U.fdtree;
        } else if(BTYPE === 2) {
            makeCodes(U.ltree, ML);  revCodes(U.ltree, ML);
            makeCodes(U.dtree, MD);  revCodes(U.dtree, MD);
            makeCodes(U.itree, MH);  revCodes(U.itree, MH);

            ltree = U.ltree;  dtree = U.dtree;

            _putsE(out, pos, numl - 257);  pos += 5;  // 286
            _putsE(out, pos, numd -  1);  pos += 5;  // 30
            _putsE(out, pos, numh -  4);  pos += 4;  // 19

            for(let i = 0; i < numh; i++) _putsE(out, pos + i * 3, U.itree[(U.ordr[i] << 1) + 1]);
            pos += 3 * numh;
            pos = _codeTiny(lset, U.itree, out, pos);
            pos = _codeTiny(dset, U.itree, out, pos);
        } else {
            throw new Error(`unknown BTYPE ${BTYPE}`);
        }

        let off = o0;
        for(let si = 0; si < li; si += 2) {
            const qb = lits[si], len = (qb >>> 23), end = off + (qb & ((1 << 23) - 1));
            while(off < end) pos = _writeLit(data[off++], ltree, out, pos);

            if(len !== 0) {
                const qc = lits[si + 1], dst = (qc >> 16), lgi = (qc >> 8) & 255, dgi = (qc & 255);
                pos = _writeLit(257 + lgi, ltree, out, pos);
                _putsE(out, pos, len - U.of0[lgi]);  pos += U.exb[lgi];

                pos = _writeLit(dgi, dtree, out, pos);
                _putsF(out, pos, dst - U.df0[dgi]);  pos += U.dxb[dgi];  off += len;
            }
        }
        pos = _writeLit(256, ltree, out, pos);
    }
    // console.log(pos-opos, fxdSize, dynSize, cstSize);
    return pos;
}

function _copyExact(data: Uint8Array, off: number, len: number, out: Uint8Array, pos: number) {
    let p8 = (pos >>> 3);
    out[p8] = (len);
    out[p8 + 1] = (len >>> 8);
    out[p8 + 2] = 255 - out[p8];
    out[p8 + 3] = 255 - out[p8 + 1];
    p8 += 4;
    out.set(new Uint8Array(data.buffer, off, len), p8);
    // for(var i=0; i<len; i++) out[p8+i]=data[off+i];
    return pos + ((len + 4) << 3);
}


/*
    Interesting facts:
    - decompressed block can have bytes, which do not occur in a Huffman tree (copied from the previous block by reference)
*/

function getTrees() {
    const ML = _hufTree(U.lhst, U.ltree, 15);
    const MD = _hufTree(U.dhst, U.dtree, 15);
    const lset: number[] = [];
    const numl = _lenCodes(U.ltree, lset);
    const dset: number[] = [];
    const numd = _lenCodes(U.dtree, dset);
    for(let i = 0; i < lset.length; i += 2) U.ihst[lset[i]]++;
    for(let i = 0; i < dset.length; i += 2) U.ihst[dset[i]]++;
    const MH = _hufTree(U.ihst, U.itree,  7);
    let numh = 19;
    while(numh > 4 && U.itree[(U.ordr[numh - 1] << 1) + 1] === 0) numh--;
    return [ML, MD, MH, numl, numd, numh, lset, dset] as const;
}

function contSize(tree: number[], hst: NumberArray) {
    let s = 0;
    for(let i = 0; i < hst.length; i++) s += hst[i] * tree[(i << 1) + 1];
    return s;
}

function _codeTiny(set: number[], tree: number[], out: Uint8Array, pos: number) {
    for(let i = 0; i < set.length; i += 2) {
        const l = set[i], rst = set[i + 1];  // console.log(l, pos, tree[(l<<1)+1]);
        pos = _writeLit(l, tree, out, pos);
        const rsl = l === 16 ? 2 : (l === 17 ? 3 : 7);
        if(l > 15) {
            _putsE(out, pos, rst);
            pos += rsl;
        }
    }
    return pos;
}


function _lenCodes(tree: number[], set: number[]) {
    let len = tree.length;
    while(len !== 2 && tree[len - 1] === 0) len -= 2;  // when no distances, keep one code with length 0
    for(let i = 0; i < len; i += 2) {
        const l = tree[i + 1], nxt = (i + 3 < len ? tree[i + 3] : -1),  nnxt = (i + 5 < len ? tree[i + 5] : -1),  prv = (i === 0 ? -1 : tree[i - 1]);
        if(l === 0 && nxt === l && nnxt === l) {
            let lz = i + 5;
            while(lz + 2 < len && tree[lz + 2] === l) lz += 2;
            const zc = Math.min((lz + 1 - i) >>> 1, 138);
            if(zc < 11) set.push(17, zc - 3);
            else set.push(18, zc - 11);
            i += zc * 2 - 2;
        } else if(l === prv && nxt === l && nnxt === l) {
            let lz = i + 5;
            while(lz + 2 < len && tree[lz + 2] === l) lz += 2;
            const zc = Math.min((lz + 1 - i) >>> 1, 6);
            set.push(16, zc - 3);
            i += zc * 2 - 2;
        } else {
            set.push(l, 0);
        }
    }
    return len >>> 1;
}

function _goodIndex(v: number, arr: number[]) {
    let i = 0;
    if(arr[i | 16] <= v) i |= 16;
    if(arr[i | 8] <= v) i |= 8;
    if(arr[i | 4] <= v) i |= 4;
    if(arr[i | 2] <= v) i |= 2;
    if(arr[i | 1] <= v) i |= 1;
    return i;
}

function _writeLit(ch: number, ltree: number[], out: Uint8Array, pos: number) {
    _putsF(out, pos, ltree[ch << 1]);
    return pos + ltree[(ch << 1) + 1];
}

function _putsE(dt: NumberArray, pos: number, val: number) {
    val = val << (pos & 7);
    const o = (pos >>> 3);
    dt[o] |= val;
    dt[o + 1] |= (val >>> 8);
}

function _putsF(dt: NumberArray, pos: number, val: number) {
    val = val << (pos & 7);
    const o = (pos >>> 3);
    dt[o] |= val;
    dt[o + 1] |= (val >>> 8);
    dt[o + 2] |= (val >>> 16);
}
