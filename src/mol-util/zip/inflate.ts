/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 */

import { NumberArray } from '../type-helpers';
import { U, makeCodes, codes2map } from './util';
import { RuntimeContext } from '../../mol-task';

function InflateContext(data: Uint8Array, buf?: Uint8Array) {
    const noBuf = buf === undefined;
    if(buf === undefined) buf = new Uint8Array((data.length >>> 2) << 3);
    return {
        data,
        buf,
        noBuf,
        BFINAL: 0,
        off: 0,
        pos: 0
    };
}
type InflateContext = ReturnType<typeof InflateContext>

function inflateBlocks(ctx: InflateContext, count: number) {
    const { data, noBuf } = ctx;
    let { buf, BFINAL, off, pos } = ctx;

    let iBlock = 0;

    while(BFINAL === 0 && iBlock < count) {
        let lmap, dmap;
        let ML = 0, MD = 0;

        BFINAL = _bitsF(data, pos, 1);
        iBlock += 1;
        const BTYPE = _bitsF(data, pos + 1, 2);
        pos += 3;

        if(BTYPE === 0) {
            // uncompressed block
            if((pos & 7) !== 0) pos += 8 - (pos & 7);
            const p8 = (pos >>> 3) + 4;
            const len = data[p8 - 4] | (data[p8 - 3] << 8);
            if(noBuf) buf = _check(buf, off + len);
            buf.set(new Uint8Array(data.buffer, data.byteOffset + p8, len), off);
            pos = ((p8 + len) << 3);
            off += len;
            continue;
        }

        // grow output buffer if not provided
        if(noBuf) buf = _check(buf, off + (1 << 17));

        if(BTYPE === 1) {
            // block compressed with fixed Huffman codes
            lmap = U.flmap;
            dmap = U.fdmap;
            ML = (1 << 9) - 1;
            MD = (1 << 5) - 1;
        } else if(BTYPE === 2) {
            // block compressed with dynamic Huffman codes
            const HLIT = _bitsE(data, pos, 5) + 257;
            const HDIST = _bitsE(data, pos + 5, 5) + 1;
            const HCLEN = _bitsE(data, pos + 10, 4) + 4;
            pos += 14;

            for(let i = 0; i < 38; i += 2) {
                U.itree[i] = 0;
                U.itree[i + 1] = 0;
            }
            let tl = 1;
            for(let i = 0; i < HCLEN; i++) {
                const l = _bitsE(data, pos + i * 3, 3);
                U.itree[(U.ordr[i] << 1) + 1] = l;
                if(l > tl) tl = l;
            }
            pos += 3 * HCLEN;
            makeCodes(U.itree, tl);
            codes2map(U.itree, tl, U.imap);

            lmap = U.lmap;  dmap = U.dmap;

            pos = _decodeTiny(U.imap, (1 << tl) - 1, HLIT + HDIST, data, pos, U.ttree);
            const mx0 = _copyOut(U.ttree,    0, HLIT, U.ltree);
            ML = (1 << mx0) - 1;
            const mx1 = _copyOut(U.ttree, HLIT, HDIST, U.dtree);
            MD = (1 << mx1) - 1;

            makeCodes(U.ltree, mx0);
            codes2map(U.ltree, mx0, lmap);

            makeCodes(U.dtree, mx1);
            codes2map(U.dtree, mx1, dmap);
        } else {
            throw new Error(`unknown BTYPE ${BTYPE}`);
        }

        while(true) {
            const code = lmap[_get17(data, pos) & ML];
            pos += code & 15;
            const lit = code >>> 4;
            if((lit >>> 8) === 0) {
                buf[off++] = lit;
            } else if(lit === 256) {
                break;
            } else {
                let end = off + lit - 254;
                if(lit > 264) {
                    const ebs = U.ldef[lit - 257];
                    end = off + (ebs >>> 3) + _bitsE(data, pos, ebs & 7);
                    pos += ebs & 7;
                }

                const dcode = dmap[_get17(data, pos) & MD];
                pos += dcode & 15;
                const dlit = dcode >>> 4;
                const dbs = U.ddef[dlit];
                const dst = (dbs >>> 4) + _bitsF(data, pos, dbs & 15);
                pos += dbs & 15;

                if(noBuf) buf = _check(buf, off + (1 << 17));
                while(off < end) {
                    buf[off] = buf[off++ - dst];
                    buf[off] = buf[off++ - dst];
                    buf[off] = buf[off++ - dst];
                    buf[off] = buf[off++ - dst];
                }
                off = end;
            }
        }
    }

    ctx.buf = buf;
    ctx.BFINAL = BFINAL;
    ctx.off = off;
    ctx.pos = pos;
}

// https://tools.ietf.org/html/rfc1951
export async function _inflate(runtime: RuntimeContext, data: Uint8Array, buf?: Uint8Array) {
    if(data[0] === 3 && data[1] === 0) return (buf ? buf : new Uint8Array(0));

    const ctx = InflateContext(data, buf);
    while(ctx.BFINAL === 0) {
        if (runtime.shouldUpdate) {
            await runtime.update({ message: 'Inflating blocks...', current: ctx.pos, max: data.length });
        }
        inflateBlocks(ctx, 100);
    }
    return ctx.buf.length === ctx.off ? ctx.buf : ctx.buf.slice(0, ctx.off);
}

function _check(buf: Uint8Array, len: number) {
    const bl = buf.length;
    if(len <= bl) return buf;
    const nbuf = new Uint8Array(Math.max(bl << 1, len));
    nbuf.set(buf, 0);
    return nbuf;
}

function _decodeTiny(lmap: NumberArray, LL: number, len: number, data: Uint8Array, pos: number, tree: number[]) {
    let i = 0;
    while(i < len) {
        const code = lmap[_get17(data, pos) & LL];
        pos += code & 15;
        const lit = code >>> 4;
        if(lit <= 15) {
            tree[i] = lit;
            i++;
        } else {
            let ll = 0, n = 0;
            if(lit === 16) {
                n = (3  + _bitsE(data, pos, 2));
                pos += 2;
                ll = tree[i - 1];
            } else if(lit === 17) {
                n = (3  + _bitsE(data, pos, 3));
                pos += 3;
            } else if(lit === 18) {
                n = (11 + _bitsE(data, pos, 7));
                pos += 7;
            }
            const ni = i + n;
            while(i < ni) {
                tree[i] = ll;
                i++;
            }
        }
    }
    return pos;
}

function _copyOut(src: number[], off: number, len: number, tree: number[]) {
    let mx = 0, i = 0;
    const tl = tree.length >>> 1;
    while(i < len) {
        let v = src[i + off];
        tree[(i << 1)] = 0;
        tree[(i << 1) + 1] = v;
        if(v > mx)mx = v;
        i++;
    }
    while(i < tl) {
        tree[(i << 1)] = 0;
        tree[(i << 1) + 1] = 0;
        i++;
    }
    return mx;
}

function _bitsE(dt: NumberArray, pos: number, length: number) {
    return ((dt[pos >>> 3] | (dt[(pos >>> 3) + 1] << 8)) >>> (pos & 7)) & ((1 << length) - 1);
}

function _bitsF(dt: NumberArray, pos: number, length: number) {
    return ((dt[pos >>> 3] | (dt[(pos >>> 3) + 1] << 8) | (dt[(pos >>> 3) + 2] << 16)) >>> (pos & 7)) & ((1 << length) - 1);
}

function _get17(dt: NumberArray, pos: number) {	// return at least 17 meaningful bytes
    return (dt[pos >>> 3] | (dt[(pos >>> 3) + 1] << 8) | (dt[(pos >>> 3) + 2] << 16) ) >>> (pos & 7);
}