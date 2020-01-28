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

export function _inflate(data: Uint8Array, buf?: Uint8Array) {
    if(data[0] === 3 && data[1] === 0) return (buf ? buf : new Uint8Array(0));
    // var F=UZIP.F, bitsF = F._bitsF, bitsE = F._bitsE, decodeTiny = F._decodeTiny, makeCodes = F.makeCodes, codes2map=F.codes2map, get17 = F._get17;
    // var U = F.U;

    const noBuf = buf === undefined;
    if(buf === undefined) buf = new Uint8Array((data.length>>>2)<<3);

    let BFINAL=0, BTYPE=0, HLIT=0, HDIST=0, HCLEN=0, ML=0, MD=0;
    let off = 0, pos = 0;
    let lmap, dmap;

    while(BFINAL === 0) {
        BFINAL = _bitsF(data, pos  , 1);
        BTYPE  = _bitsF(data, pos+1, 2);
        pos+=3;

        if(BTYPE === 0) {
            if((pos&7) !== 0) pos+=8-(pos&7);
            const p8 = (pos>>>3)+4, len = data[p8-4]|(data[p8-3]<<8);  // console.log(len);//bitsF(data, pos, 16),
            if(noBuf) buf=_check(buf, off+len);
            buf.set(new Uint8Array(data.buffer, data.byteOffset+p8, len), off);
            // for(var i=0; i<len; i++) buf[off+i] = data[p8+i];
            // for(var i=0; i<len; i++) if(buf[off+i] != data[p8+i]) throw "e";
            pos = ((p8+len)<<3);  off+=len;  continue;
        }
        if(noBuf) buf=_check(buf, off+(1<<17));  // really not enough in many cases (but PNG and ZIP provide buffer in advance)
        if(BTYPE === 1) {
            lmap = U.flmap;
            dmap = U.fdmap;
            ML = (1<<9)-1;
            MD = (1<<5)-1;
        } else if(BTYPE === 2) {
            HLIT  = _bitsE(data, pos   , 5)+257;
            HDIST = _bitsE(data, pos+ 5, 5)+  1;
            HCLEN = _bitsE(data, pos+10, 4)+  4;  pos+=14;
            // const ppos = pos;
            for(let i=0; i<38; i+=2) {
                U.itree[i]=0;
                U.itree[i+1]=0;
            }
            let tl = 1;
            for(let i=0; i<HCLEN; i++) {
                const l=_bitsE(data, pos+i*3, 3);
                U.itree[(U.ordr[i]<<1)+1] = l;
                if(l>tl) tl = l;
            }
            pos+=3*HCLEN;  // console.log(itree);
            makeCodes(U.itree, tl);
            codes2map(U.itree, tl, U.imap);

            lmap = U.lmap;  dmap = U.dmap;

            pos = _decodeTiny(U.imap, (1<<tl)-1, HLIT+HDIST, data, pos, U.ttree);
            const mx0 = _copyOut(U.ttree,    0, HLIT , U.ltree);  ML = (1<<mx0)-1;
            const mx1 = _copyOut(U.ttree, HLIT, HDIST, U.dtree);  MD = (1<<mx1)-1;

            // var ml = decodeTiny(U.imap, (1<<tl)-1, HLIT , data, pos, U.ltree); ML = (1<<(ml>>>24))-1;  pos+=(ml&0xffffff);
            makeCodes(U.ltree, mx0);
            codes2map(U.ltree, mx0, lmap);

            // var md = decodeTiny(U.imap, (1<<tl)-1, HDIST, data, pos, U.dtree); MD = (1<<(md>>>24))-1;  pos+=(md&0xffffff);
            makeCodes(U.dtree, mx1);
            codes2map(U.dtree, mx1, dmap);
        } else {
            throw new Error(`unknown BTYPE ${BTYPE}`)
        }

        // var ooff=off, opos=pos;
        while(true) {
            const code = lmap[_get17(data, pos) & ML];
            pos += code&15;
            const lit = code >>> 4;  // U.lhst[lit]++;
            if((lit >>> 8) === 0) {
                buf[off++] = lit;
            } else if(lit === 256) {
                break;
            } else {
                let end = off+lit-254;
                if(lit > 264) {
                    const ebs = U.ldef[lit-257];
                    end = off + (ebs>>>3) + _bitsE(data, pos, ebs&7);
                    pos += ebs&7;
                }
                // UZIP.F.dst[end-off]++;

                const dcode = dmap[_get17(data, pos) & MD];  pos += dcode&15;
                const dlit = dcode>>>4;
                const dbs = U.ddef[dlit], dst = (dbs>>>4) + _bitsF(data, pos, dbs&15);  pos += dbs&15;

                // var o0 = off-dst, stp = Math.min(end-off, dst);
                // if(stp>20) while(off<end) {  buf.copyWithin(off, o0, o0+stp);  off+=stp;  }  else
                // if(end-dst<=off) buf.copyWithin(off, off-dst, end-dst);  else
                // if(dst==1) buf.fill(buf[off-1], off, end);  else
                if(noBuf) buf = _check(buf, off+(1<<17));
                while(off<end) {
                    buf[off]=buf[off++-dst];
                    buf[off]=buf[off++-dst];
                    buf[off]=buf[off++-dst];
                    buf[off]=buf[off++-dst];
                }
                off=end;
                // while(off!=end) {  buf[off]=buf[off++-dst];  }
            }
        }
        // console.log(off-ooff, (pos-opos)>>>3);
    }
    // console.log(UZIP.F.dst);
    // console.log(tlen, dlen, off-tlen+tcnt);
    return buf.length === off ? buf : buf.slice(0, off);
}

function _check(buf: Uint8Array, len: number) {
    const bl = buf.length;
    if(len <= bl) return buf;
    const nbuf = new Uint8Array(Math.max(bl << 1, len));
    nbuf.set(buf, 0);
    // for(var i=0; i<bl; i+=4) {  nbuf[i]=buf[i];  nbuf[i+1]=buf[i+1];  nbuf[i+2]=buf[i+2];  nbuf[i+3]=buf[i+3];  }
    return nbuf;
}

function _decodeTiny(lmap: NumberArray, LL: number, len: number, data: Uint8Array, pos: number, tree: number[]) {
    let i = 0;
    while(i<len) {
        const code = lmap[_get17(data, pos)&LL];
        pos += code&15;
        const lit = code>>>4;
        if(lit<=15) {
            tree[i]=lit;
            i++;
        } else {
            let ll = 0, n = 0;
            if(lit === 16) {
                n = (3  + _bitsE(data, pos, 2));  pos += 2;  ll = tree[i-1];
            } else if(lit === 17) {
                n = (3  + _bitsE(data, pos, 3));  pos += 3;
            } else if(lit === 18) {
                n = (11 + _bitsE(data, pos, 7));  pos += 7;
            }
            const ni = i+n;
            while(i<ni) {
                tree[i]=ll;
                i++;
            }
        }
    }
    return pos;
}

function _copyOut(src: number[], off: number, len: number, tree: number[]) {
    let mx=0, i=0, tl=tree.length>>>1;
    while(i<len) {
        let v=src[i+off];
        tree[(i<<1)]=0;
        tree[(i<<1)+1]=v;
        if(v>mx)mx=v;
        i++;
    }
    while(i<tl ) {
        tree[(i<<1)]=0;
        tree[(i<<1)+1]=0;
        i++;
    }
    return mx;
}

function _bitsE(dt: NumberArray, pos: number, length: number) {
    return ((dt[pos>>>3] | (dt[(pos>>>3)+1]<<8))>>>(pos&7))&((1<<length)-1);
}

function _bitsF(dt: NumberArray, pos: number, length: number) {
    return ((dt[pos>>>3] | (dt[(pos>>>3)+1]<<8) | (dt[(pos>>>3)+2]<<16))>>>(pos&7))&((1<<length)-1);
}

function _get17(dt: NumberArray, pos: number) {	// return at least 17 meaningful bytes
    return (dt[pos>>>3] | (dt[(pos>>>3)+1]<<8) | (dt[(pos>>>3)+2]<<16) )>>>(pos&7);
}