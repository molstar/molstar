import * as B from 'benchmark';
import { Tuple, Segmentation, OrderedSet as OrdSet } from '../mol-data/int';
// import { ElementSet } from 'mol-model/structure'

// export namespace Iteration {
//     const U = 1000, V = 2500;

//     const control: number[] = [];
//     const sets = Object.create(null);
//     for (let i = 1; i < U; i++) {
//         const set = [];
//         for (let j = 1; j < V; j++) {
//             control[control.length] = j * j + 1;
//             set[set.length] = j * j + 1;
//         }
//         sets[i * i] = OrdSet.ofSortedArray(set);
//     }
//     const ms = ElementSet.create(sets);

//     export function native() {
//         let s = 0;
//         for (let i = 0, _i = control.length; i < _i; i++) s += control[i];
//         return s;
//     }

//     export function iterators() {
//         let s = 0;
//         const it = ElementSet.atoms(ms);
//         while (it.hasNext) {
//             const v = it.move();
//             s += Tuple.snd(v);
//         }
//         return s;
//     }

//     export function elementAt() {
//         let s = 0;
//         for (let i = 0, _i = ElementSet.atomCount(ms); i < _i; i++) s += Tuple.snd(ElementSet.atomGetAt(ms, i));
//         return s;
//     }

//     export function manual() {
//         let s = 0;
//         const keys = ElementSet.unitIds(ms);
//         for (let i = 0, _i = OrdSet.size(keys); i < _i; i++) {
//             const set = ElementSet.unitGetById(ms, OrdSet.getAt(keys, i));
//             for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
//                 s += OrdSet.getAt(set, j);
//             }
//         }
//         return s;
//     }

//     export function manual1() {
//         let s = 0;
//         for (let i = 0, _i = ElementSet.unitCount(ms); i < _i; i++) {
//             const set = ElementSet.unitGetByIndex(ms, i);
//             for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
//                 s += OrdSet.getAt(set, j);
//             }
//         }
//         return s;
//     }

//     export function run() {
//         const suite = new B.Suite();

//         // const values: number[] = [];
//         // for (let i = 0; i < 1000000; i++) values[i] = (Math.random() * 1000) | 0;

//         console.log(Iteration.native(), Iteration.iterators(), Iteration.elementAt(), Iteration.manual(), Iteration.manual1());

//         suite
//             .add('native', () => Iteration.native())
//             .add('iterators', () => Iteration.iterators())
//             .add('manual', () => Iteration.manual())
//             .add('manual1', () => Iteration.manual1())
//             .add('el at', () => Iteration.elementAt())
//             .on('cycle', (e: any) => console.log(String(e.target)))
//             .run();
//     }
// }

export namespace Union {
    function createArray(n: number) {
        const data = new Int32Array(n);
        let c = (Math.random() * 100) | 0;
        for (let i = 0; i < n; i++) {
            data[i] = c;
            c += 1 + (Math.random() * 100) | 0;
        }
        return data;
    }

    function rangeArray(o: number, n: number) {
        const data = new Int32Array(n);
        for (let i = 0; i < n; i++) {
            data[i] = o + i;
        }
        return data;
    }

    function createData(array: ArrayLike<number>) {
        const obj = Object.create(null);
        const set = new Set<number>();
        for (let i = 0; i < array.length; i++) {
            const a = array[i];
            obj[a] = 1;
            set.add(a);
        }

        return { ordSet: OrdSet.ofSortedArray(array), obj, set };
    }

    function unionOS(a: OrdSet, b: OrdSet) {
        return OrdSet.union(a, b);
    }

    function intOS(a: OrdSet, b: OrdSet) {
        return OrdSet.intersect(a, b);
    }

    function unionO(a: object, b: object) {
        const ret = Object.create(null);
        for (const k of Object.keys(a)) ret[k] = 1;
        for (const k of Object.keys(b)) ret[k] = 1;
        return ret;
    }

    function intO(a: object, b: object) {
        const ret = Object.create(null);
        for (const k of Object.keys(a)) if ((b as any)[k]) ret[k] = 1;
        return ret;
    }

    function _setAdd(this: Set<number>, x: number) { this.add(x); }
    function unionS(a: Set<number>, b: Set<number>) {
        const ret = new Set<number>();
        a.forEach(_setAdd, ret);
        b.forEach(_setAdd, ret);
        return ret;
    }

    function _setInt(this: { set: Set<number>, other: Set<number> }, x: number) { if (this.other.has(x)) this.set.add(x); }
    function intS(a: Set<number>, b: Set<number>) {
        if (a.size < b.size) {
            const ctx = { set: new Set<number>(), other: b };
            a.forEach(_setInt, ctx);
            return ctx.set;
        } else {
            const ctx = { set: new Set<number>(), other: a };
            b.forEach(_setInt, ctx);
            return ctx.set;
        }
    }

    export function run() {
        const suite = new B.Suite();

        const { ordSet: osA, set: sA, obj: oA } = createData(createArray(1000));
        const { ordSet: osB, set: sB, obj: oB } = createData(createArray(1000));

        console.log(OrdSet.size(unionOS(osA, osB)), Object.keys(unionO(oA, oB)).length, unionS(sA, sB).size);
        console.log(OrdSet.size(intOS(osA, osB)), Object.keys(intO(oA, oB)).length, intS(sA, sB).size);

        suite
            .add('u ord set', () => unionOS(osA, osB))
            .add('u obj', () => unionO(oA, oB))
            .add('u ES6 set', () => unionS(sA, sB))
            .add('i ord set', () => intOS(osA, osB))
            .add('i obj', () => intO(oA, oB))
            .add('i ES6 set', () => intS(sA, sB))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }

    export function runR() {
        const suite = new B.Suite();

        const rangeA = rangeArray(145, 1000);
        const rangeB = rangeArray(456, 1000);

        const { ordSet: osA, set: sA, obj: oA } = createData(rangeA);
        const { ordSet: osB, set: sB, obj: oB } = createData(rangeB);

        console.log(OrdSet.size(unionOS(osA, osB)), Object.keys(unionO(oA, oB)).length, unionS(sA, sB).size);
        console.log(OrdSet.size(intOS(osA, osB)), Object.keys(intO(oA, oB)).length, intS(sA, sB).size);

        suite
            .add('u ord set', () => unionOS(osA, osB))
            .add('u obj', () => unionO(oA, oB))
            .add('u ES6 set', () => unionS(sA, sB))
            .add('i ord set', () => intOS(osA, osB))
            .add('i obj', () => intO(oA, oB))
            .add('i ES6 set', () => intS(sA, sB))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

// export namespace Build {
//     function createSorted() {
//         const b = ElementSet.LinearBuilder(ElementSet.Empty);
//         for (let i = 0; i < 10; i++) {
//             for (let j = 0; j < 1000; j++) {
//                 b.add(i, j);
//             }
//         }
//         return b.getSet();
//     }

//     function createByUnit() {
//         const b = ElementSet.LinearBuilder(ElementSet.Empty);
//         for (let i = 0; i < 10; i++) {
//             b.beginUnit();
//             for (let j = 0; j < 1000; j++) {
//                 b.addToUnit(j);
//             }
//             b.commitUnit(i);
//         }
//         return b.getSet();
//     }


//     export function run() {
//         const suite = new B.Suite();
//         suite
//             .add('create sorted', () => createSorted())
//             .add('create by unit', () => createByUnit())
//             .on('cycle', (e: any) => console.log(String(e.target)))
//             .run();
//     }
// }

export namespace Tuples {
    function createData(n: number) {
        const ret: Tuple[] = new Float64Array(n) as any;
        for (let i = 0; i < n; i++) {
            ret[i] = Tuple.create(i, i * i + 1);
        }
        return ret;
    }

    function sum1(data: ArrayLike<Tuple>) {
        let s = 0;
        for (let i = 0, _i = data.length; i < _i; i++) {
            s += Tuple.fst(data[i]) + Tuple.snd(data[i]);
        }
        return s;
    }

    function sum2(data: ArrayLike<Tuple>) {
        let s = 0;
        for (let i = 0, _i = data.length; i < _i; i++) {
            const t = data[i];
            s += Tuple.fst(t) + Tuple.snd(t);
        }
        return s;
    }

    export function run() {
        const suite = new B.Suite();
        const data = createData(10000);
        suite
            .add('sum fst/snd', () => sum1(data))
            .add('sum fst/snd 1', () => sum2(data))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

export function testSegments() {
    const data = OrdSet.ofSortedArray([4, 9, 10, 11, 14, 15, 16]);
    const segs = Segmentation.create([0, 4, 10, 12, 13, 15, 25]);
    const it = Segmentation.transientSegments(segs, data);

    while (it.hasNext) {
        const s = it.move();
        for (let j = s.start; j < s.end; j++) {
            console.log(`${s.index}: ${OrdSet.getAt(data, j)}`);
        }
    }
}

export namespace ObjectVsMap {
    function objCreate(keys: string[]) {
        const m = Object.create(null);
        m.x = 0;
        delete m.x;
        for (let i = 0, _i = keys.length; i < _i; i++) {
            m[keys[i]] = i * i;
        }
        return m;
    }

    function mapCreate(keys: string[]) {
        const m = new Map<string, number>();
        for (let i = 0, _i = keys.length; i < _i; i++) {
            m.set(keys[i], i * i);
        }
        return m;
    }

    function objQuery(keys: string[], n: number, obj: any) {
        let ret = 0;
        for (let i = 0; i < n; i++) ret += obj[keys[i % n]];
        return ret;
    }

    function mapQuery(keys: string[], n: number, map: Map<string, number>) {
        let ret = 0;
        for (let i = 0; i < n; i++) ret += map.get(keys[i % n])!;
        return ret;
    }

    export function run() {
        const suite = new B.Suite();
        const keys: string[] = [];
        for (let i = 0; i < 1000; i++) keys[i] = 'k' + i;

        const N = 100000;
        const obj = objCreate(keys);
        const map = mapCreate(keys);
        suite
            .add('c obj', () => objCreate(keys))
            .add('c map', () => mapCreate(keys))
            .add('q obj', () => objQuery(keys, N, obj))
            .add('q map', () => mapQuery(keys, N, map))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

export namespace IntVsStringIndices {
    type WithKeys<K> = { keys: K[], data: { [key: number]: number } }
    type MapWithKeys = { keys: number[], map: Map<number, number> }

    function createCacheKeys(n: number): WithKeys<number> {
        const data = Object.create(null), keys = [];
        data.__ = void 0;
        delete data.__;
        for (let i = 1; i <= n; i++) {
            const k = i * (i + 1);
            keys[keys.length] = k;
            data[k] = i + 1;
        }
        return { data, keys };
    }

    function createMapKeys(n: number): MapWithKeys {
        const map = new Map<number, number>(), keys = [];
        for (let i = 1; i <= n; i++) {
            const k = i * (i + 1);
            keys[keys.length] = k;
            map.set(k, i + 1);
        }
        return { map, keys };
    }

    export function createInt(n: number) {
        const ret = Object.create(null);
        ret.__ = void 0;
        delete ret.__;
        for (let i = 1; i <= n; i++) ret[i * (i + 1)] = i + 1;
        return ret;
    }

    export function createStr(n: number) {
        const ret = Object.create(null);
        ret.__ = void 0;
        delete ret.__;
        for (let i = 1; i <= n; i++) ret['' + (i * (i + 1))] = i + 1;
        return ret;
    }

    export function createMap(n: number) {
        const map = new Map<number, number>();
        for (let i = 1; i <= n; i++) map.set(i * (i + 1), i + 1);
        return map;
    }

    export function sumInt(xs: { [key: number]: number }) {
        let s = 0;
        const keys = Object.keys(xs);
        for (let i = 0, _i = keys.length; i < _i; i++) s += xs[+keys[i]];
        return s;
    }

    export function sumStr(xs: { [key: string]: number }) {
        let s = 0;
        const keys = Object.keys(xs);
        for (let i = 0, _i = keys.length; i < _i; i++) s += xs[keys[i]];
        return s;
    }

    export function sumMap(map: Map<number, number>) {
        let s = 0;
        const values = map.keys();
        while (true) {
            const { done, value } = values.next();
            if (done) break;
            s += value;
        }
        return s;
    }

    export function sumCached(xs: WithKeys<number>) {
        let s = 0;
        const keys = xs.keys, data = xs.data;
        for (let i = 0, _i = keys.length; i < _i; i++) s += data[keys[i]];
        return s;
    }

    export function sumKeyMap(xs: MapWithKeys) {
        let s = 0;
        const keys = xs.keys, map = xs.map;
        for (let i = 0, _i = keys.length; i < _i; i++) s += map.get(keys[i])!;
        return s;
    }

    export function run() {
        const N = 1000;
        // const int = createInt(N);
        const map = createMap(N);
        // const str = createStr(N);
        const keys = createCacheKeys(N);
        const keyMap = createMapKeys(N);
        console.log(sumMap(map), sumCached(keys), sumKeyMap(keyMap));
        new B.Suite()
            // .add('c int', () => createInt(N))
            .add('q map', () => sumMap(map))
            .add('c map', () => createMap(N))
            .add('c mk', () => createMapKeys(N))
            // .add('c str', () => createStr(N))
            .add('c cc', () => createCacheKeys(N))
            // .add('q int', () => sumInt(int))
            .add('q mk', () => sumKeyMap(keyMap))
            // .add('q str', () => sumStr(str))
            .add('q cc', () => sumCached(keys))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

IntVsStringIndices.run();

// ObjectVsMap.run();

// testSegments();

// Tuples.run();

// interface AA { kind: 'a' }
// //interface BB { kind: 'b' }
// interface AB { kind: 'a' | 'b' }
// declare const a: AA;
// export const ab: AB = a;
