import * as B from 'benchmark'
import IntTuple from '../mol-base/collections/int-tuple'
import OrdSet from '../mol-base/collections/ordered-set'
import AtomSet from '../mol-data/atom-set'

export namespace Iteration {
    const U = 1000, V = 2500;

    const control: number[] = [];
    const sets = Object.create(null);
    for (let i = 1; i < U; i++) {
        const set = [];
        for (let j = 1; j < V; j++) {
            control[control.length] = j * j + 1;
            set[set.length] = j * j + 1;
        }
        sets[i * i] = OrdSet.ofSortedArray(set);
    }
    const ms = AtomSet.create(sets);

    export function native() {
        let s = 0;
        for (let i = 0, _i = control.length; i < _i; i++) s += control[i];
        return s;
    }

    export function iterators() {
        let s = 0;
        const it = AtomSet.atoms(ms);
        for (let v = it.move(); !it.done; v = it.move()) s += IntTuple.snd(v);
        return s;
    }

    export function elementAt() {
        let s = 0;
        for (let i = 0, _i = AtomSet.atomCount(ms); i < _i; i++) s += IntTuple.snd(AtomSet.getAtomAt(ms, i));
        return s;
    }

    export function manual() {
        let s = 0;
        const keys = AtomSet.units(ms);
        for (let i = 0, _i = OrdSet.size(keys); i < _i; i++) {
            const set = AtomSet.getByKey(ms, OrdSet.getAt(keys, i));
            for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
                s += OrdSet.getAt(set, j);
            }
        }
        return s;
    }

    export function manual1() {
        let s = 0;
        for (let i = 0, _i = AtomSet.unitCount(ms); i < _i; i++) {
            const set = AtomSet.getByIndex(ms, i);
            for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
                s += OrdSet.getAt(set, j);
            }
        }
        return s;
    }

    export function run() {
        const suite = new B.Suite();

        // const values: number[] = [];
        // for (let i = 0; i < 1000000; i++) values[i] = (Math.random() * 1000) | 0;

        console.log(Iteration.native(), Iteration.iterators(), Iteration.elementAt(), Iteration.manual(), Iteration.manual1());

        suite
            .add('native', () => Iteration.native())
            .add('iterators', () => Iteration.iterators())
            .add('manual', () => Iteration.manual())
            .add('manual1', () => Iteration.manual1())
            .add('el at', () => Iteration.elementAt())
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

export namespace Union {
    function createArray(n: number) {
        const data = new Int32Array(n);
        let c = (Math.random() * 100) | 0;
        for (let i = 0; i < n; i++) {
            data[i] = c;
            c += 1 + (Math.random() * 100) | 0
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

        return { ordSet: OrdSet.ofSortedArray(array), obj, set }
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

    function _setAdd(this: Set<number>, x: number) { this.add(x) }
    function unionS(a: Set<number>, b: Set<number>) {
        const ret = new Set<number>();
        a.forEach(_setAdd, ret);
        b.forEach(_setAdd, ret);
        return ret;
    }

    function _setInt(this: { set: Set<number>, other: Set<number> }, x: number) { if (this.other.has(x)) this.set.add(x) }
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

export namespace Build {
    function createSorted() {
        const b = AtomSet.SortedBuilder(AtomSet.Empty);
        for (let i = 0; i < 100; i++) {
            for (let j = 0; j < 100; j++) {
                b.add(i, j);
            }
        }
        return b.getSet();
    }

    function createByUnit() {
        const b = AtomSet.SortedBuilder(AtomSet.Empty);
        for (let i = 0; i < 100; i++) {
            b.beginUnit();
            for (let j = 0; j < 100; j++) {
                b.addToUnit(j);
            }
            b.commitUnit(i);
        }
        return b.getSet();
    }


    export function run() {
        const suite = new B.Suite();
        suite
            .add('create sorted', () => createSorted())
            .add('create by unit', () => createByUnit())
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

export namespace Tuples {
    function createData(n: number) {
        const ret: IntTuple[] = new Float64Array(n) as any;
        for (let i = 0; i < n; i++) {
            ret[i] = IntTuple.create(i, i * i + 1);
        }
        return ret;
    }

    function sum1(data: ArrayLike<IntTuple>) {
        let s = 0;
        for (let i = 0, _i = data.length; i < _i; i++) {
            s += IntTuple.fst(data[i]) + IntTuple.snd(data[i]);
        }
        return s;
    }

    function sum2(data: ArrayLike<IntTuple>) {
        let s = 0;
        for (let i = 0, _i = data.length; i < _i; i++) {
            const t = data[i];
            s += IntTuple.fst(t) + IntTuple.snd(t);
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
    const segs = OrdSet.ofSortedArray([0, 4, 10, 12, 13, 15, 25]);
    const it = OrdSet.segments(segs, data);

    for (let s = it.move(); !it.done; s = it.move()) {
        for (let j = s.start; j < s.end; j++) {
            console.log(`${s.segment}: ${OrdSet.getAt(data, j)}`);
        }
    }
}

testSegments();

//Tuples.run();

// interface AA { kind: 'a' }
// //interface BB { kind: 'b' }
// interface AB { kind: 'a' | 'b' }
// declare const a: AA;
// export const ab: AB = a;
