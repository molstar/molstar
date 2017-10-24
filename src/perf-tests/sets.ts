import * as B from 'benchmark'
import IntTuple from '../common/collections/int-tuple'
import OrdSet from '../common/collections/ordered-set'
import AtomSet from '../mol-data/atom-set'

namespace Iteration {
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
        const it = AtomSet.values(ms);
        for (let v = it.move(); !it.done; v = it.move()) s += v.snd;
        return s;
    }

    export function elementAt() {
        let s = 0;
        for (let i = 0, _i = AtomSet.size(ms); i < _i; i++) s += IntTuple.snd(AtomSet.getAt(ms, i));
        return s;
    }

    export function manual() {
        let s = 0;
        const keys = AtomSet.keys(ms);
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
        for (let i = 0, _i = AtomSet.keyCount(ms); i < _i; i++) {
            const set = AtomSet.getByIndex(ms, i);
            for (let j = 0, _j = OrdSet.size(set); j < _j; j++) {
                s += OrdSet.getAt(set, j);
            }
        }
        return s;
    }
}

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