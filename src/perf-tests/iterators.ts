import * as B from 'benchmark';
import It from '../mol-data/iterator';

function createData(n: number) {
    const data = []; // new Int32Array(n);
    let last = (15 * Math.random()) | 0;
    for (let i = 0; i < n; i++) {
        data[i] = last;
        last += (15 * Math.random()) | 0;
    }
    return data;
}

export namespace Iterators {
    const data = createData(100000);

    export function forLoop() {
        let sum = 0;
        for (let i = 0, _i = data.length; i < _i; i++) {
            sum += data[i];
        }
        return sum;
    }

    export function forOf() {
        let sum = 0;
        for (const e of data) {
            sum += e;
        }
        return sum;
    }

    export function forEach() {
        const ctx = { sum: 0 };
        data.forEach(function (this: typeof ctx, v: number) { this.sum += v; }, ctx);
        return ctx.sum;
    }

    export function forEachAllParams() {
        const ctx = { sum: 0 };
        data.forEach(function (this: typeof ctx, v: number, _: any, __: any) { this.sum += v; }, ctx);
        return ctx.sum;
    }

    export function forEachClosure() {
        let sum = 0;
        data.forEach(v => sum += v);
        return sum;
    }

    export function forEachClosureAll() {
        let sum = 0;
        data.forEach((v, _, __) => sum += v);
        return sum;
    }


    export function forEachClosureAllFunction() {
        let sum = 0;
        data.forEach(function (v, _, __) { sum += v; });
        return sum;
    }

    interface ES6Iterator {
        [Symbol.iterator](): ES6Iterator,
        done: boolean;
        value: number;
        next(): { done: boolean, value: number }
        reset(data: any[]): ES6Iterator;
    }

    class _MutableES6Iterator implements ES6Iterator {
        done = true;
        value = 0;

        private xs: any[] = void 0 as any;
        private index: number = -1;
        private length: number = 0;

        [Symbol.iterator]() { return this; };

        next() {
            const index = ++this.index;
            if (index < this.length) this.value = this.xs[index];
            else this.done = true;
            return this;
        }

        reset(xs: any[]) {
            this.value = xs[0];
            this.length = xs.length;
            this.done = false;
            this.xs = xs;
            this.index = -1;
            return this;
        }
    }

    class _ImmutableES6Iterator implements ES6Iterator {
        done = true;
        value = 0;

        private xs: any[] = void 0 as any;
        private index: number = -1;
        private length: number = 0;

        [Symbol.iterator]() { return this; };

        next() {
            const index = ++this.index;
            if (index < this.length) this.value = this.xs[index];
            else this.done = true;
            return { done: this.done, value: this.value };
        }

        reset(xs: any[]) {
            this.value = xs[0];
            this.length = xs.length;
            this.done = false;
            this.xs = xs;
            this.index = -1;
            return this;
        }
    }

    export function mutableES6Iterator() {
        const it = new _MutableES6Iterator();
        let sum = 0;
        for (let e = it.reset(data).next().value; !it.done; e = it.next().value) {
            sum += e;
        }
        return sum;
    }

    // export function mutableES6IteratorOf() {
    //     const it = new _ImmutableES6Iterator();
    //     let sum = 0;
    //     for (const e of it.reset(data)) {
    //         sum += e;
    //     }
    //     return sum;
    // }

    export function immutableES6Iterator() {
        const it = new _ImmutableES6Iterator();
        let sum = 0;
        it.reset(data);
        while (true) {
            const { value, done } = it.next();
            if (done) break;
            sum += value;
        }
        return sum;
    }

    // export function immutableES6IteratorOf() {
    //     const it = new _MutableES6Iterator();
    //     let sum = 0;
    //     for (const e of it.reset(data)) {
    //         sum += e;
    //     }
    //     return sum;
    // }

    interface MutableIterator {
        done: boolean;
        next(): number;
        start(data: any[]): number;
    }

    class _MutableIterator implements MutableIterator {
        done = true;

        private xs: any[] = void 0 as any;
        private index: number = -1;
        private length: number = 0;

        next() {
            const index = ++this.index;
            if (index < this.length) return this.xs[index];
            else {
                this.done = true;
                return 0;
            }
        }

        start(xs: any[]) {
            this.length = xs.length;
            this.done = !this.length;
            this.xs = xs;
            this.index = 0;
            return this.done ? 0 : this.xs[0];
        }
    }

    export function mutableIterator() {
        const it = new _MutableIterator();
        let sum = 0;
        for (let e = it.start(data); !it.done; e = it.next()) {
            sum += e;
        }
        return sum;
    }

    export function run() {
        const suite = new B.Suite();

        suite
            .add('for', () => Iterators.forLoop())
            .add('forOf', () => Iterators.forOf())
            .add('forEach', () => Iterators.forEach())
            .add('forEach all params', () => Iterators.forEachAllParams())
            .add('forEachClosure', () => Iterators.forEachClosure())
            .add('forEachClosure all', () => Iterators.forEachClosureAll())
            .add('forEachClosure all function', () => Iterators.forEachClosureAllFunction())
            .add('mutableIterator ES6', () => Iterators.mutableES6Iterator())
            // .add('mutableIteratorOf ES6', () => Iterators.mutableES6IteratorOf())
            .add('immutableIterator ES6', () => Iterators.immutableES6Iterator())
            // .add('immutableIteratorOf ES6', () => Iterators.immutableES6IteratorOf())
            .add('mutableIterator', () => Iterators.mutableIterator())
            .on('cycle', (e: any) => {
                console.log(String(e.target));
            })
            // .on('complete', function (this: any) {
            //     console.log('Fastest is ' + this.filter('fastest').map('name'));
            // })
            .run();
    }
}

const it = It.Array([1, 2, 3]);
while (it.hasNext) {
    console.log(it.move());
}