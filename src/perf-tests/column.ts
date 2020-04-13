import * as B from 'benchmark';
import { Column as C } from '../mol-data/db';

export namespace Column {
    function createData(n: number) {
        const ret = new Float32Array(n);
        for (let i = 0; i < n; i++) {
            ret[i] = i * i + 1;
        }
        return ret;
    }

    function raw(xs: ArrayLike<number>) {
        let sum = 0;
        for (let i = 0, _i = xs.length; i < _i; i++) {
            sum += xs[i];
        }
        return sum;
    }

    function column(col: C<number>) {
        let sum = 0;
        for (let i = 0, _i = col.rowCount; i < _i; i++) {
            sum += col.value(i);
        }
        return sum;
    }

    function column1(col: C<number>) {
        let sum = 0;
        for (let i = 0, _i = col.rowCount; i < _i; i++) {
            sum += col.value(i);
        }
        return sum;
    }

    function val(i: number) { return i * i + 1; }

    export function runMono() {
        const suite = new B.Suite();
        const data = createData(1000);
        const nativeData = [...data as any];
        const col = C.ofArray({ array: data, schema: C.Schema.float });
        const lambda = C.ofLambda({ value: val, rowCount: data.length, schema: C.Schema.float });
        const cnst = C.ofConst(10, data.length, C.Schema.float);
        suite
            .add('raw', () => raw(data))
            .add('native raw', () => raw(nativeData))
            .add('arraycol', () => column(col))
            .add('arraycol1', () => column(col))
            .add('const', () => column1(cnst))
            .add('arraycol2', () => column(col))
            .add('lambda', () => column1(lambda))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }

    export function runPoly() {
        const suite = new B.Suite();
        const data = createData(10000);
        const nativeData = [...data as any];
        const col = C.ofArray({ array: data, schema: C.Schema.float });
        const lambda = C.ofLambda({ value: val, rowCount: data.length, schema: C.Schema.float });
        const cnst = C.ofConst(10, data.length, C.Schema.float);
        suite
            .add('raw', () => raw(data))
            .add('native raw', () => raw(nativeData))
            .add('arraycol', () => column(col))
            .add('const', () => column(cnst))
            .add('arraycol2', () => column(col))
            .add('lambda', () => column(lambda))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

Column.runMono();
Column.runPoly();