import * as B from 'benchmark'
import SB from 'mol-base/utils/string-builder'

export namespace Test {
    function createData(n: number) {
        const ret: string[] = [];
        for (let i = 0; i < n; i++) {
            ret[i] = '' + ((100000000 * Math.random() + 1) | 0);
        }
        return ret;
    }

    function build(data: string[], chunkSize: number): SB {
        const sb = SB.create(chunkSize);
        for (let i = 0, _i = data.length; i < _i; i++) {
            SB.writeSafe(sb, data[i]);
        }
        return sb;
    }

    function naive(data: string[]) {
        let ret = '';
        for (let i = 0, _i = data.length; i < _i; i++) ret += data[i];
        return ret;
    }

    function join(data: string[]) {
        let ret = [];
        for (let i = 0, _i = data.length; i < _i; i++) ret[i] = data[i];
        return ret.join('');
    }

    export function run() {
        const data = createData(26 * 100000 * 2);

        const N = 512;
        const suite = new B.Suite();
        suite
            .add(`naive`, () => naive(data)) // cras
            .add(`join`, () => join(data))
            //.add(`${N} chunks`, () => SB.getChunks(build(data, N)))
            .add(`${N} str`, () => SB.getString(build(data, N)))
            .on('cycle', (e: any) => console.log(String(e.target)))
            .run();
    }
}

Test.run();