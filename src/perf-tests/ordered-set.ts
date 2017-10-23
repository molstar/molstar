// import * as B from 'benchmark'
// import OrderedSet from '../structure/collections/ordered-set'
// import OrderedSet1 from '../structure/collections/ordered-set.1'

// const range = OrderedSet.ofRange(0, 100);
// const range1 = OrderedSet1.ofRange(0, 100);
// // const pairSet = IntPair.pack1(0, 100);

// // namespace PairSet {
// //     const pair = IntPair.zero();
// //     export function has(p: number, x: number) {
// //         IntPair.unpack(p, pair);
// //         return x >= pair.fst && x <= pair.snd;
// //     }
// // }

// const suite = new B.Suite();

// const values: number[] = [];
// for (let i = 0; i < 1000000; i++) values[i] = (Math.random() * 1000) | 0;

// let idx = 0;

// suite
//     .add('range', () => range.has(idx % values.length))
//     .add('range1', () => OrderedSet1.has(range1, idx % values.length))
//     .on('cycle', (e: any) => {
//         console.log(String(e.target));
//     })
//     .run();
