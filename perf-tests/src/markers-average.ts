import * as B from 'benchmark';

function uint8ForLoop(array: Uint8Array, count: number): number {
    if (count === 0) return 0;
    let sum = 0;
    for (let i = 0; i < count; ++i) {
        sum += array[i] && 1;
    }
    return sum / count;
}

function uint32ForLoop(array: Uint8Array, count: number): number {
    if (count === 0) return 0;

    const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);
    const viewEnd = (count - 4) >> 2;
    const backStart = 4 * viewEnd;

    let sum = 0;
    for (let i = 0; i < viewEnd; ++i) {
        const v = view[i];
        sum += ((v & 0xFF) && 1) + ((v & 0xFF00) && 1) + ((v & 0xFF0000) && 1) + ((v & 0xFF000000) && 1);
    }
    for (let i = backStart; i < count; ++i) {
        sum += array[i] && 1;
    }
    return sum / count;
}

function uint32ForLoopAddOnlyBaseline(array: Uint8Array, count: number): number {
    if (count === 0) return 0;

    const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);
    const viewEnd = (count - 4) >> 2;
    const backStart = 4 * viewEnd;

    let sum = 0;
    for (let i = 0; i < viewEnd; ++i) {
        const v = view[i];
        sum += v;
    }
    for (let i = backStart; i < count; ++i) {
        sum += array[i] && 1;
    }
    return sum / count;
}

function createTypedLut() {
    const lut = new Uint8Array(0x0303 + 1);
    lut[0x0001] = 1;
    lut[0x0002] = 1;
    lut[0x0003] = 1;
    lut[0x0100] = 1;
    lut[0x0200] = 1;
    lut[0x0300] = 1;
    lut[0x0101] = 2;
    lut[0x0201] = 2;
    lut[0x0301] = 2;
    lut[0x0102] = 2;
    lut[0x0202] = 2;
    lut[0x0302] = 2;
    lut[0x0103] = 2;
    lut[0x0203] = 2;
    lut[0x0303] = 2;
    return lut;
}

function createNativeLut() {
    const lut: number[] = [];
    for (let i = 0, il = 0x0303 + 1; i < il; ++i) lut[i] = 0;
    lut[0x0000] = 0;
    lut[0x0001] = 1;
    lut[0x0002] = 1;
    lut[0x0003] = 1;
    lut[0x0100] = 1;
    lut[0x0200] = 1;
    lut[0x0300] = 1;
    lut[0x0101] = 2;
    lut[0x0201] = 2;
    lut[0x0301] = 2;
    lut[0x0102] = 2;
    lut[0x0202] = 2;
    lut[0x0302] = 2;
    lut[0x0103] = 2;
    lut[0x0203] = 2;
    lut[0x0303] = 2;
    return lut;
}

function createTypedLut2bits() {
    const lut = new Uint8Array(256);
    for (let a = 0; a < 4; ++a) {
        for (let b = 0; b < 4; ++b) {
            for (let c = 0; c < 4; ++c) {
                for (let d = 0; d < 4; ++d) {
                    const i = d | c << 2 | b << 4 | a << 6;
                    const v = (a && 1) + (b && 1) + (c && 1) + (d && 1);
                    lut[i] = v;
                    // console.log([a, b, c, d], i, v);
                }
            }
        }
    }
    return lut;
}

const lutNative = createNativeLut();
const lutTyped = createTypedLut();
const lutTyped2bits = createTypedLut2bits();

function uint32ForLoopWithLutNative(array: Uint8Array, count: number): number {
    if (count === 0) return 0;

    const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);
    const viewEnd = (count - 4) >> 2;
    const backStart = 4 * viewEnd;

    let sum = 0;
    for (let i = 0; i < viewEnd; ++i) {
        const v = view[i];
        sum += lutNative[v & 0xFFFF] + lutNative[v >> 16];
    }
    for (let i = backStart; i < count; ++i) {
        sum += array[i] && 1;
    }
    return sum / count;
}

function uint32ForLoopWithLutTyped(array: Uint8Array, count: number): number {
    if (count === 0) return 0;

    const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);
    const viewEnd = (count - 4) >> 2;
    const backStart = 4 * viewEnd;

    let sum = 0;
    for (let i = 0; i < viewEnd; ++i) {
        const v = view[i];
        sum += lutTyped[v & 0xFFFF] + lutTyped[v >> 16];
    }
    for (let i = backStart; i < count; ++i) {
        sum += array[i] && 1;
    }
    return sum / count;
}

function uint32ForLoopWithLut2bits(array: Uint8Array, count: number): number {
    if (count === 0) return 0;

    const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);
    const viewEnd = (count - 4) >> 2;
    const backStart = 4 * viewEnd;

    let sum = 0;
    for (let i = 0; i < viewEnd; ++i) {
        const v = view[i];
        sum += lutTyped2bits[(v >> 18 | v >> 12 | v >> 6 | v) & 0xFF];
    }
    for (let i = backStart; i < count; ++i) {
        sum += array[i] && 1;
    }
    return sum / count;
}

function createData(elements: number, instances: number) {
    const data = new Uint8Array(4000 * 5000);
    const start = Math.floor(instances / 2);
    data.fill(1, start, start + elements);
    return data;
};

export function run(elements: number, instances: number) {
    const suite = new B.Suite();
    const count = elements * instances;
    const data = createData(elements, instances);

    console.log('uint8ForLoop', uint8ForLoop(data, count));
    console.log('uint32ForLoop', uint32ForLoop(data, count));
    console.log('uint32ForLoopWithLutNative', uint32ForLoopWithLutNative(data, count));
    console.log('uint32ForLoopWithLutTyped', uint32ForLoopWithLutTyped(data, count));
    console.log('uint32ForLoopWithLut2bits', uint32ForLoopWithLut2bits(data, count));

    suite
        .add('uint8ForLoop', () => uint8ForLoop(data, count))
        .add('uint32ForLoop', () => uint32ForLoop(data, count))
        .add('uint32ForLoopAddOnlyBaseline', () => uint32ForLoopAddOnlyBaseline(data, count))
        .add('uint32ForLoopWithLutNative', () => uint32ForLoopWithLutNative(data, count))
        .add('uint32ForLoopWithLutTyped', () => uint32ForLoopWithLutTyped(data, count))
        .add('uint32ForLoopWithLut2bits', () => uint32ForLoopWithLut2bits(data, count))
        .on('cycle', (e: any) => console.log(String(e.target)))
        .run();
}

// console.log(createTypedLut2bits());

// run(5000, 4000);
// uint8ForLoop 0.00025
// uint32ForLoop 0.00025
// uint32ForLoopWithLutNative 0.00025
// uint32ForLoopWithLutTyped 0.00025
// uint32ForLoopWithLut2bits 0.00025
// uint8ForLoop x 49.84 ops/sec ±3.30% (64 runs sampled)
// uint32ForLoop x 97.70 ops/sec ±1.71% (71 runs sampled)
// uint32ForLoopAddOnlyBaseline x 220 ops/sec ±2.49% (85 runs sampled)
// uint32ForLoopWithLutNative x 135 ops/sec ±1.71% (76 runs sampled)
// uint32ForLoopWithLutTyped x 137 ops/sec ±1.69% (77 runs sampled)
// uint32ForLoopWithLut2bits x 111 ops/sec ±2.73% (70 runs sampled)
