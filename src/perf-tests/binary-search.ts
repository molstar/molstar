
function createData(n: number) {
    const data = [];//new Int32Array(n);
    let last = (15 * Math.random()) | 0;
    for (let i = 0; i < n; i++) {
        data[i] = last;
        last += (15 * Math.random()) | 0;
    }
    return data;
}

function binarySearchHelper(list: ArrayLike<number>, t: number) {
    let min = 0, max = list.length - 1;
    while (min <= max) {
        if (min + 11 > max) {
            for (let i = min; i <= max; i++) {
                if (t === list[i]) return i;
            }
            return -1;
        }

        const mid = (min + max) >> 1;
        const v = list[mid];
        if (t < v) max = mid - 1;
        else if (t > v) min = mid + 1;
        else return mid;
    }
    return -1;
}

function objSearch(obj: any, val: number) {
    return typeof obj[val] !== 'undefined';
}

function setSearch(set: Set<number>, val: number) {
    return set.has(val);
}

function prepare(list: ArrayLike<number>) {
    const obj = Object.create(null), set = new Set<number>();
    for (let i = 0; i < list.length; i++) {
        const v = list[i];
        obj[v] = i;
        set.add(v);
    }

    return { list, obj, set };
}

function prepareSet(list: ArrayLike<number>) {
    const set = new Set<number>();
    for (let i = 0; i < list.length; i++) {
        const v = list[i];
        set.add(v);
    }
    return set;
}

function prepareObj(list: ArrayLike<number>) {
    const obj = Object.create(null);
    for (let i = 0; i < list.length; i++) {
        const v = list[i];
        obj[v] = i;
    }

    return obj;
}

function testBinary(list: ArrayLike<number>, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (binarySearchHelper(list, points[i]) >= 0) r += points[i];
    }
    return r;
}

function testObj(obj: any, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (objSearch(obj, points[i])) r += points[i];
    }
    return r;
}

function testSet(set: Set<number>, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (setSearch(set, points[i])) r += points[i];
    }
    return r;
}

function run(f: () => number, n: number) {
    for (let i = 0; i < n; i++) f();
}

(function () {
    const size = 100000;
    const list = createData(size);
    const queryPoints = createData(size);

    let obj = prepareObj(list);
    let set = prepareSet(list);

    console.log('list', testBinary(list, queryPoints));
    console.log('obj', testObj(obj, queryPoints));
    console.log('set', testSet(set, queryPoints));

    console.time('obj');
    run(() => testObj(obj, queryPoints), 100);
    console.timeEnd('obj');

    console.time('set');
    run(() => testSet(set, queryPoints), 100);
    console.timeEnd('set');

    console.time('bin-search');
    run(() => testBinary(list, queryPoints), 100);
    console.timeEnd('bin-search');

    console.time('prepare-obj');
    run(() => prepareObj(list), 1);
    console.timeEnd('prepare-obj');

    console.time('prepare-set');
    run(() => prepareSet(list).size, 1);
    console.timeEnd('prepare-set');
}())
