// array A[] has the items to sort; array B[] is a work array
function BottomUpMergeSort(A: number[]) {
    const n = A.length;
    let src = A, target = new (A as any).constructor(n) as any;
    // Each 1-element run in A is already "sorted".
    // Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted.
    for (let width = 1; width < n; width = 2 * width) {
        // Array A is full of runs of length width.
        for (let i = 0; i < n; i = i + 2 * width) {
            // Merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[]
            // or copy A[i:n-1] to B[] ( if(i+width >= n) )
            BottomUpMerge(src, i, Math.min(i + width, n), Math.min(i + 2 * width, n), target);
        }
        // Now work array B is full of runs of length 2*width.
        // Copy array B to array A for next iteration.
        // A more efficient implementation would swap the roles of A and B.
        const t = src;
        src = target;
        target = t;
        // Now array A is full of runs of length 2*width.
    }
    if (src !== A) {
        for (let i = 0; i < n; i++) src[i] = target[i];
    }
}

//  Left run is A[iLeft :iRight-1].
// Right run is A[iRight:iEnd-1  ].
function BottomUpMerge(A: number[], iLeft: number, iRight: number, iEnd: number, B: number[]) {
    let i = iLeft, j = iRight;
    // While there are elements in the left or right runs...
    for (let k = iLeft; k < iEnd; k++) {
        // If left run head exists and is <= existing right run head.
        const u = A[i], v = A[j];
        if (i < iRight && (j >= iEnd || u <= v)) {
            B[k] = u;
            i = i + 1;
        } else {
            B[k] = v;
            j = j + 1;
        }
    }
}

function mergeSort(a: ArrayLike<number>) {
    BottomUpMergeSort(a as any);
    return a;
}

function heapify(arr: number[], n: number, i: number) {
    const l = 2 * i + 1;  // left = 2*i + 1
    const r = 2 * i + 2;  // right = 2*i + 2

    let largest = i;  // Initialize largest as root

    // If left child is larger than root
    if (l < n && arr[l] > arr[largest])
        largest = l;

    // If right child is larger than largest so far
    if (r < n && arr[r] > arr[largest])
        largest = r;

    // If largest is not root
    if (largest !== i) {
        //swap(arr[i], arr[largest]);
        const t = arr[i];
        arr[i] = arr[largest];
        arr[largest] = t;

        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}

// main function to do heap sort
function heapSort(arr: number[]) {
    const n = arr.length;
    // Build heap (rearrange array)
    for (let i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // One by one extract an element from heap
    for (let i = n - 1; i >= 0; i--) {
        // Move current root to end
        //swap(arr[0], arr[i]);
        const t = arr[0];
        arr[0] = arr[i];
        arr[i] = t;

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
    return arr;
}

function createTestData(n: number) {
    const data = new Int32Array(n); //new Array(n);
    for (let i = 0; i < n; i++) {
        data[i] = (n * Math.random());// | 0;
    }
    return data;

    //return [ 5, 9, 10, 5, 7, 6, 8, 7, 0, 0, 0, 1, 2, 3, 4 ];
}

let swapCount = 0;
function swap(arr: number[], i: number, j: number) {
    const temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

function medianPivot(arr: number[], left: number, right: number) {
    const l = arr[left], r = arr[right], m = arr[(left + right) >> 1];
    if (l > r) return l > m ? Math.max(m, r) : l;
    else return r > m ? Math.max(m, l) : r;
}

const _qsParts = [0, 0];

function partition3(xs: number[], l: number, r: number) {
    const v = medianPivot(xs, l, r)
    let pivot = l, equals = l;
    for (let i = l; i <= r; i++) {
        const t = xs[i];
        if (t < v) {
            /*if (pivot !== i)*/ swap(xs, i, pivot);
            pivot++;
        } else if (t === v) {
            swap(xs, i, pivot);
            swap(xs, pivot, equals);
            equals++;
            pivot++;
        }
    }
    for (let i = l; i < equals; i++) {
        swap(xs, i, l + pivot - i - 1)
    }
    _qsParts[0] = pivot - equals + l;
    _qsParts[1] = pivot - 1;
}

function partition3_1(xs: number[], l: number, r: number) {
    const v = medianPivot(xs, l, r);
    // console.log('P', v);
    // console.log('I', xs.slice(l, r + 1));
    let equals = l, tail = r;

    while (xs[tail] > v) --tail;
    for (let i = l; i <= tail; i++) {
        const t = xs[i];
        if (t > v) {
            swap(xs, i, tail);
            tail--;
            while (xs[tail] > v) --tail;
            i--;
        } else if (t === v) {
            swap(xs, i, equals);
            equals++;
        }
    }
    //console.log('M', xs.slice(l, r + 1), equals, tail);
    for (let i = l; i < equals; i++) {
        swap(xs, i, l + tail - i)
    }
    //console.log('F', xs.slice(l, r + 1), tail - equals + l + 1, tail);
    _qsParts[0] = tail - equals + l + 1;
    _qsParts[1] = tail;
}

function insertionSort(xs: number[], start: number, end: number) {
    for (let i = start + 1; i <= end; i++) {
        const key = xs[i];
        let j = i - 1;
        while (j >= 0 && xs[j] > key) {
            xs[j + 1] = xs[j];
            j = j - 1;
        }
        xs[j + 1] = key;
    }
}

function quickSort(xs: number[], low: number, high: number) {
    while (low < high) {
        if (high - low < 16) {
            insertionSort(xs, low, high);
            return;
        }

        partition3_1(xs, low, high);
        const li = _qsParts[0], ri = _qsParts[1];

        if (li - low < high - ri) {
            quickSort(xs, low, li - 1);
            low = ri + 1;
        } else { // Else recur for right part
            quickSort(xs, ri + 1, high);
            high = li - 1;
        }
    }
}

function checkSorted(arr: ArrayLike<number>) {
    for (let i = 0; i < arr.length - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            console.log('not sorted');
            return;
        }
    }
    console.log('sorted');
}

(function test() {
    // console.log(medianPivot([1, 2, 3], 0, 2))
    // console.log(medianPivot([1, 3, 2], 0, 2))
    // console.log(medianPivot([2, 1, 3], 0, 2))
    // console.log(medianPivot([3, 1, 2], 0, 2))
    // console.log(medianPivot([2, 3, 1], 0, 2))
    // console.log(medianPivot([3, 2, 1], 0, 2))

    const n = 10000;

    Array.prototype.sort.call(createTestData(n), (a: number, b: number) => a - b);
    mergeSort(createTestData(n));

    let sd;



    sd = createTestData(n);
    quickSort(sd as any, 0, sd.length - 1);

    //console.log(sd);

    // sd = createTestData(n);
    // console.time('merge');
    // mergeSort(sd);
    // console.timeEnd('merge');
    // checkSorted(sd);

    sd = createTestData(n);

    checkSorted(sd);
    console.time('heap');
    heapSort(sd as any);
    console.timeEnd('heap');
    checkSorted(sd);

    console.time('heap-sorted');
    heapSort(sd as any);
    console.timeEnd('heap-sorted');
    checkSorted(sd);

    console.log('--------------');

    //console.log('quick', sd);

    sd = createTestData(n);
    checkSorted(sd);
    console.time('qs');
    quickSort(sd as any, 0, sd.length - 1);
    console.timeEnd('qs');
    checkSorted(sd);

    console.time('qs-sorted');
    quickSort(sd as any, 0, sd.length - 1);
    console.timeEnd('qs-sorted');
    checkSorted(sd);

    const reverseSorted = new Int32Array(n);
    for (let i = 0; i < n; i++) {
        reverseSorted[i] = sd[n - i - 1];
    }

    console.time('qs-reverse-sorted');
    quickSort(reverseSorted as any, 0, reverseSorted.length - 1);
    console.timeEnd('qs-reverse-sorted');
    checkSorted(reverseSorted);

    console.log('swap count', swapCount);

    console.log('--------------');

    sd = createTestData(n);
    checkSorted(sd);
    console.time('native');
    sd.sort((a: number, b: number) => a - b);
    console.timeEnd('native');
    checkSorted(sd);

    console.time('native-sorted');
    sd.sort((a: number, b: number) => a - b);
    console.timeEnd('native-sorted');
    checkSorted(sd);
    //console.log(sd);

}())