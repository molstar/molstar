/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Matrix from './matrix';

export namespace EVD {
    export interface Cache {
        size: number,
        matrix: Matrix,
        eigenValues: number[],
        D: number[],
        E: number[]
    }

    export function createCache(size: number): Cache {
        return {
            size,
            matrix: Matrix.create(size, size),
            eigenValues: <any>new Float64Array(size),
            D: <any>new Float64Array(size),
            E: <any>new Float64Array(size)
        };
    }

    /**
     * Computes EVD and stores the result in the cache.
     */
    export function compute(cache: Cache): void {
        symmetricEigenDecomp(cache.size, cache.matrix.data as number[], cache.eigenValues, cache.D, cache.E);
    }
}

// The EVD code has been adapted from Math.NET, MIT license, Copyright (c) 2002-2015 Math.NET
//
// THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function symmetricEigenDecomp(order: number, matrixEv: number[], vectorEv: number[], d: number[], e: number[]) {
    for (let i = 0; i < order; i++) {
        e[i] = 0.0;
    }

    let om1 = order - 1;
    for (let i = 0; i < order; i++) {
        d[i] = matrixEv[i * order + om1];
    }

    symmetricTridiagonalize(matrixEv, d, e, order);
    symmetricDiagonalize(matrixEv, d, e, order);

    for (let i = 0; i < order; i++) {
        vectorEv[i] = d[i];
    }
}

function symmetricTridiagonalize(a: number[], d: number[], e: number[], order: number) {
    // Householder reduction to tridiagonal form.
    for (let i = order - 1; i > 0; i--) {
        // Scale to avoid under/overflow.
        let scale = 0.0;
        let h = 0.0;

        for (let k = 0; k < i; k++) {
            scale = scale + Math.abs(d[k]);
        }

        if (scale === 0.0) {
            e[i] = d[i - 1];
            for (let j = 0; j < i; j++) {
                d[j] = a[(j * order) + i - 1];
                a[(j * order) + i] = 0.0;
                a[(i * order) + j] = 0.0;
            }
        } else {
            // Generate Householder vector.
            for (let k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }

            let f = d[i - 1];
            let g = Math.sqrt(h);
            if (f > 0) {
                g = -g;
            }

            e[i] = scale * g;
            h = h - (f * g);
            d[i - 1] = f - g;

            for (let j = 0; j < i; j++) {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.
            for (let j = 0; j < i; j++) {
                f = d[j];
                a[(i * order) + j] = f;
                g = e[j] + (a[(j * order) + j] * f);

                for (let k = j + 1; k <= i - 1; k++) {
                    g += a[(j * order) + k] * d[k];
                    e[k] += a[(j * order) + k] * f;
                }

                e[j] = g;
            }

            f = 0.0;

            for (let j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }

            let hh = f / (h + h);

            for (let j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }

            for (let j = 0; j < i; j++) {
                f = d[j];
                g = e[j];

                for (let k = j; k <= i - 1; k++) {
                    a[(j * order) + k] -= (f * e[k]) + (g * d[k]);
                }

                d[j] = a[(j * order) + i - 1];
                a[(j * order) + i] = 0.0;
            }
        }

        d[i] = h;
    }

    // Accumulate transformations.
    for (let i = 0; i < order - 1; i++) {
        a[(i * order) + order - 1] = a[(i * order) + i];
        a[(i * order) + i] = 1.0;
        let h = d[i + 1];
        if (h !== 0.0) {
            for (let k = 0; k <= i; k++) {
                d[k] = a[((i + 1) * order) + k] / h;
            }

            for (let j = 0; j <= i; j++) {
                let g = 0.0;
                for (let k = 0; k <= i; k++) {
                    g += a[((i + 1) * order) + k] * a[(j * order) + k];
                }

                for (let k = 0; k <= i; k++) {
                    a[(j * order) + k] -= g * d[k];
                }
            }
        }

        for (let k = 0; k <= i; k++) {
            a[((i + 1) * order) + k] = 0.0;
        }
    }

    for (let j = 0; j < order; j++) {
        d[j] = a[(j * order) + order - 1];
        a[(j * order) + order - 1] = 0.0;
    }

    a[(order * order) - 1] = 1.0;
    e[0] = 0.0;
}

function symmetricDiagonalize(a: number[], d: number[], e: number[], order: number) {
    const maxiter = 1000;

    for (let i = 1; i < order; i++) {
        e[i - 1] = e[i];
    }

    e[order - 1] = 0.0;

    let f = 0.0;
    let tst1 = 0.0;
    let eps = Math.pow(2, -53); // DoubleWidth = 53
    for (let l = 0; l < order; l++) {
        // Find small subdiagonal element
        tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
        let m = l;
        while (m < order) {
            if (Math.abs(e[m]) <= eps * tst1) {
                break;
            }

            m++;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.
        if (m > l) {
            let iter = 0;
            do {
                iter = iter + 1; // (Could check iteration count here.)

                // Compute implicit shift
                let g = d[l];
                let p = (d[l + 1] - g) / (2.0 * e[l]);
                let r = hypotenuse(p, 1.0);
                if (p < 0) {
                    r = -r;
                }

                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);

                let dl1 = d[l + 1];
                let h = g - d[l];
                for (let i = l + 2; i < order; i++) {
                    d[i] -= h;
                }

                f = f + h;

                // Implicit QL transformation.
                p = d[m];
                let c = 1.0;
                let c2 = c;
                let c3 = c;
                let el1 = e[l + 1];
                let s = 0.0;
                let s2 = 0.0;
                for (let i = m - 1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypotenuse(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = (c * d[i]) - (s * g);
                    d[i + 1] = h + (s * ((c * g) + (s * d[i])));

                    // Accumulate transformation.
                    for (let k = 0; k < order; k++) {
                        h = a[((i + 1) * order) + k];
                        a[((i + 1) * order) + k] = (s * a[(i * order) + k]) + (c * h);
                        a[(i * order) + k] = (c * a[(i * order) + k]) - (s * h);
                    }
                }

                p = (-s) * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence. If too many iterations have been performed,
                // throw exception that Convergence Failed
                if (iter >= maxiter) {
                    throw 'SVD: Not converging.';
                }
            } while (Math.abs(e[l]) > eps * tst1);
        }

        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.
    for (let i = 0; i < order - 1; i++) {
        let k = i;
        let p = d[i];
        for (let j = i + 1; j < order; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }

        if (k !== i) {
            d[k] = d[i];
            d[i] = p;
            for (let j = 0; j < order; j++) {
                p = a[(i * order) + j];
                a[(i * order) + j] = a[(k * order) + j];
                a[(k * order) + j] = p;
            }
        }
    }
}

function hypotenuse(a: number, b: number) {
    if (Math.abs(a) > Math.abs(b)) {
        let r = b / a;
        return Math.abs(a) * Math.sqrt(1 + (r * r));
    }

    if (b !== 0.0) {
        let r = a / b;
        return Math.abs(b) * Math.sqrt(1 + (r * r));
    }

    return 0.0;
}