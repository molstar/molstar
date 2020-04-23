/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Tensor as T } from '../tensor';
import { Mat4 } from '../3d';

describe('tensor', () => {
    it('vector', () => {
        const V = T.Vector(3);
        const data = V.create();

        V.set(data, 0, 1);
        V.set(data, 1, 2);
        V.set(data, 2, 3);

        expect(data).toEqual(new Float64Array([1, 2, 3]));
        expect(V.get(data, 0)).toEqual(1);
        expect(V.get(data, 1)).toEqual(2);
        expect(V.get(data, 2)).toEqual(3);
    });

    it('matrix cm', () => {
        const M = T.ColumnMajorMatrix(3, 2, Int32Array);
        const data = M.create();

        // rows: [ [0, 1], [1, 2], [2, 3]  ]
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 2; j++) {
                M.set(data, i, j, i + j);
            }
        }

        expect(data).toEqual(new Int32Array([0, 1, 2, 1, 2, 3]));
    });

    it('matrix rm', () => {
        const M = T.RowMajorMatrix(3, 2, Int32Array);
        const data = M.create();

        // rows: [ [0, 1], [1, 2], [2, 3]  ]
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 2; j++) {
                M.set(data, i, j, i + j);
            }
        }
        expect(data).toEqual(new Int32Array([0, 1, 1, 2, 2, 3]));
    });

    it('mat4 equiv', () => {
        const M = T.ColumnMajorMatrix(4, 4);
        const data = M.create();
        const m = Mat4.zero();

        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                const v = (i + 1) * (j + 2);
                M.set(data, i, j, v);
                Mat4.setValue(m, i, j, v);
            }
        }

        for (let i = 0; i < 16; i++) expect(data[i]).toEqual(m[i]);
    });

    it('2d ij', () => {
        const M = T.Space([3, 4], [0, 1]);
        const data = M.create();
        const exp = new Float64Array(3 * 4);

        let o = 0;
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 4; j++) {
                M.set(data, i, j, o);
                expect(M.get(data, i, j)).toBe(o);
                exp[o] = o;
                o++;
            }
        }

        expect(data).toEqual(exp);
    });

    it('2d ji', () => {
        const M = T.Space([3, 4], [1, 0]);
        const data = M.create();
        const exp = new Float64Array(3 * 4);

        let o = 0;
        for (let j = 0; j < 4; j++) {
            for (let i = 0; i < 3; i++) {
                M.set(data, i, j, o);
                expect(M.get(data, i, j)).toBe(o);
                exp[o] = o;
                o++;
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d ijk', () => {
        const M = T.Space([3, 4, 5], [0, 1, 2]);
        const data = M.create();
        const exp = new Float64Array(3 * 4 * 5);

        let o = 0;
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 4; j++) {
                for (let k = 0; k < 5; k++) {
                    M.set(data, i, j, k, o);
                    expect(M.get(data, i, j, k)).toBe(o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d ikj', () => {
        const M = T.Space([3, 3, 3], [0, 2, 1]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3);

        let o = 0;
        for (let i = 0; i < 3; i++) {
            for (let k = 0; k < 3; k++) {
                for (let j = 0; j < 3; j++) {
                    M.set(data, i, j, k, o);
                    expect(M.get(data, i, j, k)).toBe(o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d jik', () => {
        const M = T.Space([3, 3, 3], [1, 0, 2]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3);

        let o = 0;
        for (let j = 0; j < 3; j++) {
            for (let i = 0; i < 3; i++) {
                for (let k = 0; k < 3; k++) {
                    M.set(data, i, j, k, o);
                    expect(M.get(data, i, j, k)).toBe(o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });
    it('3d jki', () => {
        const M = T.Space([3, 3, 3], [1, 2, 0]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3);

        let o = 0;
        for (let j = 0; j < 3; j++) {
            for (let k = 0; k < 3; k++) {
                for (let i = 0; i < 3; i++) {
                    M.set(data, i, j, k, o);
                    expect(M.get(data, i, j, k)).toBe(o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d kij', () => {
        const M = T.Space([3, 3, 3], [2, 0, 1]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3);

        let o = 0;
        for (let k = 0; k < 3; k++) {
            for (let i = 0; i < 3; i++) {
                for (let j = 0; j < 3; j++) {
                    M.set(data, i, j, k, o);
                    expect(M.get(data, i, j, k)).toBe(o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('3d kji', () => {
        const M = T.Space([3, 3, 3], [2, 1, 0]);
        const data = M.create();
        const exp = new Float64Array(3 * 3 * 3);

        let o = 0;
        for (let k = 0; k < 3; k++) {
            for (let j = 0; j < 3; j++) {
                for (let i = 0; i < 3; i++) {
                    M.set(data, i, j, k, o);
                    expect(M.get(data, i, j, k)).toBe(o);
                    exp[o] = o;
                    o++;
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('4d jikl', () => {
        const M = T.Space([2, 3, 4, 5], [1, 0, 2, 3]);
        const data = M.create();
        const exp = new Float64Array(2 * 3 * 4 * 5);

        let o = 0;
        for (let j = 0; j < 3; j++) {
            for (let i = 0; i < 2; i++) {
                for (let k = 0; k < 4; k++) {
                    for (let l = 0; l < 5; l++) {
                        M.set(data, i, j, k, l, o);
                        expect(M.get(data, i, j, k, l)).toBe(o);
                        exp[o] = o;
                        o++;
                    }
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('4d jilk', () => {
        const M = T.Space([2, 3, 4, 5], [1, 0, 3, 2]);
        const data = M.create();
        const exp = new Float64Array(2 * 3 * 4 * 5);

        let o = 0;
        for (let j = 0; j < 3; j++) {
            for (let i = 0; i < 2; i++) {
                for (let l = 0; l < 5; l++) {
                    for (let k = 0; k < 4; k++) {
                        M.set(data, i, j, k, l, o);
                        expect(M.get(data, i, j, k, l)).toBe(o);
                        exp[o] = o;
                        o++;
                    }
                }
            }
        }

        expect(data).toEqual(exp);
    });

    it('indexing', () => {
        function permutations<T>(inputArr: T[]): T[][] {
            let result: T[][] = [];
            function permute(arr: any, m: any = []) {
                if (arr.length === 0) {
                    result.push(m);
                } else {
                    for (let i = 0; i < arr.length; i++) {
                        let curr = arr.slice();
                        let next = curr.splice(i, 1);
                        permute(curr.slice(), m.concat(next));
                    }
                }
            }
            permute(inputArr);

            return result;
        }

        for (let dim = 1; dim <= 5; dim++) {
            const axes = [], dims: number[] = [];
            const u: number[] = [], v: number[] = [];

            for (let i = 0; i < dim; i++) {
                axes.push(i);
                dims.push(3);
                u.push(0);
                v.push(0);
            }

            const forEachDim = (space: T.Space, d: number): boolean => {
                if (d === dim) {
                    const o = space.dataOffset(...u);
                    space.getCoords(o, v);

                    for (let e = 0; e < dims.length; e++) {
                        expect(u[e]).toEqual(v[e]);
                        return false;
                    }
                } else {
                    for (let i = 0; i < dims[d]; i++) {
                        u[d] = i;
                        if (!forEachDim(space, d + 1)) return false;
                    }
                }
                return true;
            };

            for (const ao of permutations(axes)) {
                const space = T.Space(dims, ao);
                if (!forEachDim(space, 0)) break;
            }
        }
    });
});