/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 * 
 * QCProt: Quaternion-based characteristic polynomial
 * 
 * This code is based on
 * https://web.archive.org/web/20240203104803/https://theobald.brandeis.edu/qcp/main.c
 * Douglas L. Theobald (2005) "Rapid calculation of RMSD using a
 * quaternion-based characteristic polynomial." Acta Crystallographica A
 * 61(4):478-480.
 *
 * Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009) "Fast
 * determination of the optimal rotational matrix for macromolecular
 * superpositions." Journal of Computational Chemistry 31(7):1561-1563.
 *
 * Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer. Redistributions
 * in binary form must reproduce the above copyright notice, this list of
 * conditions and the following disclaimer in the documentation and/or other
 * materials provided with the distribution. Neither the name of the
 * <ORGANIZATION> nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior
 * written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
export class QCProt {
    // The quaternion components.
    q1: number;
    q2: number;
    q3: number;
    q4: number;
    // For storing the root mean squared error.
    rmsd: number;
    // For storing the weighted root mean squared error.
    wrmsd: number;
    // For storing the rotation matrix.
    rotmat: number[];
    // For storing the weights.
    weight: number[] | undefined;
    // The coordinates to be aligned.
    newX: number[];
    newY: number[];
    newZ: number[];
    // The coordinates to align with.
    x: number[];
    y: number[];
    z: number[];
    // The number of points.
    len: number;

    constructor(
        newX: number[],
        newY: number[],
        newZ: number[],
        x: number[],
        y: number[],
        z: number[]) {
        this.rotmat = [],
        this.newX = newX;
        this.newY = newY;
        this.newZ = newZ;
        this.x = x;
        this.y = y;
        this.z = z;
        this.len = newX.length;
        //console.log('newX.length:', newX.length);
        //console.log('newY.length:', newY.length);
        //console.log('newZ.length:', newZ.length);
        //console.log('x.length:', x.length);
        //console.log('y.length:', y.length);
        //console.log('z.length:', z.length);
        if (this.len !== newY.length || this.len !== newZ.length ||
            this.len !== x.length || this.len !== y.length || 
            this.len !== z.length) {
            throw new Error('Coordinate arrays must have the same length');
        }
        let A: number[] = [];
        let E0: number = this.innerProduct(A);
        this.fastCalcRMSDAndRotation(A, E0, -1);
    }

    getRotatedCoordinates(): { xR: number[]; yR: number[]; zR: number[] } {
        let xR: number[] = [];
        let yR: number[] = [];
        let zR: number[] = [];
        for (let i = 0; i < this.len; i++) {
            let x = this.newX[i];
            let y = this.newY[i];
            let z = this.newZ[i];
            let xRot = this.rotmat[0] * x + this.rotmat[1] * y + this.rotmat[2] * z;
            let yRot = this.rotmat[3] * x + this.rotmat[4] * y + this.rotmat[5] * z;
            let zRot = this.rotmat[6] * x + this.rotmat[7] * y + this.rotmat[8] * z;
            xR.push(xRot);
            yR.push(yRot);
            zR.push(zRot);
        }
        // Calculate euclidean distance between original and rotated coordinates
        let sdist = 0.0;
        let swdist = 0.0;
        let wsum = 0.0;
        for (let i = 0; i < this.len; i++) {
            wsum += this.weight![i];
            let dx = this.x[i] - xR[i];
            let dy = this.y[i] - yR[i];
            let dz = this.z[i] - zR[i];
            let dist = dx * dx + dy * dy + dz * dz;
            swdist += this.weight![i] * dist;
            sdist += dist;
        }
        this.wrmsd = Math.sqrt(swdist / wsum);
        this.rmsd = Math.sqrt(sdist / this.len);
        return { xR, yR, zR };
    }

    fastCalcRMSDAndRotation(A : number[], E0 : number, minScore: number): number {
        let Sxx: number;
        let Sxy: number;
        let Sxz: number;
        let Syx: number;
        let Syy: number;
        let Syz: number;
        let Szx: number;
        let Szy: number;
        let Szz: number;
        let Szz2: number;
        let Syy2: number;
        let Sxx2: number;
        let Sxy2: number;
        let Syz2: number;
        let Sxz2: number;
        let Syx2: number;
        let Szy2: number;
        let Szx2: number;
        let SyzSzymSyySzz2: number;
        let Sxx2Syy2Szz2Syz2Szy2: number;
        let Sxy2Sxz2Syx2Szx2: number;
        let SxzpSzx: number;
        let SyzpSzy: number;
        let SxypSyx: number;
        let SyzmSzy: number;
        let SxzmSzx: number;
        let SxymSyx: number;
        let SxxpSyy: number;
        let SxxmSyy: number;
        let C: number[] = [];
        let i: number;
        let mxEigenV: number;
        let oldg: number;
        let b: number;
        let a: number;
        let delta: number;
        let qsqr: number;
        let normq: number;
        let a11: number;
        let a12: number;
        let a13: number;
        let a14: number;
        let a21: number;
        let a22: number;
        let a23: number;
        let a24: number;
        let a31: number;
        let a32: number;
        let a33: number;
        let a34: number;
        let a41: number;
        let a42: number;
        let a43: number;
        let a44: number;
        let a2: number;
        let x2: number;
        let y2: number;
        let z2: number;
        let xy: number;
        let az: number;
        let zx: number;
        let ay: number;
        let yz: number;
        let ax: number;
        let a3344_4334: number;
        let a3244_4234: number;
        let a3243_4233: number;
        let a3143_4133: number;
        let a3144_4134: number;
        let a3142_4132: number;
        let evecprec: number = 1.0e-6;
        let evalprec: number = 1.0e-11;
        
        Sxx = A[0];
        Sxy = A[1];
        Sxz = A[2];
        Syx = A[3];
        Syy = A[4];
        Syz = A[5];
        Szx = A[6];
        Szy = A[7];
        Szz = A[8];

        Sxx2 = Sxx * Sxx;
        Syy2 = Syy * Syy;
        Szz2 = Szz * Szz;

        Sxy2 = Sxy * Sxy;
        Syz2 = Syz * Syz;
        Sxz2 = Sxz * Sxz;

        Syx2 = Syx * Syx;
        Szy2 = Szy * Szy;
        Szx2 = Szx * Szx;

        SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz);
        Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

        C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
        C[1] = 8.0 * (Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx - Sxx * Syy * Szz - Syz * Szx * Sxy - Szy * Syx * Sxz);

        SxzpSzx = Sxz + Szx;
        SyzpSzy = Syz + Szy;
        SxypSyx = Sxy + Syx;
        SyzmSzy = Syz - Szy;
        SxzmSzx = Sxz - Szx;
        SxymSyx = Sxy - Syx;
        SxxpSyy = Sxx + Syy;
        SxxmSyy = Sxx - Syy;
        Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

        C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
                + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
                + (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz)) * (-(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz))
                + (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz)) * (-(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz))
                + (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz)) * (-(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz))
                + (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz)) * (-(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz));

        /* Newton-Raphson */
        mxEigenV = E0;
        for (i = 0; i < 50; i++) {
            oldg = mxEigenV;
            x2 = mxEigenV * mxEigenV;
            b = (x2 + C[2]) * mxEigenV;
            a = b + C[1];
            delta = ((a * mxEigenV + C[0]) / (2.0 * x2 * mxEigenV + b + a));
            mxEigenV -= delta;
            //rmsd = Math.sqrt(Math.abs(2d * (E0 - mxEigenV) / len));
            //System.out.println("Iteration " + i + " rmsd " + rmsd);
            if (Math.abs(mxEigenV - oldg) < Math.abs(evalprec * mxEigenV)) {
                break;
            }
        }

        if (i == 50) {
            console.log("More than " + i + " iterations needed!");
        }
        /* the Math.abs() is to guard against extremely small, but *negative* numbers due to floating point error */
        let score: number = Math.sqrt(Math.abs(2.0 * (E0 - mxEigenV) / this.len));

        //System.out.println("rmsd " + rmsd);
        if (minScore > 0) {
            if (score < minScore) {
                return -1; // Don't bother with rotation.
            }
        }

        a11 = SxxpSyy + Szz - mxEigenV;
        a12 = SyzmSzy;
        a13 = -SxzmSzx;
        a14 = SxymSyx;
        a21 = SyzmSzy;
        a22 = SxxmSyy - Szz - mxEigenV;
        a23 = SxypSyx;
        a24 = SxzpSzx;
        a31 = a13;
        a32 = a23;
        a33 = Syy - Sxx - Szz - mxEigenV;
        a34 = SyzpSzy;
        a41 = a14;
        a42 = a24;
        a43 = a34;
        a44 = Szz - SxxpSyy - mxEigenV;
        a3344_4334 = a33 * a44 - a43 * a34;
        a3244_4234 = a32 * a44 - a42 * a34;
        a3243_4233 = a32 * a43 - a42 * a33;
        a3143_4133 = a31 * a43 - a41 * a33;
        a3144_4134 = a31 * a44 - a41 * a34;
        a3142_4132 = a31 * a42 - a41 * a32;
        this.q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
        this.q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
        this.q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
        this.q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;

        qsqr = this.q1 * this.q1 + this.q2 * this.q2 + this.q3 * this.q3 + this.q4 * this.q4;

        /**
         * The following code tries to calculate another column in the adjoint
         * matrix when the norm of the current column is too small. Usually this
         * block will never be activated. To be absolutely safe this should be
         * uncommented, but it is most likely unnecessary.
         */
        if (qsqr < evecprec) {
            this.q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
            this.q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
            this.q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
            this.q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;
            qsqr = this.q1 * this.q1 + this.q2 * this.q2 + this.q3 * this.q3 + this.q4 * this.q4;

            if (qsqr < evecprec) {
                let a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
                let a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
                let a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

                this.q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
                this.q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
                this.q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
                this.q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
                qsqr = this.q1 * this.q1 + this.q2 * this.q2 + this.q3 * this.q3 + this.q4 * this.q4;

                if (qsqr < evecprec) {
                    this.q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                    this.q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                    this.q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                    this.q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                    qsqr = this.q1 * this.q1 + this.q2 * this.q2 + this.q3 * this.q3 + this.q4 * this.q4;

                    if (qsqr < evecprec) {
                        /* if qsqr is still too small, return the identity matrix. */
                        this.rotmat[0] = this.rotmat[4] = this.rotmat[8] = 1.0;
                        this.rotmat[1] = this.rotmat[2] = this.rotmat[3] = this.rotmat[5] = this.rotmat[6] = this.rotmat[7] = 0.0;

                        return 0;
                    }
                }
            }
        }

        normq = Math.sqrt(qsqr);
        this.q1 /= normq;
        this.q2 /= normq;
        this.q3 /= normq;
        this.q4 /= normq;

        a2 = this.q1 * this.q1;
        x2 = this.q2 * this.q2;
        y2 = this.q3 * this.q3;
        z2 = this.q4 * this.q4;

        xy = this.q2 * this.q3;
        az = this.q1 * this.q4;
        zx = this.q4 * this.q2;
        ay = this.q1 * this.q3;
        yz = this.q3 * this.q4;
        ax = this.q1 * this.q2;

        this.rotmat[0] = a2 + x2 - y2 - z2;
        this.rotmat[1] = 2 * (xy + az);
        this.rotmat[2] = 2 * (zx - ay);
        this.rotmat[3] = 2 * (xy - az);
        this.rotmat[4] = a2 - x2 + y2 - z2;
        this.rotmat[5] = 2 * (yz + ax);
        this.rotmat[6] = 2 * (zx + ay);
        this.rotmat[7] = 2 * (yz - ax);
        this.rotmat[8] = a2 - x2 - y2 + z2;

        return 1;
    }

    innerProduct(A: number[]): number {
        let x1: number;
        let x2: number;
        let y1: number;
        let y2: number;
        let z1: number;
        let z2: number;
        let fx1: number[] = this.x;
        let fy1: number[] = this.y;
        let fz1: number[] = this.z;
        let fx2: number[] = this.newX;
        let fy2: number[] = this.newY;
        let fz2: number[] = this.newZ;
        let G1 = 0.0;
        let G2 = 0.0;
        A[0] = 0.0;
        A[1] = 0.0;
        A[2] = 0.0;
        A[3] = 0.0;
        A[4] = 0.0;
        A[5] = 0.0;
        A[6] = 0.0;
        A[7] = 0.0;
        A[8] = 0.0;
        if (!this.weight) {
            for (let i = 0; i < this.len; i++) {
                x1 = fx1[i];
                y1 = fy1[i];
                z1 = fz1[i];
                G1 += x1 * x1 + y1 * y1 + z1 * z1;
                x2 = fx2[i];
                y2 = fy2[i];
                z2 = fz2[i];
                G2 += x2 * x2 + y2 * y2 + z2 * z2;
                A[0] += (x1 * x2);
                A[1] += (x1 * y2);
                A[2] += (x1 * z2);
                A[3] += (y1 * x2);
                A[4] += (y1 * y2);
                A[5] += (y1 * z2);
                A[6] += (z1 * x2);
                A[7] += (z1 * y2);
                A[8] += (z1 * z2);
            }
            this.weight = [];
            for (let i = 0; i < this.len; i++) {
                this.weight.push(i + 1.0);
            }
        } else {
            for (let i = 0; i < this.len; i++) {
                x1 = this.weight[i] * fx1[i];
                y1 = this.weight[i] * fy1[i];
                z1 = this.weight[i] * fz1[i];
                G1 += x1 * fx1[i] + y1 * fy1[i] + z1 * fz1[i];
                x2 = fx2[i];
                y2 = fy2[i];
                z2 = fz2[i];
                G2 += this.weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);
                A[0] += (x1 * x2);
                A[1] += (x1 * y2);
                A[2] += (x1 * z2);
                A[3] += (y1 * x2);
                A[4] += (y1 * y2);
                A[5] += (y1 * z2);
                A[6] += (z1 * x2);
                A[7] += (z1 * y2);
                A[8] += (z1 * z2);
            }
        }
        return (G1 + G2) / 2.0;
    }

}