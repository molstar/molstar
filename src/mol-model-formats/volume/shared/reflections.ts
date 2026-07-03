/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from '../../../mol-task';
import { Spacegroup } from '../../../mol-math/geometry';
import { fft3dHermitian, nextPow2 } from '../../../mol-math/fourier-transform/fft';

/**
 * Core electron density computation shared by all structure-factor → volume pathways.
 *
 * Given arrays of Miller indices (h, k, l), structure-factor amplitudes, and phases
 * (in degrees), this function:
 *   1. Sizes a reciprocal-space grid (next power-of-2 that accommodates all indices),
 *   2. Expands reflections through all spacegroup symmetry operators,
 *   3. Applies Friedel's law,
 *   4. Computes the 3D forward FFT using the Hermitian half-spectrum shortcut,
 *   5. Returns the real-valued density grid and its dimensions.
 *
 * @param h         Miller index h array (length = number of reflections)
 * @param k         Miller index k array
 * @param l         Miller index l array
 * @param amp       Structure-factor amplitudes (same length)
 * @param phi       Structure-factor phases in degrees (same length)
 * @param sg        Spacegroup (symmetry operators in fractional coordinates)
 * @param ctx       RuntimeContext for progress reporting
 * @returns         Density Float32Array in row-major (N0, N1, N2) order plus grid dims
 */
export async function computeElectronDensityFromReflections(
    h: ArrayLike<number>,
    k: ArrayLike<number>,
    l: ArrayLike<number>,
    amp: ArrayLike<number>,
    phi: ArrayLike<number>,
    sg: Spacegroup,
    ctx: RuntimeContext,
): Promise<{ density: Float32Array, N0: number, N1: number, N2: number }> {
    const count = h.length;

        // 3. Determine grid dimensions from reflection index ranges
        let maxH = 0, maxK = 0, maxL = 0;
        for (let i = 0; i < count; i++) {
            const ah = Math.abs(h[i]), ak = Math.abs(k[i]), al = Math.abs(l[i]);
            if (ah > maxH) maxH = ah;
            if (ak > maxK) maxK = ak;
            if (al > maxL) maxL = al;
        }

        // Grid must contain indices from -max to +max: minimum size = 2*max + 1
        // Round up to next power of 2 for FFT efficiency
        const N0 = nextPow2(2 * maxH + 2);
        const N1 = nextPow2(2 * maxK + 2);
        const N2 = nextPow2(2 * maxL + 2);

        // 4. Expand reflections using spacegroup operators and fill reciprocal-space grid
        await ctx.update({ message: 'Expanding structure factors...' });

        const totalSize = N0 * N1 * N2;
        // Only the non-redundant lower half-spectrum (l = 0..N2/2) is stored; the
        // upper planes follow from Hermitian symmetry and are never materialized,
        // which halves the reciprocal-space grid memory.
        const M = N2 >> 1;
        const P = M + 1;
        const halfSize = N0 * N1 * P;
        const gridRe = new Float32Array(halfSize);
        const gridIm = new Float32Array(halfSize);
        // Track which grid positions have been filled (avoids overwriting with lower-occupancy data)
        const filled = new Uint8Array(halfSize);

        const operators = sg.operators;

        /**
         * Place a single structure factor F = amp * exp(i * phiDeg * π/180)
         * at grid index (rh, rk, rl) using wrap-around for negative indices.
         * Only the lower half-spectrum (rl = 0..N2/2) is kept: a reflection in the
         * upper half is folded onto its Hermitian mate (-rh,-rk,-rl) as conj(F).
         */
        const setRefl = (rh: number, rk: number, rl: number, rAmp: number, phiDeg: number) => {
            let hi = ((rh % N0) + N0) % N0;
            let ki = ((rk % N1) + N1) % N1;
            let li = ((rl % N2) + N2) % N2;
            let phiRad = phiDeg * (Math.PI / 180);
            if (li > M) {
                // Fold onto the Hermitian mate (-rh,-rk,-rl), storing conj(F).
                hi = (N0 - hi) % N0;
                ki = (N1 - ki) % N1;
                li = N2 - li;
                phiRad = -phiRad;
            }
            const idx = hi * N1 * P + ki * P + li;
            if (filled[idx]) return;
            gridRe[idx] = rAmp * Math.cos(phiRad);
            gridIm[idx] = rAmp * Math.sin(phiRad);
            filled[idx] = 1;
        };

        for (let i = 0; i < count; i++) {
            const hVal = h[i], kVal = k[i], lVal = l[i];
            const aVal = amp[i], phiDeg = phi[i];
            if (!isFinite(aVal) || !isFinite(phiDeg)) continue;

            // Apply each spacegroup symmetry operator to generate equivalent reflections.
            // For real-space operator M (fractional): r' = R*r + t
            // The reciprocal-space equivalent: h' = R^T * h
            // Phase correction: δφ = 360° * (h · t)
            // Mat4 is column-major: element(row i, col j) = op[j*4+i]
            for (let oi = 0; oi < operators.length; oi++) {
                const op = operators[oi];
                // Rotation part R^T applied to (hVal, kVal, lVal):
                // h'[0] = Σ_j R[j,0]*h[j]  = op[0]*h + op[1]*k + op[2]*l
                // h'[1] = Σ_j R[j,1]*h[j]  = op[4]*h + op[5]*k + op[6]*l
                // h'[2] = Σ_j R[j,2]*h[j]  = op[8]*h + op[9]*k + op[10]*l
                const hp = Math.round(op[0] * hVal + op[1] * kVal + op[2] * lVal);
                const kp = Math.round(op[4] * hVal + op[5] * kVal + op[6] * lVal);
                const lp = Math.round(op[8] * hVal + op[9] * kVal + op[10] * lVal);
                // Translation in fractional coordinates: t = (op[12], op[13], op[14])
                const dphiDeg = 360 * (hVal * op[12] + kVal * op[13] + lVal * op[14]);
                const newPhi = phiDeg + dphiDeg;
                setRefl(hp, kp, lp, aVal, newPhi);
                // Friedel mate: F(-h') = aVal * exp(-i*newPhi) = conj(F(h'))
                setRefl(-hp, -kp, -lp, aVal, -newPhi);
            }
        }

        // 5. Compute 3D forward FFT
        // The forward DFT gives: ρ[n0,n1,n2] = Σ_{h,k,l} F[h,k,l] * exp(-2πi*(h*n0/N0+k*n1/N1+l*n2/N2))
        // which is exactly the electron density formula. Because the structure
        // factors are Hermitian the density is real; we pass only the stored lower
        // half-spectrum to the specialized transform, which writes the full density.
        await ctx.update({ message: 'Computing Fourier transform...' });
        const density = new Float32Array(totalSize);
        fft3dHermitian(gridRe, gridIm, N0, N1, N2, density);

    return { density, N0, N1, N2 };
}