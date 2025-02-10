/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

import { parseKin } from '../kin/parser';

const kinString = `@kinemage 1
@caption probe.2.26.021123, run Tue Apr 23 14:49:17 2024
 command: C:\tmp\cctbx_phenix\build\probe\exe\probe.exe -kin -mc -het -once -wat2wat -onlybadout -stdbonds water all 1ssxFH.pdb
@group dominant {dots}
@subgroup dominant {once dots}
@master {bad overlap}
@pointmaster 'O' {Hets contacts}
@vectorlist {x} color=red master={bad overlap}
{ O   HOH 319  A}hotpink P  'O' 31.146,32.100,-1.425 {"}hotpink   'O' 31.015,32.234,-1.324
{"}hotpink P  'O' 31.607,32.750,-1.156 {"}hotpink   'O' 31.410,32.784,-1.097
{"}hotpink P  'O' 31.263,32.074,-1.185 {"}hotpink   'O' 31.117,32.209,-1.122
{ O  BHOH 338  A}hotpink P  'O' 32.540,45.631,10.833 {"}hotpink   'O' 32.430,45.771,10.977
{"}hotpink P  'O' 32.316,45.500,10.828 {"}hotpink   'O' 32.230,45.689,10.998
{"}hotpink P  'O' 32.068,45.424,10.824 {"}hotpink   'O' 32.034,45.604,10.975
{"}hotpink P  'O' 32.729,45.605,11.052 {"}hotpink   'O' 32.572,45.765,11.173
`;

/// @todo Replace with more complex kinemage
const kinComplexString = `@kinemage 1
@caption probe.2.26.021123, run Tue Apr 23 14:49:17 2024
 command: C:\tmp\cctbx_phenix\build\probe\exe\probe.exe -kin -mc -het -once -wat2wat -onlybadout -stdbonds water all 1ssxFH.pdb
@group dominant {dots}
@subgroup dominant {once dots}
@master {bad overlap}
@pointmaster 'O' {Hets contacts}
@vectorlist {x} color=red master={bad overlap}
{ O   HOH 319  A}hotpink P  'O' 31.146,32.100,-1.425 {"}hotpink   'O' 31.015,32.234,-1.324
{"}hotpink P  'O' 31.607,32.750,-1.156 {"}hotpink   'O' 31.410,32.784,-1.097
{"}hotpink P  'O' 31.263,32.074,-1.185 {"}hotpink   'O' 31.117,32.209,-1.122
{ O  BHOH 338  A}hotpink P  'O' 32.540,45.631,10.833 {"}hotpink   'O' 32.430,45.771,10.977
{"}hotpink P  'O' 32.316,45.500,10.828 {"}hotpink   'O' 32.230,45.689,10.998
{"}hotpink P  'O' 32.068,45.424,10.824 {"}hotpink   'O' 32.034,45.604,10.975
{"}hotpink P  'O' 32.729,45.605,11.052 {"}hotpink   'O' 32.572,45.765,11.173
`;

describe('kin reader', () => {
    it('basic', async () => {
        const parsed = await parseKin(kinString).run();
        if (parsed.isError) return;
        const kinemage = parsed.result;

        const vectors = kinemage.vectorLists;
        expect(vectors.length).toEqual(1);

        const element = vectors[0];
        expect(element.name).toEqual('x');
        expect(element.position1Array.length).toEqual(7);

        /// @todo Add more tests

        expect.assertions(3);
    });

    it('complex', async () => {
        const parsed = await parseKin(kinComplexString).run();
        if (parsed.isError) return;

        /// @todo Add more complex tests

    });
});