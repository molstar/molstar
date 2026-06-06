/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

import { parseKin } from '../reader/parser';

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

// Complex kinemage with multiple features: animate groups, pointmasters, various list types
const kinComplexString = `@kinemage 1
@caption Complex test kinemage with multiple features
@text
This is a comprehensive test kinemage file that includes:
- Multiple groups with animate and 2animate
- Pointmasters with tags
- All list types: dots, vectors, balls, spheres, ribbons, triangles
@master {main} on
@master {secondary} off
@master {alternate}
@pointmaster 'ABC' {Primary atoms} on
@pointmaster 'XY' {Secondary atoms} off
@group {Structure} animate dominant
@subgroup {Backbone}
@vectorlist {CA trace} color=blue master={main}
{CA ALA 1}blue P 'A' 10.0,20.0,30.0 {CA ALA 2}blue 'A' 11.0,21.0,31.0
{"}blue 'A' 12.0,22.0,32.0 {CA ALA 3}blue 'A' 13.0,23.0,33.0
@dotlist {H-bonds} color=yellow master={main}
{HN ALA 2}yellow 'B' 10.5,20.5,30.5
{HN ALA 3}yellow 'B' 11.5,21.5,31.5
{HN ALA 4}yellow 'B' 12.5,22.5,32.5
@subgroup {Sidechains}
@balllist {CB atoms} color=green master={secondary} radius=0.5
{CB ARG 1}green r=0.5 'C' 9.0,19.0,29.0
{CB ARG 2}green r=0.6 'C' 10.0,20.0,30.0
{CB ARG 3}green r=0.55 'C' 11.0,21.0,31.0
@group {Alternate conformations} 2animate
@subgroup {Alt A}
@spherelist {Waters A} color=cyan master={alternate} radius=1.0
{HOH 101}cyan r=1.0 'X' 15.0,25.0,35.0
{HOH 102}cyan r=1.2 'X' 16.0,26.0,36.0
@subgroup {Alt B}
@spherelist {Waters B} color=magenta master={alternate} radius=1.0
{HOH 101}magenta r=1.0 'Y' 15.2,25.2,35.2
{HOH 102}magenta r=1.1 'Y' 16.1,26.1,36.1
@group {Surface} off
@subgroup {Ribbons}
@ribbonlist {Alpha helix} color=red master={main}
{ASP 5}red 14.0,24.0,34.0
{GLU 6}red 15.0,25.0,35.0
{LYS 7}red 16.0,26.0,36.0
{ARG 8}red 17.0,27.0,37.0
{THR 9}red P 18.0,28.0,38.0
{VAL 10}red 19.0,29.0,39.0
@subgroup {Triangles}
@trianglelist {Surface patch} color=sky master={secondary}
{Tri 1}sky 20.0,30.0,40.0
{Tri 1}sky 21.0,30.0,40.0
{Tri 1}sky 20.5,31.0,40.0
{Tri 2}sky X 22.0,32.0,42.0
{Tri 2}sky 23.0,32.0,42.0
{Tri 2}sky 22.5,33.0,42.0
@group {Contacts} animate
@subgroup {Clashes}
@vectorlist {Bad overlaps} color=hotpink master={main} width=4
{O HOH 319 A}hotpink P 31.146,32.100,-1.425 {O HOH 320 A}hotpink 31.015,32.234,-1.324
{"}hotpink P 31.607,32.750,-1.156 {"}hotpink 31.410,32.784,-1.097
`;

describe('kin reader', () => {
    it('basic', async () => {
        const parsed = await parseKin(kinString).run();
        if (parsed.isError) {
            console.error('Parse error:', parsed);
            fail('Parse should not error');
        }
        if (parsed.result.length !== 1) {
            fail(`Expected 1 kinemage, got ${parsed.result.length}`);
        }
        const kinemage = parsed.result[0];

        const vectors = kinemage.vectorLists;
        expect(vectors.length).toEqual(1);

        const element = vectors[0];
        expect(element.name).toEqual('x');
        expect(element.position1Array.length).toEqual(7*3);

        // Test that colors are parsed correctly
        expect(element.color1Array.length).toEqual(7);

        // Test masters are set up
        expect(element.masterArray).toContain('bad overlap');

        expect.assertions(5);
    });

    it('complex', async () => {
        const parsed = await parseKin(kinComplexString).run();
        if (parsed.isError) {
            fail('Parse should not error');
        }

        expect(parsed.result.length).toBeGreaterThan(0);
        const kinemage = parsed.result[0];

        // Verify structure is valid
        expect(kinemage.vectorLists).toBeDefined();
        expect(kinemage.masterDict).toBeDefined();
        expect(kinemage.groupDict).toBeDefined();
        expect(kinemage.pointmasterDict).toBeDefined();

        // Test animate groups
        expect(kinemage.groupsAnimate.length).toEqual(2);
        expect(kinemage.groupsAnimate).toContain('Structure');
        expect(kinemage.groupsAnimate).toContain('Contacts');
        expect(kinemage.activeAnimateGroup).toEqual(0);

        // Test 2animate groups
        expect(kinemage.groupsAnimate2.length).toEqual(1);
        expect(kinemage.groupsAnimate2).toContain('Alternate conformations');
        expect(kinemage.activeAnimateGroup2).toEqual(0);

        // Test pointmasters
        expect(Object.keys(kinemage.pointmasterDict).length).toBeGreaterThan(0);
        expect(kinemage.pointmasterDict['A']).toEqual('Primary atoms');
        expect(kinemage.pointmasterDict['B']).toEqual('Primary atoms');
        expect(kinemage.pointmasterDict['X']).toEqual('Secondary atoms');

        // Test masters
        expect(kinemage.masterDict['main']).toBeDefined();
        expect(kinemage.masterDict['main'].visible).toEqual(true);
        expect(kinemage.masterDict['secondary']).toBeDefined();
        expect(kinemage.masterDict['secondary'].visible).toEqual(false);

        // Test list types
        expect(kinemage.vectorLists.length).toEqual(2);
        expect(kinemage.dotLists.length).toEqual(1);
        expect(kinemage.ballLists.length).toEqual(3); // 1 balllist + 2 spherelists
        expect(kinemage.ribbonLists.length).toEqual(2); // 1 ribbonlist + 1 trianglelist

        // Test specific list properties
        const caTrace = kinemage.vectorLists.find(v => v.name === 'CA trace');
        expect(caTrace).toBeDefined();
        expect(caTrace?.masterArray).toContain('main');

        const hBonds = kinemage.dotLists[0];
        expect(hBonds.name).toEqual('H-bonds');
        expect(hBonds.positionArray.length).toEqual(9); // 3 dots * 3 coords

        const cbAtoms = kinemage.ballLists.find(b => b.name === 'CB atoms');
        expect(cbAtoms).toBeDefined();
        expect(cbAtoms?.radiusArray.length).toEqual(3);

        const helix = kinemage.ribbonLists.find(r => r.name === 'Alpha helix');
        expect(helix).toBeDefined();
        expect(helix?.pairTriangleNormals).toEqual(true); // ribbonlist

        const surface = kinemage.ribbonLists.find(r => r.name === 'Surface patch');
        expect(surface).toBeDefined();
        expect(surface?.pairTriangleNormals).toEqual(false); // trianglelist

        // Test groups
        expect(Object.keys(kinemage.groupDict).length).toEqual(4);
        expect(kinemage.groupDict['Structure'].animate).toEqual(true);
        expect(kinemage.groupDict['Alternate conformations']['2animate']).toEqual(true);
        expect(kinemage.groupDict['Surface'].off).toEqual(true);

        expect.assertions(38);
    });

});