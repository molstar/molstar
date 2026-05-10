/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * Based on earlier kin-parser.ts file from the NGL project (see second author notice below).
 * @file Ported NGL-based Kinemage file parser
 * @author ReliaSolve <russ@reliasolve.com>
 * @private
 */

/**
 * file Kin Parser
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Kinemage, RibbonObject } from './schema';
import { Hsv } from '../../../mol-util/color/spaces/hsv';
import { Color } from '../../../mol-util/color';

const ColorDict: { [k: string]: Color } = {
    red: Hsv.toColor(Hsv.fromArray([0, 100, 100])),
    orange: Hsv.toColor(Hsv.fromArray([20, 100, 100])),
    gold: Hsv.toColor(Hsv.fromArray([40, 100, 100])),
    yellow: Hsv.toColor(Hsv.fromArray([60, 100, 100])),
    lime: Hsv.toColor(Hsv.fromArray([80, 100, 100])),
    green: Hsv.toColor(Hsv.fromArray([120, 80, 100])),
    sea: Hsv.toColor(Hsv.fromArray([150, 100, 100])),
    cyan: Hsv.toColor(Hsv.fromArray([180, 100, 85])),
    sky: Hsv.toColor(Hsv.fromArray([210, 75, 95])),
    blue: Hsv.toColor(Hsv.fromArray([240, 70, 100])),
    purple: Hsv.toColor(Hsv.fromArray([275, 75, 100])),
    magenta: Hsv.toColor(Hsv.fromArray([300, 95, 100])),
    hotpink: Hsv.toColor(Hsv.fromArray([335, 100, 100])),
    pink: Hsv.toColor(Hsv.fromArray([350, 55, 100])),
    peach: Hsv.toColor(Hsv.fromArray([25, 75, 100])),
    lilac: Hsv.toColor(Hsv.fromArray([275, 55, 100])),
    pinktint: Hsv.toColor(Hsv.fromArray([340, 30, 100])),
    peachtint: Hsv.toColor(Hsv.fromArray([25, 50, 100])),
    yellowtint: Hsv.toColor(Hsv.fromArray([60, 50, 100])),
    greentint: Hsv.toColor(Hsv.fromArray([135, 40, 100])),
    bluetint: Hsv.toColor(Hsv.fromArray([220, 40, 100])),
    lilactint: Hsv.toColor(Hsv.fromArray([275, 35, 100])),
    white: Hsv.toColor(Hsv.fromArray([0, 0, 100])),
    gray: Hsv.toColor(Hsv.fromArray([0, 0, 50])),
    brown: Hsv.toColor(Hsv.fromArray([20, 45, 75])),
    deadwhite: Hsv.toColor(Hsv.fromArray([0, 0, 100])),
    deadblack: Hsv.toColor(Hsv.fromArray([0, 0, 0])),
    invisible: Hsv.toColor(Hsv.fromArray([0, 0, 0]))
};

const reWhitespaceComma = /[\s,]+/;
const reCurlyWhitespace = /[^{}\s]*{[^{}]+}|[^{}\s]+/g;
const reTrimCurly = /^{+|}+$/g;
const reCollapseEqual = /\s*=\s*/g;

function parseListDef(line: string, localColorDict: { [k: string]: Color }) {
    let name;
    let defaultColor: Color = localColorDict['white']; // Default color is white, but it can be overridden by the list definition
    let radius;
    let nobutton = false;
    const master = [];
    let width = 2; // Default width is 2, but it can be overridden by the list definition

    line = line.replace(reCollapseEqual, '=');

    const lm = line.match(reCurlyWhitespace) as RegExpMatchArray;
    for (let j = 1; j < lm.length; ++j) {
        const e = lm[j];
        if (e[0] === '{') {
            name = e.substring(1, e.length - 1);
        } else {
            const es = e.split('=');
            if (es.length === 2) {
                if (es[0] === 'color') {
                    const colorName = parseStr(es[1]);
                    defaultColor = localColorDict[colorName];
                } else if (es[0] === 'width') {
                    width = parseInt(es[1]);
                } else if (es[0] === 'master') {
                    master.push(es[1].replace(reTrimCurly, ''));
                } else if (es[0] === 'radius') {
                    radius = parseFloat(es[1]);
                } else {
                    console.log('Kinemage: Unknown list definition term found: ' + es[0]);
                }
            } else if (e === 'nobutton') {
                nobutton = true;
            } else {
                console.log('Kinemage: Unknown list definition term found: ' + e);
            }
        }
    }

    return {
        listName: name,
        listColor: defaultColor,
        listMasters: master,
        listWidth: width,
        listRadius: radius,
        nobutton: nobutton
    };
}

function parseListElm(line: string, localColorDict: { [k: string]: Color }) {
    line = line.trim();

    const idx1 = line.indexOf('{');
    const idx2 = line.indexOf('}');
    const ls = line.substr(idx2 + 1).split(reWhitespaceComma);

    const label = line.substr(idx1 + 1, idx2 - 1);
    const position = [
        parseFloat(ls[ls.length - 3]),
        parseFloat(ls[ls.length - 2]),
        parseFloat(ls[ls.length - 1])
    ];
    let color, width, radius;
    let lineBreak = false;
    let triangleBreak = false;
    const pointMasters: string[] = [];
    for (let lsindex = 4; lsindex <= ls.length; lsindex++) {
        const literal = ls[ls.length - lsindex];
        if (literal in localColorDict) {
            color = localColorDict[ls[ls.length - lsindex]];
        }
        if (literal.startsWith('width')) {
            width = parseInt(literal.substring(5));
        }
        if (literal.startsWith('r=')) {
            radius = parseFloat(literal.split('=')[1]);
        }
        if (literal.startsWith('P')) {
            lineBreak = true;
        }
        if (literal.startsWith('X')) {
            triangleBreak = true;
        }
        if (literal.startsWith("'") && literal.endsWith("'")) {
            // Handle single-character tags by putting each character into a pointMaster tag, e.g. 'ab' would be two tags, 'a' and 'b'
            const tagString: string = literal.substring(1, literal.length - 1);
            for (let i = 0; i < tagString.length; i++) {
                pointMasters.push(tagString[i]);
            }
        }
    }

    return {
        label: label,
        position: position,
        color: color,
        radius: radius,
        width: width,
        isLineBreak: lineBreak,
        isTriangleBreak: triangleBreak,
        pointMasters: pointMasters
    };
}

function parseStr(line: string) {
    const start = line.indexOf('{');
    const end = line.indexOf('}');
    return line.substring(
        start !== -1 ? start + 1 : 0,
        end !== -1 ? end : undefined
    ).trim();
}

function parseFlag(line: string) {
    const end = line.indexOf('}');
    return end === -1 ? undefined : line.substr(end + 1).trim();
}

function parseGroup(line: string) {
    let name: string = '';
    const master: string[] = [];
    const flags: { [k: string]: string | boolean } = {};

    line = line.replace(reCollapseEqual, '=');

    const lm = line.match(reCurlyWhitespace) as RegExpMatchArray;
    for (let j = 1; j < lm.length; ++j) {
        const e = lm[j];
        if (e[0] === '{') {
            name = e.substring(1, e.length - 1);
        } else {
            const es = e.split('=');
            if (es.length === 2) {
                if (es[0] === 'master') {
                    master.push(es[1].replace(reTrimCurly, ''));
                } else {
                    flags[es[0]] = es[1].replace(reTrimCurly, '');
                }
            } else {
                flags[es[0]] = true;
            }
        }
    }

    return {
        groupName: name,
        groupFlags: flags,
        groupMasters: master,
    };
}

function parsePointmaster(line: string) {
    let name: string = '';
    const tags: string[] = [];
    let on: boolean | undefined = undefined;

    // Find the string name between curly braces, or print an error if not found
    const nameMatch = line.match(/{([^}]+)}/);
    if (nameMatch) {
        name = nameMatch[1];

        // Find all characters between the pair of single quotes, which are the tags, and add them to the tags array
        const tagMatch = line.match(/'([^']+)'/);
        if (tagMatch) {
            const tagString: string = tagMatch[1];
            for (let i = 0; i < tagString.length; i++) {
                tags.push(tagString[i]);
            }

            // See if the line contains the word "on" or "off" and set the on variable accordingly
            if (line.includes(' on')) {
                on = true;
            } else if (line.includes(' off')) {
                on = false;
            }

        } else {
            console.log('Kinemage: Pointmaster definition missing tags: ' + line);
        }
    } else {
        console.log('Kinemage: Pointmaster definition missing name: ' + line);
    }

    return {
        name: name,
        tags: tags,
        on: on
    };
}

function convertKinTriangleArrays(ribbonObject: RibbonObject) {
    // have to convert ribbons/triangle lists from stripdrawmode to normal drawmode
    // index           [ 0 1 2 3 4 5 6 7 8 91011 ]
    // label/color/ptm [ 0 1 2 3 4 5 ] to [ 0 1 2 1 2 3 2 3 4 3 4 5 ]
    // convertedindex                                       [ 0 1 2 3 4 5 6 7 8 91011121314151617181920212223242526 ]
    // index           [ 0 1 2 3 4 5 6 7 8 91011121314 ]    [ 0 1 2 3 4 5 6 7 8 3 4 5 6 7 8 91011 6 7 8 91011121314 ]
    // position        [ 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 ] to [ 0 0 0 1 1 1 2 2 2 1 1 1 2 2 2 3 3 3 2 2 2 3 3 3 4 4 4 ]
    const { labelArray, positionArray, colorArray, breakArray } = ribbonObject;
    const convertedLabels = [];
    for (let i = 0; i < (labelArray.length - 2) * 3; ++i) {
        convertedLabels[i] = labelArray[i - Math.floor(i / 3) * 2];
    }
    const convertedColors = [];
    for (let i = 0; i < (colorArray.length - 2) * 3; ++i) {
        convertedColors[i] = colorArray[i - Math.floor(i / 3) * 2];
    }
    const convertedPMs = [];
    for (let i = 0; i < (ribbonObject.pointmasterArray.length - 2) * 3; ++i) {
        convertedPMs[i] = ribbonObject.pointmasterArray[i - Math.floor(i / 3) * 2];
    }
    const convertedBreaks = [];
    for (let i = 0; i < (breakArray.length - 2) * 3; ++i) {
        convertedBreaks[i] = breakArray[i - Math.floor(i / 3) * 2];
    }
    const convertedPositions = [];
    for (let i = 0; i < (positionArray.length / 3 - 2) * 9; ++i) {
        convertedPositions[i] = positionArray[i - Math.floor(i / 9) * 6];
    }
    const vector3Positions = [];
    for (let i = 0; i < (convertedPositions.length) / 3; ++i) {
        vector3Positions.push([convertedPositions[i * 3], convertedPositions[i * 3] + 1, convertedPositions[i * 3] + 2]);
    }
    return {
        group: ribbonObject.group,
        subgroup: ribbonObject.subgroup,
        name: ribbonObject.name,
        masterArray: ribbonObject.masterArray,
        pointmasterArray: convertedPMs,
        nobutton: ribbonObject.nobutton,
        labelArray: convertedLabels,
        positionArray: convertedPositions,
        breakArray: convertedBreaks,
        colorArray: convertedColors,
        pairTriangleNormals: ribbonObject.pairTriangleNormals
    };
}

function removePointBreaksTriangleArrays(convertedRibbonObject: RibbonObject) {
    // after converting ribbon/triangle arrys to drawmode, removed point break triangles
    // label/color [ 0 1 2 3 4 5 ] to [ 0 1 2 1 2 3 2 3 4 3 4 5 ]
    // position    [ 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 ] to [ 0 0 0 1 1 1 2 2 2 1 1 1 2 2 2 3 3 3 2 2 2 3 3 3 4 4 4 ]
    const { labelArray, positionArray, colorArray, breakArray } = convertedRibbonObject;
    const editedLabels = [];
    const editedPositions = [];
    const editedColors = [];
    const editedPMs = [];
    const editedBreaks = [];
    for (let i = 0; i < breakArray.length / 3; i++) {
        const breakPointer = i * 3;
        const positionPointer = i * 9;
        if (!breakArray[breakPointer + 1] && !breakArray[breakPointer + 2]) {
            editedLabels.push(labelArray[breakPointer]);
            editedLabels.push(labelArray[breakPointer + 1]);
            editedLabels.push(labelArray[breakPointer + 2]);
            editedBreaks.push(breakArray[breakPointer]);
            editedBreaks.push(breakArray[breakPointer + 1]);
            editedBreaks.push(breakArray[breakPointer + 2]);
            editedPositions.push(positionArray[positionPointer]);
            editedPositions.push(positionArray[positionPointer + 1]);
            editedPositions.push(positionArray[positionPointer + 2]);
            editedPositions.push(positionArray[positionPointer + 3]);
            editedPositions.push(positionArray[positionPointer + 4]);
            editedPositions.push(positionArray[positionPointer + 5]);
            editedPositions.push(positionArray[positionPointer + 6]);
            editedPositions.push(positionArray[positionPointer + 7]);
            editedPositions.push(positionArray[positionPointer + 8]);
            editedColors.push(colorArray[breakPointer]);
            editedColors.push(colorArray[breakPointer + 1]);
            editedColors.push(colorArray[breakPointer + 2]);
            editedPMs.push(convertedRibbonObject.pointmasterArray[breakPointer]);
            editedPMs.push(convertedRibbonObject.pointmasterArray[breakPointer + 1]);
            editedPMs.push(convertedRibbonObject.pointmasterArray[breakPointer + 2]);
        } else {
            // console.log('X triangle break found')
            // console.log('skipping: '+positionArray[positionPointer]+','+positionArray[positionPointer+1]+','+positionArray[positionPointer+2]+','
            //                        +positionArray[positionPointer+3]+','+positionArray[positionPointer+4]+','+positionArray[positionPointer+5]+','
            //                        +positionArray[positionPointer+6]+','+positionArray[positionPointer+7]+','+positionArray[positionPointer+8])
        }
    }
    return {
        group: convertedRibbonObject.group,
        subgroup: convertedRibbonObject.subgroup,
        name: convertedRibbonObject.name,
        masterArray: convertedRibbonObject.masterArray,
        pointmasterArray: editedPMs,
        nobutton: convertedRibbonObject.nobutton,
        labelArray: editedLabels,
        positionArray: editedPositions,
        breakArray: editedBreaks,
        colorArray: editedColors,
        pairTriangleNormals: convertedRibbonObject.pairTriangleNormals
    };
}

class KinParser {
    // @brief Property that is filled in by the constructor as it parses the file. Read by the caller.
    kinemage: Kinemage;

    // @brief Constructor for the KinParser class.
    // @param data The string data to be parsed, including all lines in the file.
    constructor(data: string) {
        this._parse(data);
    }

    private _parse(data: string) {
        // http://kinemage.biochem.duke.edu/software/king.php

        const kinemage: Kinemage = {
            comments: [],
            kinemage: undefined,
            onewidth: undefined,
            viewDict: {},
            pdbfile: undefined,
            texts: [],
            text: '',
            captions: [],
            caption: '',
            groupDict: {},
            subgroupDict: {},
            masterDict: {},
            pointmasterDict: {},
            dotLists: [],
            vectorLists: [],
            ballLists: [],
            ribbonLists: [],
            groupsAnimate: [],
            activeAnimateGroup: -1,
            groupsAnimate2: [],
            activeAnimateGroup2: -1
        };
        this.kinemage = kinemage;

        // Keep a local copy of the ColorDict that we can update with new colors defined in the file.
        const localColorDict: { [k: string]: Color } = Object.assign({}, ColorDict);

        let currentGroup: string = '';
        let currentGroupMasters: string[];
        let currentSubgroup: string = '';
        let currentSubgroupMasters: string[];

        let isDotList = false;
        let prevDotLabel = '';
        let dotDefaultColor: Color;
        let dotLabel: string[], dotPosition: number[], dotColor: Color[], dotPointMasters: string[][];

        let isVectorList = false;
        let prevVecLabel = '';
        let prevVecPosition: number[] | null = null;
        let prevVecColor: Color | null = null;
        let vecDefaultColor: Color, vecDefaultWidth: number;
        let vecLabel1: string[], vecLabel2: string[], vecPosition1: number[], vecPosition2: number[], vecColor1: Color[], vecColor2: Color[];
        let vecWidth: number[], vecPointMasters: string[][];

        let isBallList = false;
        let prevBallLabel = '';
        let ballRadius: number[], ballDefaultColor: Color, ballDefaultRadius: number;
        let ballLabel: string[], ballPosition: number[], ballColor: Color[], ballPointMasters: string[][];

        let isRibbonList = false;
        let ribbonIsTriangles = false;
        let prevRibbonPointLabel = '';

        let ribbonListDefaultColor: Color = localColorDict['white'];
        let ribbonPointLabelArray: string[], ribbonPointPositionArray: number[], ribbonPointBreakArray: boolean[], ribbonPointColorArray: Color[];
        let ribbonPointMasters: string[][];

        let isText = false;
        let isCaption = false;

        let foundAnimate = false;
        let found2Animate = false;

        function _parseChunkOfLines(_i: number, _n: number, lines: string[]) {
            for (let i = _i; i < _n; ++i) {
                const line = lines[i];

                if (line[0] === '@') {
                    isDotList = false;
                    isVectorList = false;
                    isBallList = false;
                    isRibbonList = false;
                    isText = false;
                    isCaption = false;
                }

                if (!line) {
                    isDotList = false;
                    isVectorList = false;
                    isBallList = false;
                    isRibbonList = false;
                } else if (line.startsWith('@dot') /* dot or dotlist */) {
                    // @dotlist {x} color=white master={vdw contact} master={dots}

                    let { listColor, listName, listMasters, nobutton } = parseListDef(line, localColorDict);

                    isDotList = true;
                    prevDotLabel = '';
                    dotLabel = [];
                    dotPosition = [];
                    dotColor = [];
                    dotPointMasters = [];
                    dotDefaultColor = listColor;

                    if (currentGroupMasters) {
                        listMasters = listMasters.concat(currentGroupMasters);
                    }
                    if (currentSubgroupMasters) {
                        listMasters = listMasters.concat(currentSubgroupMasters);
                    }

                    kinemage.dotLists.push({
                        group: currentGroup,
                        subgroup: currentSubgroup,
                        name: listName,
                        masterArray: listMasters,
                        pointmasterArray: dotPointMasters,
                        nobutton: nobutton,
                        labelArray: dotLabel,
                        positionArray: dotPosition,
                        colorArray: dotColor
                    });
                } else if (line.startsWith('@vector') /* vector or vectorlist */) {
                    // @vectorlist {x} color=white master={small overlap} master={dots}

                    let { listMasters, listName, listWidth, listColor, nobutton } = parseListDef(line, localColorDict);

                    if (listMasters) {
                        listMasters.forEach(function (name: string) {
                            if (!kinemage.masterDict[name]) {
                                kinemage.masterDict[name] = {
                                    indent: false,
                                    visible: true
                                };
                            }
                        });
                    }

                    isVectorList = true;
                    prevVecLabel = '';
                    prevVecPosition = null;
                    prevVecColor = null;
                    vecLabel1 = [];
                    vecLabel2 = [];
                    vecPosition1 = [];
                    vecPosition2 = [];
                    vecColor1 = [];
                    vecColor2 = [];
                    vecWidth = [];
                    vecDefaultColor = listColor;
                    vecPointMasters = [];
                    vecDefaultWidth = 2;
                    if (listWidth) {
                        vecDefaultWidth = listWidth;
                    }

                    if (currentGroupMasters) {
                        listMasters = listMasters.concat(currentGroupMasters);
                    }
                    if (currentSubgroupMasters) {
                        listMasters = listMasters.concat(currentSubgroupMasters);
                    }

                    kinemage.vectorLists.push({
                        group: currentGroup,
                        subgroup: currentSubgroup,
                        name: listName,
                        masterArray: listMasters,
                        pointmasterArray: vecPointMasters,
                        nobutton: nobutton,
                        label1Array: vecLabel1,
                        label2Array: vecLabel2,
                        position1Array: vecPosition1,
                        position2Array: vecPosition2,
                        color1Array: vecColor1,
                        color2Array: vecColor2,
                        width: vecWidth
                    });
                } else if (line.startsWith('@ball') /* ball or balllist*/ || line.startsWith('@sphere') /* sphere or spherelist */) {
                    let { listName, listColor, listMasters, listRadius, nobutton } = parseListDef(line, localColorDict);

                    if (listMasters) {
                        listMasters.forEach(function (name: string) {
                            if (!kinemage.masterDict[name]) {
                                kinemage.masterDict[name] = {
                                    indent: false,
                                    visible: true
                                };
                            }
                        });
                    }

                    isBallList = true;

                    prevBallLabel = '';
                    ballLabel = [];
                    ballRadius = [];
                    ballPosition = [];
                    ballColor = [];
                    ballPointMasters = [];
                    ballDefaultColor = listColor;
                    ballDefaultRadius = listRadius !== undefined ? listRadius : 1;

                    if (currentGroupMasters) {
                        listMasters = listMasters.concat(currentGroupMasters);
                    }
                    if (currentSubgroupMasters) {
                        listMasters = listMasters.concat(currentSubgroupMasters);
                    }

                    kinemage.ballLists.push({
                        group: currentGroup,
                        subgroup: currentSubgroup,
                        name: listName,
                        masterArray: listMasters,
                        pointmasterArray: ballPointMasters,
                        nobutton: nobutton,
                        labelArray: ballLabel,
                        radiusArray: ballRadius,
                        positionArray: ballPosition,
                        colorArray: ballColor
                    });
                } else if (line.startsWith('@ribbon') /* ribbon or ribbonlist */ || line.startsWith('@triangle') /* triangle or trianglelist */) {
                    let { listMasters, listName, listColor, nobutton } = parseListDef(line, localColorDict);

                    if (listMasters) {
                        listMasters.forEach(function (name: string) {
                            if (!kinemage.masterDict[name]) {
                                kinemage.masterDict[name] = {
                                    indent: false,
                                    visible: true
                                };
                            }
                        });
                    }
                    isRibbonList = true;
                    ribbonIsTriangles = line.startsWith('@triangle'); /* triangle or trianglelist */
                    prevRibbonPointLabel = '';
                    ribbonPointLabelArray = [];
                    ribbonPointPositionArray = [];
                    ribbonPointBreakArray = [];
                    ribbonPointColorArray = [];
                    ribbonListDefaultColor = listColor;
                    ribbonPointMasters = [];

                    if (currentGroupMasters) {
                        listMasters = listMasters.concat(currentGroupMasters);
                    }
                    if (currentSubgroupMasters) {
                        listMasters = listMasters.concat(currentSubgroupMasters);
                    }

                    kinemage.ribbonLists.push({
                        group: currentGroup,
                        subgroup: currentSubgroup,
                        name: listName,
                        masterArray: listMasters,
                        pointmasterArray: ribbonPointMasters,
                        nobutton: nobutton,
                        labelArray: ribbonPointLabelArray,
                        positionArray: ribbonPointPositionArray,
                        breakArray: ribbonPointBreakArray,
                        colorArray: ribbonPointColorArray,
                        pairTriangleNormals: !ribbonIsTriangles
                    });
                } else if (line.startsWith('@text')) {
                    isText = true;
                    kinemage.texts.push(line.substr(5));
                } else if (line.startsWith('@caption')) {
                    isCaption = true;
                    kinemage.captions.push(line.substr(8));
                } else if (isDotList) {
                    // { CB  THR   1  A}sky  'P' 18.915,14.199,5.024

                    let { label, color, position, pointMasters } = parseListElm(line, localColorDict);

                    if (label === '"') {
                        label = prevDotLabel;
                    } else {
                        prevDotLabel = label;
                    }

                    if (color === undefined) {
                        color = dotDefaultColor;
                    }

                    dotLabel.push(label);
                    dotPosition.push(...position);
                    dotColor.push(color);
                    dotPointMasters.push(pointMasters);
                } else if (isVectorList) {
                    // { n   thr A   1  B13.79 1crnFH} P 17.047, 14.099, 3.625 { n   thr A   1  B13.79 1crnFH} L 17.047, 14.099, 3.625

                    const doubleLine = line.replace(/(?!^){/g, '\n{');
                    const splitLine = doubleLine.split(/\n/);

                    for (let i2 = 0; i2 < splitLine.length; i2++) {
                        const singlePointLine = splitLine[i2];
                        let { label, color, width, position, isLineBreak, pointMasters } = parseListElm(singlePointLine, localColorDict);

                        if (label === '"') {
                            label = prevVecLabel;
                        }

                        if (color === undefined) {
                            color = vecDefaultColor;
                        }

                        if (!isLineBreak) {
                            if (prevVecPosition !== null) {
                                if (width === undefined) {
                                    width = vecDefaultWidth;
                                }

                                vecLabel1.push(prevVecLabel);
                                vecPosition1.push(...prevVecPosition);
                                vecColor1.push(prevVecColor ? prevVecColor : vecDefaultColor);

                                vecLabel2.push(label);
                                vecPosition2.push(...position);
                                vecColor2.push(color);
                                vecWidth.push(width);

                                vecPointMasters.push(pointMasters);
                            }
                        }

                        prevVecLabel = label;
                        prevVecPosition = position;
                        prevVecColor = color;
                    }
                } else if (isBallList) {
                    // {cb arg A   1   1.431 -106.80} r=1.431  39.085, 8.083, 22.182

                    let { label, radius, color, position, pointMasters } = parseListElm(line, localColorDict);

                    if (label === '"') {
                        label = prevBallLabel;
                    } else {
                        prevBallLabel = label;
                    }

                    if (radius === undefined) {
                        radius = ballDefaultRadius;
                    }

                    if (color === undefined) {
                        color = ballDefaultColor;
                    }

                    ballLabel.push(label);
                    ballRadius.push(radius);
                    ballPosition.push(...position);
                    ballColor.push(color);
                    ballPointMasters.push(pointMasters);
                } else if (isRibbonList) {
                    let { label, color, position, isTriangleBreak, pointMasters } = parseListElm(line, localColorDict);

                    if (label === '"') {
                        label = prevRibbonPointLabel;
                    } else {
                        prevRibbonPointLabel = label;
                    }

                    if (color === undefined) {
                        color = ribbonListDefaultColor;
                    }

                    ribbonPointLabelArray.push(label);
                    ribbonPointPositionArray.push(...position);
                    ribbonPointBreakArray.push(isTriangleBreak);
                    ribbonPointColorArray.push(color);
                    ribbonPointMasters.push(pointMasters);
                } else if (isText) {
                    kinemage.texts.push(line);
                } else if (isCaption) {
                    kinemage.captions.push(line);
                } else if (line.startsWith('@kinemage')) {
                    kinemage.kinemage = parseInt(line.substr(9).trim());
                } else if (line.startsWith('@onewidth')) {
                    kinemage.onewidth = true;
                } else if (line.startsWith('@pdbfile')) {
                    kinemage.pdbfile = parseStr(line);
                } else if (line.startsWith('@group')) {
                    const { groupName, groupFlags, groupMasters } = parseGroup(line);
                    if (!kinemage.groupDict[groupName as string]) {
                        kinemage.groupDict[groupName as string] = {
                            dominant: false,
                            // If the groupFlags include animate or 2animate, set those to true in the groupDict. Otherwise, set them to false.
                            animate: groupFlags['animate'] ? true : false,
                            '2animate': groupFlags['2animate'] ? true : false,
                            // If the foundAnimate or found2Animate flags are true, set off to true; otherwise set it to the flags value.
                            off: (foundAnimate || found2Animate) ? true : groupFlags['off'] ? true : false
                        };
                        // If the animate or 2animate flags are found in the groupFlags, set foundAnimate
                        // or found2Animate to true, respectively. Also update the list and index.
                        if (groupFlags['animate']) {
                            foundAnimate = true;
                            kinemage.groupsAnimate.push(groupName as string);
                            kinemage.activeAnimateGroup = 0;
                        }
                        if (groupFlags['2animate']) {
                            found2Animate = true;
                            kinemage.groupsAnimate2.push(groupName as string);
                            kinemage.activeAnimateGroup2 = 0;
                        }
                        currentGroupMasters = groupMasters;
                    }
                    currentGroup = groupName;

                    if (currentGroupMasters) {
                        currentGroupMasters.forEach(function (master) {
                            if (!kinemage.masterDict[master]) {
                                kinemage.masterDict[master] = {
                                    indent: false,
                                    visible: true
                                };
                            }
                        });
                    }

                    for (const key in groupFlags as { [k: string]: boolean }) {
                        kinemage.groupDict[groupName as string][key] = (groupFlags as { [k: string]: boolean })[key];
                    }
                } else if (line.startsWith('@subgroup')) {
                    const { groupName, groupFlags, groupMasters } = parseGroup(line);

                    const combinedName = currentGroup + ':' + groupName as string;
                    if (!kinemage.subgroupDict[combinedName]) {
                        kinemage.subgroupDict[combinedName] = {
                            dominant: false,
                            // If the groupFlag includes "off", set off to true; otherwise, set it to false.
                            off: groupFlags['off'] ? true : false,
                            group: currentGroup
                        };
                        currentSubgroupMasters = groupMasters;
                    }
                    currentSubgroup = combinedName;

                    if (currentSubgroupMasters) {
                        currentSubgroupMasters.forEach(function (master) {
                            if (!kinemage.masterDict[master]) {
                                kinemage.masterDict[master] = {
                                    indent: false,
                                    visible: true
                                };
                            }
                        });
                    }

                    for (const key in groupFlags as { [k: string]: boolean }) {
                        kinemage.subgroupDict[combinedName as string][key] = (groupFlags as { [k: string]: boolean })[key];
                    }
                } else if (line.startsWith('@master')) {
                    const name = parseStr(line);
                    const flag = parseFlag(line);

                    if (!kinemage.masterDict[name]) {
                        kinemage.masterDict[name] = {
                            indent: false,
                            visible: true
                        };
                    }

                    // TODO: There can be more than one flag on a @master line: indent, off, nobutton
                    if (flag === 'on') {
                        kinemage.masterDict[name].visible = true;
                    } else if (flag === 'off') {
                        kinemage.masterDict[name].visible = false;
                    } else if (flag === 'indent') {
                        kinemage.masterDict[name].indent = true;
                    } else if (!flag) {
                        // nothing to do
                    }
                } else if (line.startsWith('@pointmaster')) {
                    const { name, tags, on } = parsePointmaster(line);
                    if (name.length > 0 && tags.length > 0) {

                        // Ensure that we have a masterDict entry for this pointmaster name, even though it doesn't have any flags of its own.
                        if (!kinemage.masterDict[name]) {
                            kinemage.masterDict[name] = {
                                indent: false,
                                visible: on !== false // If the on variable is explicitly false, set visible to false. Otherwise, set it to true.
                            };
                        }

                        // Add the mapping to point each single-character tag to the pointmaster name in the pointmasterDict.
                        for (let i = 0; i < tags.length; i++) {
                            kinemage.pointmasterDict[tags[i]] = name;
                        }
                    }
                } else if (line.startsWith('@colorset')) {
                    // We have a string inside curly brackets {} followed by the name of an existing dictionary color.
                    const colorName = parseStr(line);
                    const colorReference = parseFlag(line);
                    if (colorReference && colorReference in localColorDict) {
                        localColorDict[colorName] = localColorDict[colorReference];
                    }
                } else if (/^@(\d*)viewid\b/.test(line)) {
                    const m = line.match(/^@(\d*)viewid\b/);
                    const viewCount = (m && m[1] && m[1].length > 0) ? parseInt(m[1], 10) : 1;
                    if (!kinemage.viewDict[viewCount]) kinemage.viewDict[viewCount] = {};
                    kinemage.viewDict[viewCount].name = parseStr(line);
                } else if (/^@(\d*)center\b/.test(line)) {
                    // Match all of the line after center as another string.
                    const m = line.match(/^@(\d*)center\b\s*(.*)$/);
                    const viewCount = (m && m[1] && m[1].length > 0) ? parseInt(m[1], 10) : 1;
                    // Pull out the three whitespace-separated numbers after the keyword.  Parse each as a float and
                    // add them to a length-3 list of numbers.
                    const rest = (m && m[2]) ? m[2].trim() : '';
                    // Split on whitespace and take the first three tokens, parsed as floating-point numbers, as the center coordinates.
                    const parts = rest.length > 0 ? rest.split(/\s+/).filter(Boolean) : [];
                    const centerTokens = parts.slice(0, 3).map(parseFloat);
                    // If the length is 3 and all are valid numbers, add the list of three numbers to the view dictionary.
                    if (centerTokens.length === 3 && centerTokens.every(num => !isNaN(num))) {
                        if (!kinemage.viewDict[viewCount]) kinemage.viewDict[viewCount] = {};
                        kinemage.viewDict[viewCount].center = centerTokens;
                    }
                } else if (/^@(\d*)matrix\b/.test(line)) {
                    // Match all of the line after matrix as another string.
                    const m = line.match(/^@(\d*)matrix\b\s*(.*)$/);
                    const viewCount = (m && m[1] && m[1].length > 0) ? parseInt(m[1], 10) : 1;
                    // Pull out the nine whitespace-separated numbers after the keyword.  Parse each as a float and
                    // add them to a length-9 list of numbers.
                    const rest = (m && m[2]) ? m[2].trim() : '';
                    // Split on whitespace and take the first nine tokens, parsed as floating-point numbers, as the matrix values.
                    const parts = rest.length > 0 ? rest.split(/\s+/).filter(Boolean) : [];
                    const matrixTokens = parts.slice(0, 9).map(parseFloat);
                    // If the length is 9 and all are valid numbers, add the list of nine numbers to the view dictionary.
                    if (matrixTokens.length === 9 && matrixTokens.every(num => !isNaN(num))) {
                        if (!kinemage.viewDict[viewCount]) kinemage.viewDict[viewCount] = {};
                        kinemage.viewDict[viewCount].matrix = matrixTokens;
                    }
                } else if (/^@(\d*)span\b/.test(line)) {
                    // Match all of the line after span as another string.
                    const m = line.match(/^@(\d*)span\b\s*(.*)$/);
                    const viewCount = (m && m[1] && m[1].length > 0) ? parseInt(m[1], 10) : 1;
                    // Pull out the remainder of the line and parse it as a float.
                    const rest = (m && m[2]) ? m[2].trim() : '';
                    const spanValue = parseFloat(rest);
                    // If it is a valid number, add it to the view dictionary.
                    if (!isNaN(spanValue)) {
                        if (!kinemage.viewDict[viewCount]) kinemage.viewDict[viewCount] = {};
                        kinemage.viewDict[viewCount].span = spanValue;
                    }
                } else if (/^@(\d*)zoom\b/.test(line)) {
                    // Match all of the line after zoom as another string.
                    const m = line.match(/^@(\d*)zoom\b\s*(.*)$/);
                    const viewCount = (m && m[1] && m[1].length > 0) ? parseInt(m[1], 10) : 1;
                    // Pull out the remainder of the line and parse it as a float.
                    const rest = (m && m[2]) ? m[2].trim() : '';
                    const zoomValue = parseFloat(rest);
                    // If it is a valid number, add it to the view dictionary.
                    if (!isNaN(zoomValue)) {
                        if (!kinemage.viewDict[viewCount]) kinemage.viewDict[viewCount] = {};
                        kinemage.viewDict[viewCount].zoom = zoomValue;
                    }
                } else if (/^@(\d*)zslab\b/.test(line)) {
                    // Match all of the line after zslab as another string.
                    const m = line.match(/^@(\d*)zslab\b\s*(.*)$/);
                    const viewCount = (m && m[1] && m[1].length > 0) ? parseInt(m[1], 10) : 1;
                    // Pull out the remainder of the line and parse it as a float.
                    const rest = (m && m[2]) ? m[2].trim() : '';
                    const zslabValue = parseFloat(rest);
                    // If it is a valid number, add it to the view dictionary.
                    if (!isNaN(zslabValue)) {
                        if (!kinemage.viewDict[viewCount]) kinemage.viewDict[viewCount] = {};
                        kinemage.viewDict[viewCount].zslab = zslabValue;
                    }
                } else {
                    console.log('Kinemage: Unrecognized line: ' + line);
                }
            }
        }

        // Break the file into a list of lines and then parse them all.
        const lines = data.split(/\r?\n/);
        _parseChunkOfLines(0, lines.length, lines);

        kinemage.text = kinemage.texts.join('\n').trim();
        kinemage.caption = kinemage.captions.join('\n').trim();
        if (kinemage.ribbonLists) {
            const convertedLists: RibbonObject[] = [];
            kinemage.ribbonLists.forEach(function (listObject) {
                convertedLists.push(removePointBreaksTriangleArrays(convertKinTriangleArrays(listObject)));
            });
            kinemage.ribbonLists = convertedLists;
        }

    }
}

export { KinParser };
