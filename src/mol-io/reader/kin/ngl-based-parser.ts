/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * Based on earlier kin-parser.ts file from the NGL project (see second author notice below).
 * @file Ported NGL-based Kinemage file parser
 * @author ReliaSolve <russ@reliasolve.com>
 * @private
 */

/**
 * file Kin Parser
 * author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Debug, Log, ParserRegistry } from '../globals'
// import { Vector3 } from 'three'
// import Parser from './parser'

/// @todo Fill in commments

import { Kinemage, RibbonObject } from './schema';

function hsvToRgb (h: number, s: number, v: number) {
  h /= 360
  s /= 100
  v /= 100
  let r, g, b
  const i = Math.floor(h * 6)
  const f = h * 6 - i
  const p = v * (1 - s)
  const q = v * (1 - f * s)
  const t = v * (1 - (1 - f) * s)
  switch (i % 6) {
    case 0: r = v; g = t; b = p; break
    case 1: r = q; g = v; b = p; break
    case 2: r = p; g = v; b = t; break
    case 3: r = p; g = q; b = v; break
    case 4: r = t; g = p; b = v; break
    case 5: r = v; g = p; b = q; break
  }
  return [ r, g, b ] as number []
}

const ColorDict: {[k: string]: number[]} = {
  red: hsvToRgb(0, 100, 100),
  orange: hsvToRgb(20, 100, 100),
  gold: hsvToRgb(40, 100, 100),
  yellow: hsvToRgb(60, 100, 100),
  lime: hsvToRgb(80, 100, 100),
  green: hsvToRgb(120, 80, 100),
  sea: hsvToRgb(150, 100, 100),
  cyan: hsvToRgb(180, 100, 85),
  sky: hsvToRgb(210, 75, 95),
  blue: hsvToRgb(240, 70, 100),
  purple: hsvToRgb(275, 75, 100),
  magenta: hsvToRgb(300, 95, 100),
  hotpink: hsvToRgb(335, 100, 100),
  pink: hsvToRgb(350, 55, 100),
  peach: hsvToRgb(25, 75, 100),
  lilac: hsvToRgb(275, 55, 100),
  pinktint: hsvToRgb(340, 30, 100),
  peachtint: hsvToRgb(25, 50, 100),
  yellowtint: hsvToRgb(60, 50, 100),
  greentint: hsvToRgb(135, 40, 100),
  bluetint: hsvToRgb(220, 40, 100),
  lilactint: hsvToRgb(275, 35, 100),
  white: hsvToRgb(0, 0, 100),
  gray: hsvToRgb(0, 0, 50),
  brown: hsvToRgb(20, 45, 75),
  deadwhite: [ 1, 1, 1 ],
  deadblack: [ 0, 0, 0 ],
  invisible: [ 0, 0, 0 ]
}

const reWhitespaceComma = /[\s,]+/
const reCurlyWhitespace = /[^{}\s]*{[^{}]+}|[^{}\s]+/g
const reTrimCurly = /^{+|}+$/g
const reTrimQuotes = /^['"]+|['"]+$/g
const reCollapseEqual = /\s*=\s*/g

function parseListDef (line: string) {
  let name
  let defaultColor
  let master = []
  let width

  line = line.replace(reCollapseEqual, '=')

  const lm = line.match(reCurlyWhitespace) as RegExpMatchArray
  for (let j = 1; j < lm.length; ++j) {
    const e = lm[ j ]
    if (e[ 0 ] === '{') {
      name = e.substring(1, e.length - 1)
    } else {
      const es = e.split('=')
      if (es.length === 2) {
        if (es[ 0 ] === 'color') {
          defaultColor = ColorDict[ es[ 1 ] ]
        } else if (es[ 0 ] === 'width') {
          width = parseInt(es[ 1 ])
        } else if (es[ 0 ] === 'master') {
          master.push(es[ 1 ].replace(reTrimCurly, ''))
        }
      }
    }
  }

  return {
    listName: name,
    listColor: defaultColor,
    listMasters: master,
    listWidth: width
  }
}

function parseListElm (line: string) {
  line = line.trim()

  const idx1 = line.indexOf('{')
  const idx2 = line.indexOf('}')
  const ls = line.substr(idx2 + 1).split(reWhitespaceComma)

  const label = line.substr(idx1 + 1, idx2 - 1)
  const position = [
    parseFloat(ls[ ls.length - 3 ]),
    parseFloat(ls[ ls.length - 2 ]),
    parseFloat(ls[ ls.length - 1 ])
  ]
  let color, width, radius
  let lineBreak = false
  let triangleBreak = false
  for (let lsindex = 4; lsindex <= ls.length; lsindex++) {
    const literal = ls[ ls.length - lsindex ]
    if (literal in ColorDict) {
      color = ColorDict[ ls[ ls.length - lsindex ] ]
    }
    if (literal.startsWith('width')) {
      width = parseInt(literal.substring(5))
    }
    if (literal.startsWith('r=')) {
      radius = parseFloat(literal.split('=')[1])
    }
    if (literal.startsWith('P')) {
      lineBreak = true
    }
    if (literal.startsWith('X')) {
      triangleBreak = true
    }
  }
  // const color = line[ idx2 + 1 ] === ' ' ? undefined : ColorDict[ ls[ 0 ] ]

  return {
    label: label,
    position: position,
    color: color,
    radius: radius,
    width: width,
    isLineBreak: lineBreak,
    isTriangleBreak: triangleBreak
  }
}

function parseStr (line: string) {
  const start = line.indexOf('{')
  const end = line.indexOf('}')
  return line.substring(
    start !== -1 ? start + 1 : 0,
    end !== -1 ? end : undefined
  ).trim()
}

function parseFlag (line: string) {
  const end = line.indexOf('}')
  return end === -1 ? undefined : line.substr(end + 1).trim()
}

function parseGroup (line: string) {
  let name:string = ''
  let master:string[] = []
  let flags: {[k: string]: string|boolean} = {}

  line = line.replace(reCollapseEqual, '=')

  const lm = line.match(reCurlyWhitespace)  as RegExpMatchArray
  for (let j = 1; j < lm.length; ++j) {
    const e = lm[ j ]
    if (e[ 0 ] === '{') {
      name = e.substring(1, e.length - 1)
    } else {
      const es = e.split('=')
      if (es.length === 2) {
        if (es[ 0 ] === 'master') {
          master.push(es[ 1 ].replace(reTrimCurly, ''))
        } else {
          flags[ es[ 0 ] ] = es[ 1 ].replace(reTrimCurly, '')
        }
      } else {
        flags[ es[ 0 ] ] = true
      }
    }
  }

  return { groupName: name,
           groupFlags: flags,
           groupMasters: master
  }
}
function convertKinTriangleArrays (ribbonObject: RibbonObject) {
  // have to convert ribbons/triangle lists from stripdrawmode to normal drawmode
  // index                    [ 0 1 2 3 4 5 6 7 8 91011 ]
  // label [ 0 1 2 3 4 5 ] to [ 0 1 2 1 2 3 2 3 4 3 4 5 ]
  // convertedindex                                      [ 0 1 2 3 4 5 6 7 8 91011121314151617181920212223242526 ]
  // index          [ 0 1 2 3 4 5 6 7 8 91011121314 ]    [ 0 1 2 3 4 5 6 7 8 3 4 5 6 7 8 91011 6 7 8 91011121314 ]
  // position/color [ 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 ] to [ 0 0 0 1 1 1 2 2 2 1 1 1 2 2 2 3 3 3 2 2 2 3 3 3 4 4 4 ]
  let { labelArray, positionArray, colorArray, breakArray } = ribbonObject
  let convertedLabels = []
  for (let i = 0; i < (labelArray.length - 2) * 3; ++i) {
    convertedLabels[i] = labelArray[i - Math.floor(i / 3) * 2]
  }
  let convertedBreaks = []
  for (let i = 0; i < (breakArray.length - 2) * 3; ++i) {
    convertedBreaks[i] = breakArray[i - Math.floor(i / 3) * 2]
  }
  let convertedPositions = []
  for (let i = 0; i < (positionArray.length / 3 - 2) * 9; ++i) {
    convertedPositions[i] = positionArray[i - Math.floor(i / 9) * 6]
  }
  let convertedColors = []
  for (let i = 0; i < (colorArray.length / 3 - 2) * 9; ++i) {
    convertedColors[i] = colorArray[i - Math.floor(i / 9) * 6]
  }
  let vector3Positions = []
  for (let i = 0; i < (convertedPositions.length) / 3; ++i) {
    vector3Positions.push([convertedPositions[i * 3], convertedPositions[i * 3] + 1, convertedPositions[i * 3] + 2])
  }
  //let normals = []
  //for (let i = 0; i < vector3Positions.length - 1; ++i) {
  //  let normalVec3 = vector3Positions[i].cross(vector3Positions[i + 1])
  //  normals.push(normalVec3.x)
  //  normals.push(normalVec3.y)
  //  normals.push(normalVec3.z)
  //}
  return {
    name: ribbonObject.name,
    masterArray: ribbonObject.masterArray,
    labelArray: convertedLabels,
    positionArray: convertedPositions,
    breakArray: convertedBreaks,
    colorArray: convertedColors
  }
}

function removePointBreaksTriangleArrays (convertedRibbonObject: RibbonObject) {
  // after converting ribbon/triangle arrys to drawmode, removed point break triangles
  // label [ 0 1 2 3 4 5 ] to [ 0 1 2 1 2 3 2 3 4 3 4 5 ]
  // position/color [ 0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 ] to [ 0 0 0 1 1 1 2 2 2 1 1 1 2 2 2 3 3 3 2 2 2 3 3 3 4 4 4 ]
  let { labelArray, positionArray, colorArray, breakArray } = convertedRibbonObject
  let editedLabels = []
  let editedPositions = []
  let editedColors = []
  let editedBreaks = []
  for (let i = 0; i < breakArray.length / 3; i++) {
    let breakPointer = i * 3
    let positionPointer = i * 9
    if (!breakArray[breakPointer+1]&&!breakArray[breakPointer+2]) {
      editedLabels.push(labelArray[breakPointer])
      editedLabels.push(labelArray[breakPointer+1])
      editedLabels.push(labelArray[breakPointer+2])
      editedBreaks.push(breakArray[breakPointer])
      editedBreaks.push(breakArray[breakPointer+1])
      editedBreaks.push(breakArray[breakPointer+2])
      editedPositions.push(positionArray[positionPointer])
      editedPositions.push(positionArray[positionPointer+1])
      editedPositions.push(positionArray[positionPointer+2])
      editedPositions.push(positionArray[positionPointer+3])
      editedPositions.push(positionArray[positionPointer+4])
      editedPositions.push(positionArray[positionPointer+5])
      editedPositions.push(positionArray[positionPointer+6])
      editedPositions.push(positionArray[positionPointer+7])
      editedPositions.push(positionArray[positionPointer+8])
      editedColors.push(colorArray[positionPointer])
      editedColors.push(colorArray[positionPointer+1])
      editedColors.push(colorArray[positionPointer+2])
      editedColors.push(colorArray[positionPointer+3])
      editedColors.push(colorArray[positionPointer+4])
      editedColors.push(colorArray[positionPointer+5])
      editedColors.push(colorArray[positionPointer+6])
      editedColors.push(colorArray[positionPointer+7])
      editedColors.push(colorArray[positionPointer+8])
    } else {
      //console.log('X triangle break found')
      //console.log('skipping: '+positionArray[positionPointer]+','+positionArray[positionPointer+1]+','+positionArray[positionPointer+2]+','
      //                        +positionArray[positionPointer+3]+','+positionArray[positionPointer+4]+','+positionArray[positionPointer+5]+','
      //                        +positionArray[positionPointer+6]+','+positionArray[positionPointer+7]+','+positionArray[positionPointer+8])
    }
  }
  return {
    name: convertedRibbonObject.name,
    masterArray: convertedRibbonObject.masterArray,
    labelArray: editedLabels,
    positionArray: editedPositions,
    breakArray: editedBreaks,
    colorArray: editedColors
  }
}

class KinParser {
  /// @brief Property that is filled in by the constructor as it parses the file. Read by the caller.
  kinemage: Kinemage

  /// @brief Constructor for the KinParser class.
  /// @param data The string data to be parsed, including all lines in the file.
  constructor(data: string) {
    this._parse(data);
  }

  private _parse (data: string) {
    // http://kinemage.biochem.duke.edu/software/king.php

    const kinemage: Kinemage = {
      comments: [],
      kinemage: undefined,
      onewidth: undefined,
      '1viewid': undefined,
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
      ribbonLists: []
    }
    this.kinemage = kinemage

    let currentGroupMasters: string[]
    let currentSubgroupMasters: string[]

    let isDotList = false
    let prevDotLabel = ''
    let dotDefaultColor: number[]
    let dotLabel: string[], dotPosition: number[], dotColor: number[]

    let isVectorList = false
    let prevVecLabel = ''
    let prevVecPosition: number[]|null = null
    let prevVecColor: number[]|null = null
    let vecDefaultColor: number[], vecDefaultWidth: number[]
    let vecLabel1: string[], vecLabel2: string[], vecPosition1: number[], vecPosition2: number[], vecColor1: number[], vecColor2: number[]

    let isBallList = false
    let prevBallLabel = ''
    let ballRadius: number[], ballDefaultColor: number[]
    let ballLabel: string[], ballPosition: number[], ballColor: number[]

    let isRibbonList = false
    let prevRibbonPointLabel = ''

    let ribbonListDefaultColor: number[]
    let ribbonPointLabelArray: string[], ribbonPointPositionArray: number[], ribbonPointBreakArray: boolean[], ribbonPointColorArray: number[]

    let isText = false
    let isCaption = false

    // @vectorlist {mc} color= white  master= {mainchain}
    // { n   thr A   1  B13.79 1crnFH} P 17.047, 14.099, 3.625 { n   thr A   1  B13.79 1crnFH} L 17.047, 14.099, 3.625

    // @dotlist {x} color=white master={vdw contact} master={dots}
    // { CB  THR   1  A}sky  'P' 18.915,14.199,5.024

    function _parseChunkOfLines (_i: number, _n: number, lines: string[]) {
      for (let i = _i; i < _n; ++i) {
        const line = lines[ i ]

        if (line[ 0 ] === '@') {
          isDotList = false
          isVectorList = false
          isBallList = false
          isRibbonList = false
          isText = false
          isCaption = false
        }

        if (!line) {
          isDotList = false
          isVectorList = false
          isBallList = false
          isRibbonList = false
        } else if (line.startsWith('@dotlist')) {
          // @dotlist {x} color=white master={vdw contact} master={dots}

          let { listColor, listName, listMasters } = parseListDef(line)

          isDotList = true
          prevDotLabel = ''
          dotLabel = []
          dotPosition = []
          dotColor = []
          dotDefaultColor = listColor as number[]

          if (currentGroupMasters) {
            listMasters = listMasters.concat(currentGroupMasters)
          }
          if (currentSubgroupMasters) {
            listMasters = listMasters.concat(currentSubgroupMasters)
          }

          kinemage.dotLists.push({
            name: listName,
            masterArray: listMasters,
            labelArray: dotLabel,
            positionArray: dotPosition,
            colorArray: dotColor
          })
        } else if (line.startsWith('@vectorlist')) {
          // @vectorlist {x} color=white master={small overlap} master={dots}

          let { listMasters, listName, listWidth, listColor } = parseListDef(line)

          if (listMasters) {
            listMasters.forEach(function (name: string) {
              if (!kinemage.masterDict[ name ]) {
                kinemage.masterDict[ name ] = {
                  indent: false,
                  visible: false
                }
              }
            })
          }

          isVectorList = true
          prevVecLabel = ''
          prevVecPosition = null
          prevVecColor = null
          vecLabel1 = []
          vecLabel2 = []
          vecPosition1 = []
          vecPosition2 = []
          vecColor1 = []
          vecColor2 = []
          vecDefaultColor = listColor as number[]
          vecDefaultWidth = []
          if (listWidth) {
            vecDefaultWidth.push(listWidth)
          }

          if (currentGroupMasters) {
            listMasters = listMasters.concat(currentGroupMasters)
          }
          if (currentSubgroupMasters) {
            listMasters = listMasters.concat(currentSubgroupMasters)
          }

          kinemage.vectorLists.push({
            name: listName,
            masterArray: listMasters,
            label1Array: vecLabel1,
            label2Array: vecLabel2,
            position1Array: vecPosition1,
            position2Array: vecPosition2,
            color1Array: vecColor1,
            color2Array: vecColor2,
            width: vecDefaultWidth
          })
        } else if (line.startsWith('@balllist')) {
          let { listName, listColor, listMasters } = parseListDef(line)

          if (listMasters) {
            listMasters.forEach(function (name: string) {
              if (!kinemage.masterDict[ name ]) {
                kinemage.masterDict[ name ] = {
                  indent: false,
                  visible: false
                }
              }
            })
          }

          isBallList = true

          prevBallLabel = ''
          ballLabel = []
          ballRadius = []
          ballPosition = []
          ballColor = []
          ballDefaultColor = listColor as number[]

          if (currentGroupMasters) {
            listMasters = listMasters.concat(currentGroupMasters)
          }
          if (currentSubgroupMasters) {
            listMasters = listMasters.concat(currentSubgroupMasters)
          }

          kinemage.ballLists.push({
            name: listName,
            masterArray: listMasters,
            labelArray: ballLabel,
            radiusArray: ballRadius,
            positionArray: ballPosition,
            colorArray: ballColor
          })
        } else if (line.startsWith('@ribbonlist')||line.startsWith('@trianglelist')) {
          let { listMasters, listName, listColor } = parseListDef(line)

          if (listMasters) {
            listMasters.forEach(function (name: string) {
              if (!kinemage.masterDict[ name ]) {
                kinemage.masterDict[ name ] = {
                  indent: false,
                  visible: false
                }
              }
            })
          }
          isRibbonList = true
          prevRibbonPointLabel = ''
          ribbonPointLabelArray = []
          ribbonPointPositionArray = []
          ribbonPointBreakArray = []
          ribbonPointColorArray = []
          ribbonListDefaultColor = listColor as number[]

          if (currentGroupMasters) {
            listMasters = listMasters.concat(currentGroupMasters)
          }
          if (currentSubgroupMasters) {
            listMasters = listMasters.concat(currentSubgroupMasters)
          }

          kinemage.ribbonLists.push({
            name: listName,
            masterArray: listMasters,
            labelArray: ribbonPointLabelArray,
            positionArray: ribbonPointPositionArray,
            breakArray: ribbonPointBreakArray,
            colorArray: ribbonPointColorArray
          })
        } else if (line.startsWith('@text')) {
          isText = true
          kinemage.texts.push(line.substr(5))
        } else if (line.startsWith('@caption')) {
          isCaption = true
          kinemage.captions.push(line.substr(8))
        } else if (isDotList) {
          // { CB  THR   1  A}sky  'P' 18.915,14.199,5.024

          let { label, color, position } = parseListElm(line)

          if (label === '"') {
            label = prevDotLabel
          } else {
            prevDotLabel = label
          }

          if (color === undefined) {
            color = dotDefaultColor
          }

          dotLabel.push(label)
          dotPosition.push(...position)
          dotColor.push(...color)
        } else if (isVectorList) {
          // { n   thr A   1  B13.79 1crnFH} P 17.047, 14.099, 3.625 { n   thr A   1  B13.79 1crnFH} L 17.047, 14.099, 3.625

          let doubleLine = line.replace(/(?!^){/g, '\n{')
          let splitLine = doubleLine.split(/\n/)

          for (var i2 = 0; i2 < splitLine.length; i2++) {
            let singlePointLine = splitLine[i2]
            let { label, color, width, position, isLineBreak } = parseListElm(singlePointLine)

            if (label === '"') {
              label = prevVecLabel
            } else {
              prevVecLabel = label
            }

            if (color === undefined) {
              color = vecDefaultColor
            }

            if (!isLineBreak) {
              if (prevVecPosition !== null) {
                if (width) {
                  vecDefaultWidth.push(width)
                }

                vecLabel1.push(prevVecLabel)
                vecPosition1.push(...prevVecPosition)
                vecColor1.push(...prevVecColor as number[])

                vecLabel2.push(label)
                vecPosition2.push(...position)
                vecColor2.push(...color)

              }
            }

            prevVecLabel = label
            prevVecPosition = position
            prevVecColor = color
          }
        } else if (isBallList) {
          // {cb arg A   1   1.431 -106.80} r=1.431  39.085, 8.083, 22.182

          let { label, radius, color, position } = parseListElm(line)

          if (label === '"') {
            label = prevBallLabel
          } else {
            prevBallLabel = label
          }

          if (radius === undefined) {
            radius = 1 // temporary default radius
          }

          if (color === undefined) {
            color = ballDefaultColor
          }

          ballLabel.push(label)
          ballRadius.push(radius)
          ballPosition.push(...position)
          ballColor.push(...color)
        } else if (isRibbonList) {
          let { label, color, position, isTriangleBreak } = parseListElm(line)

          if (label === '"') {
            label = prevRibbonPointLabel
          } else {
            prevRibbonPointLabel = label
          }

          if (color === undefined) {
            color = ribbonListDefaultColor
          }

          ribbonPointLabelArray.push(label)
          ribbonPointPositionArray.push(...position)
          ribbonPointBreakArray.push(isTriangleBreak)
          ribbonPointColorArray.push(...color)
        } else if (isText) {
          kinemage.texts.push(line)
        } else if (isCaption) {
          kinemage.captions.push(line)
        } else if (line.startsWith('@kinemage')) {
          kinemage.kinemage = parseInt(line.substr(9).trim())
        } else if (line.startsWith('@onewidth')) {
          kinemage.onewidth = true
        } else if (line.startsWith('@1viewid')) {
          kinemage[ '1viewid' ] = parseStr(line)
        } else if (line.startsWith('@pdbfile')) {
          kinemage.pdbfile = parseStr(line)
        } else if (line.startsWith('@group')) {
          let { groupName, groupFlags, groupMasters } = parseGroup(line)
          if (!kinemage.groupDict[ groupName as string ]) {
            kinemage.groupDict[ groupName as string ] = {
              dominant: false,
              animate: false
            }
            currentGroupMasters = groupMasters
          }

          if (currentGroupMasters) {
            currentGroupMasters.forEach(function (master) {
              if (!kinemage.masterDict[ master ]) {
                kinemage.masterDict[ master ] = {
                  indent: false,
                  visible: false
                }
              }
            })
          }

          for (let key in groupFlags as {[k: string]: boolean}) {
            kinemage.groupDict[ groupName as string ][ key ] = (groupFlags as {[k: string]: boolean})[ key ]
          }
        } else if (line.startsWith('@subgroup')) {
          const { groupName, groupFlags, groupMasters } = parseGroup(line)

          if (!kinemage.subgroupDict[ groupName as string ]) {
            kinemage.subgroupDict[ groupName as string ] = {
              dominant: false,
              animate: false
            }
            currentSubgroupMasters = groupMasters
          }

          if (currentSubgroupMasters) {
            currentSubgroupMasters.forEach(function (master) {
              if (!kinemage.masterDict[ master ]) {
                kinemage.masterDict[ master ] = {
                  indent: false,
                  visible: false
                }
              }
            })
          }

          for (let key in groupFlags as {[k: string]: boolean}) {
            kinemage.subgroupDict[ groupName as string ][ key ] = (groupFlags as {[k: string]: boolean})[ key ]
          }
        } else if (line.startsWith('@master')) {
          const name = parseStr(line)
          const flag = parseFlag(line)

          if (!kinemage.masterDict[ name ]) {
            kinemage.masterDict[ name ] = {
              indent: false,
              visible: false
            }
          }

          if (flag === 'on') {
            kinemage.masterDict[ name ].visible = true
          } else if (flag === 'off') {
            kinemage.masterDict[ name ].visible = false
          } else if (flag === 'indent') {
            kinemage.masterDict[ name ].indent = true
          } else if (!flag) {
            // nothing to do
          }
        } else if (line.startsWith('@pointmaster')) {
          const { groupName, groupFlags } = parseGroup(line)

          kinemage.pointmasterDict[ groupName as string] = {
            id: Object.keys(groupFlags as {[k: string]: boolean})[ 0 ].replace(reTrimQuotes, '')
          }
        } else {
          console.log(line)
        }
      }
    }

    // Break the file into a list of lines and then parse them all.
    const lines = data.split(/\r?\n/)
    _parseChunkOfLines(0, lines.length, lines);

    kinemage.text = kinemage.texts.join('\n').trim()
    kinemage.caption = kinemage.captions.join('\n').trim()
    if (kinemage.ribbonLists) {
      let convertedLists: RibbonObject[] = []
      kinemage.ribbonLists.forEach(function (listObject) {
        convertedLists.push(removePointBreaksTriangleArrays(convertKinTriangleArrays(listObject)))
      })
      kinemage.ribbonLists = convertedLists
    }

  }
}

export default KinParser
