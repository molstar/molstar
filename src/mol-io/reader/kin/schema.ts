/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

export interface Kinemage {
  readonly comments: ReadonlyArray<string>
  kinemage?: number,
  onewidth?: any,
  '1viewid'?: string,
  pdbfile?: string,
  text: string,
  texts: string[],
  captions: string[],
  caption: string,
  groupDict: { [k: string]: { [k: string]: boolean } },
  subgroupDict: { [k: string]: any },
  masterDict: { [k: string]: { indent: boolean, visible: boolean } },
  pointmasterDict: { [k: string]: any },
  dotLists: DotList[],
  vectorLists: VectorList[],
  ballLists: BallList[],
  ribbonLists: RibbonObject[]
}

export interface DotList {
  name?: string,                  ///< Optional name of the whole List
  masterArray: any[],             ///< Array of master names per List, not per element
  labelArray: string[],           ///< Array of labels per List, not per element
  positionArray: number[],        ///< Catenation of x, y, z for each element, 3x as many as elements
  colorArray: number[]            ///< Catenation of r, g, b for each element, 3x as many as elements
}

export interface BallList {
  name?: string,                  ///< Optional name of the whole List
  masterArray: any[],             ///< Array of master names per List, not per element
  labelArray: string[],           ///< Array of labels per List, not per element
  positionArray: number[],        ///< Catenation of x, y, z for each element, 3x as many as elements
  colorArray: number[],           ///< Catenation of r, g, b for each element, 3x as many as elements
  radiusArray: number[]           ///< A single radius per element
}

export interface RibbonObject {
  name?: string,                  ///< Optional name of the whole Ribbon
  masterArray: any[],             ///< Array of master names per Ribbon, not per element
  labelArray: string[],           ///< Array of labels per Ribbon, not per element
  positionArray: number[],        ///< Catenation of x, y, z for each element, 3x as many as elements
  colorArray: number[],           ///< Catenation of r, g, b for each element, 3x as many as elements
  breakArray: boolean[]           ///< A single boolean per element indicating if there is a break there
}

export interface VectorList {
  name?: string,                  ///< Optional name of the whole List
  masterArray: any[],             ///< Array of master names per List, not per element
  label1Array: string[],          ///< Array of labels per List, not per element
  label2Array: string[],          ///< Array of labels per List, not per element
  position1Array: number[],       ///< Catenation of x, y, z for each element, 3x as many as elements
  position2Array: number[],       ///< Catenation of x, y, z for each element, 3x as many as elements
  color1Array: number[],          ///< Catenation of r, g, b for each element, 3x as many as elements
  color2Array: number[],          ///< Catenation of r, g, b for each element, 3x as many as elements
  width: number[]                 ///< A single width per element
}
