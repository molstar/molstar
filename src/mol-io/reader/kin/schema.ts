/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

export interface Kinemage {
  readonly comments: ReadonlyArray<string>
  kinemage?: number,
  onewidth?: any,
  viewDict: { [id: number]: View },
  pdbfile?: string,
  text: string,
  texts: string[],
  captions: string[],
  caption: string,
  groupDict: { [k: string]: { [k: string]: boolean } },
  subgroupDict: { [k: string]: any },   ///< Subgroup key is "GroupName:SubgroupName" to preserve tree structure
  masterDict: { [k: string]: { indent: boolean, visible: boolean } },
  pointmasterDict: { [k: string]: any },
  dotLists: DotList[],
  vectorLists: VectorList[],
  ballLists: BallList[],
  ribbonLists: RibbonObject[]
}

/** Common base for all list-like objects in a kinemage */
export interface KinListBase {
  group: string,                  ///< Name of the group this List belongs to
  subgroup: string,               ///< Name of the subgroup this List belongs to
  name?: string,                  ///< Optional name of the whole List
  masterArray: any[]              ///< Array of master names per List, not per element
}

export interface DotList extends KinListBase {
  labelArray: string[],           ///< Array of labels per element
  positionArray: number[],        ///< Catenation of x, y, z for each element, 3x as many as elements
  colorArray: number[]            ///< Catenation of r, g, b for each element, 3x as many as elements
}

export interface BallList extends KinListBase {
  labelArray: string[],           ///< Array of labels per element
  positionArray: number[],        ///< Catenation of x, y, z for each element, 3x as many as elements
  colorArray: number[],           ///< Catenation of r, g, b for each element, 3x as many as elements
  radiusArray: number[]           ///< A single radius per element
}

export interface RibbonObject extends KinListBase {
  labelArray: string[],           ///< Array of labels per element
  positionArray: number[],        ///< Catenation of x, y, z for each element, 9x as many as triangles (3 vertices per triangle)
  colorArray: number[],           ///< Catenation of r, g, b for each element, 9x as many as triangles (3 colors per triangle)
  breakArray: boolean[],          ///< A single boolean per element indicating if there is a break there
  pairTriangleNormals: boolean    ///< Whether to pair every other triangle normal for lighting (true for ribbons, false for triangles)
}

export interface VectorList extends KinListBase {
  label1Array: string[],          ///< Array of labels per element
  label2Array: string[],          ///< Array of labels per element
  position1Array: number[],       ///< Catenation of x, y, z for each element, 3x as many as elements
  position2Array: number[],       ///< Catenation of x, y, z for each element, 3x as many as elements
  color1Array: number[],          ///< Catenation of r, g, b for each element, 3x as many as elements
  color2Array: number[],          ///< Catenation of r, g, b for each element, 3x as many as elements
  width: number[]                 ///< A single width per element
}

export interface View {
  name?: string,                  ///< Optional name of the View
  center?: number[],              ///< X, Y, Z of the center of the view; the model rotates arond this point
  matrix?: number[],              ///< Specifies and orthonormal rotation matrix defining view orientation
  span?: number,                  ///< Specifies the (smaller of) width or height of the view in world coordinates at the center
  zoom?: number,                  ///< Alternate zoom specification, indicates how much of the model is visible, 1=all, 2=half
  zslab?: number                  ///< Distance from the center to the near and far clipping planes, 200 means same as span (half is percent of half span)
}
