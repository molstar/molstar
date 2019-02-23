// /**
//  * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author Alexander Rose <alexander.rose@weirdbyte.de>
//  */

// import { ParamDefinition as PD } from 'mol-util/param-definition';
// import { Loci as _Loci, EmptyLoci } from 'mol-model/loci';
// import { Vec3 } from 'mol-math/linear-algebra';
// import { Mesh } from 'mol-geo/geometry/mesh/mesh';
// import { Representation, RepresentationContext } from 'mol-repr/representation';
// import { Subject } from 'rxjs';
// import { RenderObject, MeshRenderObject, createMeshRenderObject } from 'mol-gl/render-object';
// import { createEmptyTheme } from 'mol-theme/theme';
// import { LocationIterator } from 'mol-geo/util/location-iterator';
// import { Task } from 'mol-task';
// import { Geometry } from 'mol-geo/geometry/geometry';
// import { createIdentityTransform } from 'mol-geo/geometry/transform-data';
// import { PickingId } from 'mol-geo/geometry/picking';
// import { MarkerAction } from 'mol-geo/geometry/marker-data';
// import { Lines } from 'mol-geo/geometry/lines/lines';
// import { Text } from 'mol-geo/geometry/text/text';
// import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
// import { LinesBuilder } from 'mol-geo/geometry/lines/lines-builder';
// import { TextBuilder } from 'mol-geo/geometry/text/text-builder';



// export interface Measurements {
//     readonly inputs: Measurements.Inputs
// }

// export namespace Measurements {
//     export interface Loci {
//         readonly kind: 'measurements-loci',
//         readonly data: any,
//         readonly tag: string
//         readonly indices: OrderedSet<number>
//     }

//     export interface DistanceInput {
//         kind: 'distance-input'
//         lociA: _Loci,
//         lociB: _Loci
//     }
//     export function createDistanceInput(lociA: _Loci, lociB: _Loci): DistanceInput {
//         return { kind: 'distance-input', lociA, lociB }
//     }

//     export interface AngleInput {
//         kind: 'angle-input'
//         lociA: _Loci,
//         lociB: _Loci,
//         lociC: _Loci
//     }

//     export interface DihedralInput {
//         kind: 'dihedral-input'
//         lociA: _Loci,
//         lociB: _Loci,
//         lociC: _Loci,
//         lociD: _Loci
//     }

//     // export type MeasurementInput = [Vec3, Vec3] | [Vec3, Vec3, Vec3] | [Vec3, Vec3, Vec3, Vec3]
//     export type Input = DistanceInput | AngleInput | DihedralInput
//     export type Inputs = Input[]
// }

// // const measureVecA = Vec3.zero()
// // const measureVecB = Vec3.zero()
// // export function measure(input: MeasurementInput): number | null {
// //     switch (input.length) {
// //         case 2: return Vec3.distance(input[0], input[1])
// //         case 3: return Vec3.angle(Vec3.sub(measureVecA, input[1], input[0]), Vec3.sub(measureVecB, input[2], input[1]))
// //         default:
// //             // TODO
// //             return null
// //     }
// // }

// //

// export interface MeasurementsMeshProps {
//     detail: number,
//     sizeFactor: number
// }

// interface MeasurementsGeometries { mesh?: Mesh, lines?: Lines, text?: Text }

// export function createMeasurementsGeometries(inputs: Measurements.Inputs, props: MeasurementsMeshProps, geometries: MeasurementsGeometries): MeasurementsGeometries {
//     const { mesh, lines, text } = geometries
//     const { detail, sizeFactor } = props

//     const vertexCount = inputs.length * 10
//     const meshBuilderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh)
//     const linesBuilder = LinesBuilder.create(vertexCount, vertexCount / 2, lines)
//     const textBuilder = TextBuilder.create(props, vertexCount, vertexCount / 2, text)

//     const vA = Vec3.zero()
//     const vB = Vec3.zero()
//     const vC = Vec3.zero()
//     const vD = Vec3.zero()

//     for (let i = 0, il = inputs.length; i < il; i++) {
//         const input = inputs[i]
//         if (input.lociA) Loci.getCenter(input.lociA, vA)

//         meshBuilderState.currentGroup = i
//         linesBuilder.add()
//     }

//     return {
//         mesh: MeshBuilder.getMesh(meshBuilderState),
//         lines: linesBuilder.getLines(),
//         text: textBuilder.getText()
//     }
// }

// //

// export interface MeasurementsRepresentation<P extends MeasurementsParams> extends Representation<MeasurementInputs, P> { }

// export const MeasurementsParams = {
//     ...Mesh.Params,
// }
// export type MeasurementsParams = typeof MeasurementsParams

// export function MeasuresRepresentation<P extends MeasurementsParams>(ctx: RepresentationContext): MeasurementsRepresentation<P> {
//     let version = 0
//     const updated = new Subject<number>()
//     const _state = Representation.createState()
//     const renderObjects: RenderObject[] = []
//     let _renderObject: MeshRenderObject | undefined
//     let _inputs: Measurements.Inputs

//     let _mesh: Mesh
//     let _lines: Lines
//     let _text: Text

//     let _theme = createEmptyTheme()
//     let currentProps: PD.Values<P> = PD.getDefaultValues(MeasurementsParams) as PD.Values<P>
//     let currentParams: P
//     let locationIt: LocationIterator

//     function createOrUpdate(props: Partial<PD.Values<P>> = {}, inputs?: MeasurementInputs) {
//         currentProps = Object.assign({}, currentProps, props)
//         if (inputs) _inputs = inputs

//         return Task.create('MeasurementRepresentation.create', async runtime => {
//             renderObjects.length = 0

//             if (!_inputs) return

//             const mesh = _shape.mesh
//             locationIt = ShapeGroupIterator.fromShape(_shape)
//             const transform = createIdentityTransform()

//             const values = Mesh.createValues(mesh, transform, locationIt, _theme, currentProps)
//             const state = Geometry.createRenderableState(currentProps)

//             _renderObject = createMeshRenderObject(values, state)
//             renderObjects.push(_renderObject)
//             updated.next(version++)
//         });
//     }

//     return {
//         label: 'Measures mesh',
//         get groupCount () { return locationIt ? locationIt.count : 0 },
//         get renderObjects () { return renderObjects },
//         get props () { return currentProps },
//         get params () { return currentParams },
//         get state() { return _state },
//         get theme() { return _theme },
//         updated,
//         createOrUpdate,
//         getLoci(pickingId: PickingId) { return EmptyLoci },
//         mark(loci: Loci, action: MarkerAction) { return false },
//         setState(state: Partial<Representation.State>) {
//             if (state.visible !== undefined) renderObjects.forEach(ro => ro.state.visible = state.visible!)

//             Representation.updateState(_state, state)
//         },
//         setTheme(theme: Theme) {
//             _theme = theme
//         },
//         destroy() {
//             // TODO
//             renderObjects.length = 0
//             _renderObject = undefined
//         }
//     }
// }

// export namespace ShapeGroupIterator {
//     export function fromShape(shape: Shape): LocationIterator {
//         const { groupCount } = shape
//         const instanceCount = 1
//         const location = Shape.Location(shape)
//         const getLocation = (groupIndex: number) => {
//             location.group = groupIndex
//             return location
//         }
//         return LocationIterator(groupCount, instanceCount, getLocation)
//     }
// }