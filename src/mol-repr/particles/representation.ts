/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParticleList } from '../../mol-model/particles/particle-list';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Representation, RepresentationProvider } from '../representation';

export interface ParticleRepresentation<P extends PD.Params = PD.Params> extends Representation<ParticleList, P> { }

export type ParticleRepresentationProvider<P extends PD.Params, Id extends string = string> = RepresentationProvider<ParticleList, P, Representation.State, Id>
export function ParticleRepresentationProvider<P extends PD.Params, Id extends string>(p: ParticleRepresentationProvider<P, Id>): ParticleRepresentationProvider<P, Id> { return p; }
