/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetry } from '../../../../mol-model-props/rcsb/assembly-symmetry';
import { AttachModelProperty } from '../../property-provider';
import { ajaxGet } from '../../../../mol-util/data-source';
import { GraphQLClient } from '../../../../mol-util/graphql-client';
import { getParam } from '../../../common/util';

export const RCSB_assemblySymmetry: AttachModelProperty = ({ model, params }) => {
    const url = getApiUrl(params, 'assembly_symmetry', `https:${AssemblySymmetry.GraphQLEndpointURL}`)
    const client = new GraphQLClient(url, ajaxGet)
    return AssemblySymmetry.attachFromCifOrAPI(model, client)
}

function getApiUrl(params: any, name: string, fallback: string) {
    const path = getParam<string>(params, 'RCSB', 'API', name);
    if (!path) return fallback;
    return path;
}