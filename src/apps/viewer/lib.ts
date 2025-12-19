import * as Structure from '../../mol-model/structure';
import { DataLoci, EveryLoci, Loci } from '../../mol-model/loci';
import { Volume } from '../../mol-model/volume';
import { Shape, ShapeGroup } from '../../mol-model/shape';
import * as LinearAlgebra3D from '../../mol-math/linear-algebra/3d';

export const lib = {
    structure: {
        ...Structure,
    },
    volume: {
        Volume,
    },
    shape: {
        Shape,
        ShapeGroup,
    },
    loci: {
        Loci,
        DataLoci,
        EveryLoci,
    },
    math: {
        LinearAlgebra: {
            ...LinearAlgebra3D,
        }
    }
};