#!/usr/bin/env node
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import cluster from 'cluster';
import { runChild } from '@molstar/model-server/preprocess/parallel';

if (cluster.isPrimary) {
    await import('@molstar/model-server/preprocess/master');
} else {
    runChild();
}

// example:
// node build\node_modules\servers\model\preprocess -i e:\test\Quick\1cbs_updated.cif -oc e:\test\mol-star\model\1cbs.cif -ob e:\test\mol-star\model\1cbs.bcif
