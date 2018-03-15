/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Database from './db/database'
import Table from './db/table'
import Column from './db/column'
import * as ColumnHelpers from './db/column-helpers'

type DatabaseCollection = { [name: string]: Database<Database.Schema> }

export { DatabaseCollection, Database, Table, Column, ColumnHelpers }