// Generated in 2019-01-30T16:38:09-08:00
export type Maybe<T> = T | null;

/** Built-in scalar representing an instant in time */
export type Date = any;

/** Unrepresentable type */
export type Unrepresentable = any;

// ====================================================
// Documents
// ====================================================

export namespace AssemblySymmetry {
  export type Variables = {
    pdbId: string;
  };

  export type Query = {
    __typename?: "Query";

    assemblies: Maybe<(Maybe<Assemblies>)[]>;
  };

  export type Assemblies = {
    __typename?: "CoreAssembly";

    pdbx_struct_assembly: Maybe<PdbxStructAssembly>;

    rcsb_struct_symmetry: Maybe<(Maybe<RcsbStructSymmetry>)[]>;
  };

  export type PdbxStructAssembly = {
    __typename?: "PdbxStructAssembly";

    id: string;
  };

  export type RcsbStructSymmetry = {
    __typename?: "RcsbStructSymmetry";

    clusters: (Maybe<Clusters>)[];

    kind: string;

    oligomeric_state: string;

    rotation_axes: Maybe<(Maybe<RotationAxes>)[]>;

    stoichiometry: (Maybe<string>)[];

    symbol: string;

    type: string;
  };

  export type Clusters = {
    __typename?: "RcsbStructSymmetryClusters";

    avg_rmsd: Maybe<number>;

    members: (Maybe<Members>)[];
  };

  export type Members = {
    __typename?: "ClustersMembers";

    asym_id: string;

    pdbx_struct_oper_list_ids: Maybe<(Maybe<string>)[]>;
  };

  export type RotationAxes = {
    __typename?: "RcsbStructSymmetryRotationAxes";

    start: (Maybe<number>)[];

    end: (Maybe<number>)[];

    order: Maybe<number>;
  };
}
