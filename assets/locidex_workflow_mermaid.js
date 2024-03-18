flowchart LR
    subgraph Merge [Merging results for genetic distance calculations]
    direction LR
    P1[/"Profile\n(JSON)"/] --> M[Merge]
    M -.-> M1[/"Merged profile\n(TSV)"/]
    M1 --> G[GrapeTree]
    M1 --> cg[cgMLST dists]
    M1 --> pd[profile dists]
    end

    subgraph Search [Searching a sequence database for gene-by-gene analysis]
    direction LR
    db1[(Locidex DB)] --> S[Search]
    Seq1[/"Sequences\n(Genbank | fasta)"/] --> S
    S -.-> SS[/"seq_store\n(JSON)"/]
    SS --> R[Report]
    R -.-> P[/"Profile\n(JSON)"/]
    end

    subgraph Extract [Extracting loci]
    direction LR
    db2[(Locidex DB)] --> E[Extract]
    Seq2[/"Sequences\n(Genbank | fasta)"/] --> E
    E -.-> E1[/"Extracted\nseqs (fasta)"/]
    end

    subgraph Build [Building a sequence database]
    direction LR
    Seq[/"Sequences\n(fasta)"/] --> For[Format]
    For --> B[Build]
    B --> db[(Locidex DB)]
    end 