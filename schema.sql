PRAGMA foreign_keys=ON;

CREATE TABLE assembly (
    id TEXT PRIMARY KEY,
    species TEXT NOT NULL,
    cultivar TEXT NOT NULL,
    coordinate_system TEXT NOT NULL
);

CREATE TABLE source (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    url TEXT NOT NULL,
    release TEXT NOT NULL,
    sha256 TEXT NOT NULL
);

CREATE TABLE gene (
    id TEXT PRIMARY KEY,
    assembly_id TEXT NOT NULL REFERENCES assembly(id),
    annotation_id TEXT NOT NULL REFERENCES source(id),
    chromosome INTEGER NOT NULL CHECK(chromosome BETWEEN 1 AND 12),
    start INTEGER NOT NULL CHECK(start>0),
    end INTEGER NOT NULL CHECK(end>=start),
    strand TEXT NOT NULL CHECK(strand IN ('+','-')),
    symbol TEXT,
    description TEXT
);

CREATE INDEX gene_region ON gene(assembly_id,chromosome,start,end);

CREATE TABLE alias (
    gene_id TEXT REFERENCES gene(id),
    name TEXT COLLATE NOCASE NOT NULL,
    kind TEXT NOT NULL,
    source_id TEXT NOT NULL REFERENCES source(id),
    PRIMARY KEY(gene_id,name,kind)
);

CREATE INDEX alias_search ON alias(name COLLATE NOCASE);

CREATE TABLE publication (
    pmid TEXT PRIMARY KEY,
    title TEXT,
    journal TEXT,
    author TEXT,
    published TEXT,
    url TEXT NOT NULL
);

CREATE TABLE study (
    id TEXT PRIMARY KEY,
    source_study_id TEXT NOT NULL,
    publication_id TEXT REFERENCES publication(pmid),
    model TEXT,
    population TEXT,
    sample_size TEXT,
    environment TEXT,
    tissue TEXT
);

CREATE TABLE trait (
    id TEXT PRIMARY KEY,
    label TEXT NOT NULL,
    ontology_id TEXT
);

CREATE TABLE association (
    id TEXT PRIMARY KEY,
    source_id TEXT NOT NULL REFERENCES source(id),
    source_row INTEGER NOT NULL,
    study_id TEXT NOT NULL REFERENCES study(id),
    trait_id TEXT NOT NULL REFERENCES trait(id),
    assembly_id TEXT NOT NULL REFERENCES assembly(id),
    chromosome INTEGER NOT NULL CHECK(chromosome BETWEEN 1 AND 12),
    position INTEGER NOT NULL CHECK(position>0),
    p_value REAL NOT NULL CHECK(p_value>0 AND p_value<=1),
    p_value_text TEXT NOT NULL,
    effect_allele TEXT,
    effect_allele_frequency TEXT,
    effect_size REAL,
    other_allele TEXT,
    r2_percent TEXT,
    source_reference TEXT NOT NULL,
    reported_location TEXT,
    reported_gene TEXT,
    raw_json TEXT NOT NULL
);

CREATE INDEX association_region ON association(assembly_id,chromosome,position);

CREATE INDEX association_trait ON association(trait_id);

CREATE TABLE quarantine (
    source_id TEXT REFERENCES source(id),
    source_row INTEGER,
    reason TEXT NOT NULL,
    raw_json TEXT NOT NULL,
    PRIMARY KEY(source_id,source_row)
);
