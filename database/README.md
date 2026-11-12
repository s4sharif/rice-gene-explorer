# Database definitions

`schema.sql` defines the gene and GWAS database. Build it from the repository root with `python3 -m pipelines.build_database`. Generated SQLite files live in `data/databases/`; source inputs live in `data/raw/`. Expression tables are created by `pipelines/build_expression.py`.
