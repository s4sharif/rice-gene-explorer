# Working on the project

Use Python 3.11 or newer for the application. Node.js and the tools below are needed only when editing or checking source code.

```sh
python3 -m pip install -r requirements-dev.txt
npm ci
```

Keep a one-sentence explanatory comment immediately above every function, including test methods, nested functions and JavaScript event handlers. Use named functions so the behavior is easy to follow. Python uses four-space indentation; JavaScript, HTML and CSS use two spaces.

Format and check changes before committing:

```sh
black .
ruff check .
npm run format
black --check .
npm run format:check
python3 scripts/check_function_comments.py
npm run check:comments
npm run check:syntax
python3 -m unittest discover -s tests -v
```

The integration tests require a built database. Follow the README to download source snapshots and build it first. When changing imports or statistical interpretation, also run `python3 validate_data.py` and review the relevant documentation.

Preserve original source fields, document genome compatibility and retain missing values explicitly. Do not invent allele information, interpret a display window as an LD region or label nearby genes as causal. Keep downloaded data, database files, credentials and personal correspondence out of commits. Source release labels and checksums are necessary scientific provenance.

GitHub Actions checks formatting, function comments, source checksums, database imports and integration tests. Its data steps require the upstream download servers to be available and to retain the recorded snapshots; failures from a changed source should be investigated before updating checksums.
