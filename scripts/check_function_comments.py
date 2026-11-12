"""Check that every Python function has a short explanatory comment above it."""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


# Reject Python functions whose definition lacks an immediately preceding sentence comment.
def main():
    failures = []
    count = 0
    paths = sorted((ROOT / "backend").glob("*.py"))
    paths += sorted((ROOT / "pipelines").glob("*.py"))
    paths += sorted((ROOT / "deployment").glob("*.py"))
    paths += sorted((ROOT / "scripts").glob("*.py"))
    paths += sorted((ROOT / "tests").glob("*.py"))
    for path in paths:
        source = path.read_text()
        lines = source.splitlines()
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.Lambda):
                failures.append(
                    f"{path.relative_to(ROOT)}:{node.lineno}: use a documented named function"
                )
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            count += 1
            first_line = min(
                [node.lineno] + [item.lineno for item in node.decorator_list]
            )
            comment = lines[first_line - 2].strip() if first_line > 1 else ""
            if not comment.startswith("# ") or not comment.endswith("."):
                failures.append(
                    f"{path.relative_to(ROOT)}:{node.lineno}: {node.name} needs a sentence comment"
                )
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"Verified comments above {count} Python functions.")


if __name__ == "__main__":
    main()
