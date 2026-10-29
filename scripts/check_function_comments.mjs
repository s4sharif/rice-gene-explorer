import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { resolve } from "node:path";
import { parse } from "acorn";

const root = fileURLToPath(new URL("../", import.meta.url));
const functionTypes = new Set([
  "FunctionDeclaration",
  "FunctionExpression",
  "ArrowFunctionExpression",
]);

// Visit every syntax node and verify a sentence comment immediately before each function.
function inspectNode(node, comments, source, filename, failures) {
  let count = 0;
  if (functionTypes.has(node.type)) {
    count++;
    let precedingComment;
    for (const comment of comments) {
      if (comment.end <= node.start) precedingComment = comment;
    }
    const hasSentence = precedingComment?.value.trim().endsWith(".");
    const isAdjacent =
      precedingComment && source.slice(precedingComment.end, node.start).trim() === "";
    if (!hasSentence || !isAdjacent) {
      failures.push(
        `${filename}:${node.loc.start.line}: function needs a sentence comment`,
      );
    }
  }
  for (const value of Object.values(node)) {
    const children = Array.isArray(value) ? value : [value];
    for (const child of children) {
      if (child && typeof child.type === "string") {
        count += inspectNode(child, comments, source, filename, failures);
      }
    }
  }
  return count;
}

// Parse project JavaScript files and report any missing function comments.
function main() {
  const failures = [];
  let count = 0;
  for (const filename of ["static/app.js", "scripts/check_function_comments.mjs"]) {
    const source = readFileSync(resolve(root, filename), "utf8");
    const comments = [];
    const tree = parse(source, {
      ecmaVersion: "latest",
      sourceType: "module",
      locations: true,
      onComment: comments,
    });
    count += inspectNode(tree, comments, source, filename, failures);
  }
  if (failures.length) throw new Error(failures.join("\n"));
  console.log(`Verified comments above ${count} JavaScript functions.`);
}

main();
