# Reverting updates

Updates are grouped into focused commits:

| Commit    | Update                                                    |
| --------- | --------------------------------------------------------- |
| `9971d3d` | Gene search, GWAS evidence labels and RiceLENS branding   |
| `bfbcaea` | Shared expression framework and Robertson control dataset |
| `3c66bc9` | Wang normal-root dataset                                  |
| `39e9a9e` | Li curated diploid pistil dataset                         |

To remove one paper addition while keeping the history, run:

```sh
git revert 3c66bc9
git push origin main
```

Replace the commit ID with `39e9a9e` to remove the Li addition. Restart the Python server after changing code or study manifests, then refresh the browser.

Each paper's tracked manifest enables its dataset. Reverting the paper commit removes that manifest, so its data are hidden even if an ignored local SQLite file remains. Raw downloads and databases are not committed; reverting does not delete those local files.

The shared expression framework is a prerequisite for later paper additions. Remove dependent paper additions before reverting the framework. Older website changes may also overlap later changes, so Git may require resolving conflicts when reverting those commits.
