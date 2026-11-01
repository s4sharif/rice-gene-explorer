# Free private review deployment

RiceLENS is prepared for a **Render Free Python web service**, using the private GitHub repository and the included `render.yaml`. This is a reviewer prototype, not an always-available production service. The first review service is now live at [ricelens-review.onrender.com](https://ricelens-review.onrender.com) on the Free plan. The hosted health check returns ready; unauthenticated page, asset, search, expression and CSV requests are blocked. Reviewer credentials are configured in Render and are not stored in this repository.

## What is included

- Gunicorn serves the existing read-only routes with one worker and two threads.
- The hosted entry point requires a reviewer username and a password of at least 16 characters. Missing credentials fail closed. The page, JavaScript, CSS, JSON APIs and CSV exports use the same gate; only `/healthz` returns a minimal readiness response without login.
- Render provides the HTTPS address. Use that HTTPS address for reviewer access; HTTP Basic credentials must not be sent over a plaintext public connection.
- All five serving databases are bundled as approximately **73 MB compressed**, restoring to approximately **386 MB**. They include the GWAS explorer and four expression dataset versions representing three papers. This does not add Jing or Liu expression measurements.
- The expression snapshots omit only the offline `sample_summary` tables. The website uses the retained feature, cell-group, pooled-summary and metadata tables. Original research databases and raw matrices remain untouched on the Mac.
- The service restores verified snapshots into `/tmp/ricelens-data` at each startup. This supports an ephemeral filesystem without a paid disk, external database, data-download token or atlas preprocessing on the host.
- Automatic deployments are disabled; later source pushes do not immediately publish new versions.

The compressed snapshots contain public-source research data, not email correspondence, uploaded PDFs, original XLSX files or credentials. They are stored in the existing private repository for reproducible review deployment. Upstream reuse conditions still apply; this does not make the dataset public or authorize a public data release.

## Create the free service

1. Open [Render](https://dashboard.render.com/) and sign in or create an account. Keep the workspace on its free plan.
2. Choose **New → Blueprint** and connect GitHub. Grant access to `s4sharif/rice-gene-explorer`; broad access to other repositories is unnecessary.
3. Select the repository and branch `main`. Render reads `render.yaml`, which explicitly sets `plan: free` and creates no paid disk or separate database.
4. Set the prompted secrets:
   - `RICELENS_USERNAME`: a reviewer login name, for example `reviewer`.
   - `RICELENS_PASSWORD`: a new password of at least 16 characters, saved in your password manager. Do not use a university or GitHub password, commit it to Git, or paste it into the chat.
5. Review the plan and resource summary: one Free web service, no paid resources. Deploy the Blueprint.
6. After deployment succeeds, open the provided `https://...onrender.com` address. Confirm it requests reviewer authentication, then test gene search, GWAS filtering, each expression dataset, neighboring-gene comparison and CSV download.
7. Open a private/incognito window: both the homepage and an API/export URL must request authentication. Share the address and reviewer credentials with Lisa yourself; no message has been sent.

If you configure the service manually instead of using a Blueprint:

| Setting             | Value                                                                                            |
| ------------------- | ------------------------------------------------------------------------------------------------ |
| Runtime             | Python                                                                                           |
| Instance plan       | Free                                                                                             |
| Build command       | `pip install -r requirements-hosting.txt`                                                        |
| Start command       | `python3 scripts/restore_hosting_data.py && gunicorn wsgi:application --config gunicorn.conf.py` |
| Health check        | `/healthz`                                                                                       |
| `RICELENS_DATA_DIR` | `/tmp/ricelens-data`                                                                             |
| Reviewer secrets    | Same two variables above                                                                         |

## Free-tier behavior

Render Free sleeps after 15 minutes without inbound traffic; waking takes about one minute. It has 512 MB RAM and 750 free instance-hours per workspace per month. Build and outbound bandwidth limits also apply. Free instances do not provide persistent disks and are intended for prototypes. This deployment restores its fixed read-only dataset after restarts; it does not depend on retaining runtime writes. Check the dashboard's free usage limits before broad sharing. Do not upgrade or add paid resources to resolve a limit under the current free-only budget.

Official references: [Free service limits](https://render.com/docs/free), [Blueprint configuration](https://render.com/docs/blueprint-spec), [Gunicorn settings](https://gunicorn.org/reference/settings/).

## Updating or reverting data

After completing and validating a local dataset update:

```sh
python3 scripts/prepare_hosting_data.py
```

Commit the regenerated `deployment/data` snapshots and manifest separately from scientific source changes when practical. The generator checks retained row counts and SQLite integrity. Snapshots include source and metadata checksums. Restoration rejects corrupted data or a changed study manifest until the snapshot is regenerated. A reverted study whose manifest is absent stays hidden even if its old database remains in the hosting snapshot.

Push the commits, then use **Manual Deploy → Deploy latest commit** in Render. Roll back the hosting snapshot and corresponding source version together when reverting a data release. Keep source preprocessing and independent MEX validation on the research machine; these compact snapshots cannot replace the offline validation inputs.

## Local checks completed

The production WSGI interface is tested for authentication on HTML/assets/API/CSV routes, missing-secret rejection, malformed credentials, read-only methods and HEAD behavior. Snapshot tests cover corruption, path safety, metadata mismatch and restart restoration. The actual Gunicorn server was exercised against restored snapshots, including all four expression dataset versions and CSV provenance. The remote build, HTTPS readiness and unauthenticated access checks have now passed. Logs confirm restoration of all five snapshots and a successful live deployment. An authenticated browser walkthrough and a future idle-wake/restart check remain to be completed by a reviewer using the configured credentials.
