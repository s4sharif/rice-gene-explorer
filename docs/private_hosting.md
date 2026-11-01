# Private RiceLENS hosting

The selected first-review budget is **free only**. A Render Free deployment package is now prepared with bundled read-only snapshots and mandatory reviewer authentication; see [free deployment steps](free_hosting.md). The first review service is live at [ricelens-review.onrender.com](https://ricelens-review.onrender.com), with reviewer authentication.

The current localhost address works only on Sharif's Mac. A private GitHub repository does not make a deployed website private: advisor access must be enforced for the website, API and CSV endpoints.

| Option                                                           | Potential cost                                                                                                                                                                        | Suitability                                                                                                                        |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| University or lab-managed server with institutional login or VPN | Requires an internal quote; existing capacity might avoid an additional hosting bill, but no allocation is confirmed                                                                  | Preferred long-term ownership and access management if the lab can provide a server                                                |
| Railway web service with an authentication gate                  | Hobby subscription minimum US$5/month, including US$5 usage; actual compute, memory, storage and network usage may exceed this. Pro is US$20/month minimum if the lab needs that plan | Practical review deployment; budget for usage and choose an account plan appropriate to lab ownership                              |
| Render Free web service with an authentication gate              | US$0 compute within free limits; additional usage restrictions apply                                                                                                                  | Temporary review only: sleeps after 15 minutes idle, no persistent disk, filesystem changes are lost on restart/redeploy/spin-down |
| Paid Render service with authentication and persistent disk      | Obtain the current service/disk quote before choosing; no fixed total verified here                                                                                                   | Alternative for an always-available SQLite deployment                                                                              |

Pricing is in USD, before applicable taxes or currency conversion. Railway's minimum is not a fixed all-inclusive quote. The Render Free service has been deployed. No paid subscription or paid storage was added.

## Earlier paid/institutional options

Use a lab-owned university service if available; otherwise start a small cloud review deployment with a spending limit and access restricted to named reviewers. Put an institutional/OIDC access gateway in front of every route, or use an equivalent authenticated reverse proxy. Test that an unauthenticated visitor cannot retrieve HTML, search JSON, source records or CSV exports. A hard-to-guess URL and robots exclusion are not authentication.

The local standard-library server remains a local prototype. A separate authenticated Gunicorn entry point, PORT configuration, external database directory and reproducible compact snapshots are now prepared for Render Free. The initial deployment is live; health checks and unauthenticated access protection have been verified. The complete serving snapshot restores to approximately 386 MB, including current expression data. Keep raw single-cell matrices in offline processing storage; serve compact gene-by-cell-type summaries rather than loading an atlas into each web request. Verify restart behavior, authenticated search and exports before sharing the address with Lisa.

## Official references

- [Railway plans](https://docs.railway.com/pricing/plans)
- [Railway billing](https://docs.railway.com/pricing/understanding-your-bill)
- [Railway volumes](https://docs.railway.com/volumes/reference)
- [Render Free limitations](https://render.com/docs/free)
- [Render pricing](https://render.com/pricing)
