# gene-explorer-js

JavaScript/Vite port of [`../gene-explorer`](../gene-explorer) (the R/Shiny gene-lookup app),
built so it can be deployed as a static site on Vercel instead of requiring an R server.

Plain HTML/CSS/JS — no UI framework, no runtime dependencies. All gene stats are precomputed
and shipped as a static JSON file (`public/data/gene_stats.json`); there is no backend.

## Regenerating the data

If `../gene-explorer/gene_stats.rds` changes (re-run `../compute_gene_stats.R`), regenerate
the JSON the app reads:

```
Rscript scripts/export_gene_stats.R
```

This requires R with the `jsonlite` package installed. It writes
`public/data/gene_stats.json` (~3-4 MB).

## Develop

```
npm install
npm run dev
```

## Build / deploy

```
npm run build     # outputs static site to dist/
npm run preview   # sanity-check the production build locally
```

On Vercel: create a project pointed at this repo, set **Root Directory** to
`hypomap-viz/gene-explorer-js`. Framework preset "Vite" is auto-detected
(build command `vite build`, output directory `dist`) — no other config needed.

## Status

Currently ports the core lookup flow: gene search/suggestions, global summary stats, and
the by-cell-type table + bar charts. Not yet ported from the R app: the "About the dataset"
details panel and citation links.
